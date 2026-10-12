module nlevel_cooling_mod
!---------------------------------------------------------------------------
! MoCHII: line cooling from level populations at the local electron density.
!
! For an ion the level populations normalized to the ion density depend only
! on (T_e, n_e) (electron collisions + spontaneous decay, optically thin, no
! photo-excitation), so the collisional line cooling COEFFICIENT
!   Lambda_SE(T, ne) = (sum_k emis_k(T, ne)) / ne   [erg cm^3 s^-1, over the
!                       NRAD radiative transitions; the volume cooling is
!                       ne n_ion Lambda_SE, matching cooling_eval]
! is a function of (T, n_e) alone and is tabulated at setup from the n-level
! atom files (data/atomic/nlevel_<el>_<stage>.txt). Store absolute Lambda_SE
! values, including positive subnormal values, at every grid point. The
! n_e = 10^-6 cm^-3 column is the numerical low-density reference. Cooling
! need not decrease monotonically with density for every atom and temperature.
!
! Runtime cooling combines the existing low-density fit lam_fit with the
! registered atom. A source-converged atom uses the density ratio
! Lambda = lam_fit * Lambda_SE(T,ne)/Lambda_SE(T,ne_floor) when the reference
! is positive; an exact zero reference uses absolute Lambda_SE(T,ne).
! A partial atom retains the positive low-density residual when lam_fit is
! at least the registered reference, otherwise it scales the registered
! density dependence by lam_fit/reference. The coverage flag requires
! independent source convergence before it is enabled.
! Ions with no n-level model get idx = 0 -> lam_fit returned unchanged.
!
! The table is built here from the SAME solver (nlevel_mod) that produces the
! emergent line diagnostics, so the thermal balance and the line list are one
! physics. Interpolation uses log Lambda linear in 1/T for positive endpoints
! and log Lambda linear in log n_e for positive column values. Zero endpoints
! use linear interpolation. Density and temperatures below the grid are
! clamped. Temperature-table atoms end at their last source temperature.
!---------------------------------------------------------------------------
  use define
  use nlevel_mod, only : nlevel_atom_type, nlevel_load, &
                         nlevel_cooling_prepT, nlevel_cooling_solve
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  public :: nlevel_cooling_add, nlevel_cooling_apply
  public :: nlevel_cooling_lam_se, nlevel_cooling_report
  public :: reg,col_cooling

  !--- (log T, log n_e) grid: log10 T = 2..6 step 0.01 (401 points); log10 n_e
  !--- = -6..10 step 0.1 (161 points).  ne_floor = 10^-6 is the reference
  !--- (low-density) column.  This is finer in T than the initial D4 spec
  !--- (0.05/0.2 dex): the gate G-a criterion (worst 0.5% vs the direct solve)
  !--- is set by the low-T/high-n_e corner where a high-Tx line and a saturated
  !--- low-n_crit line become comparable. Runtime lookup cost is independent
  !--- of the grid size. Full-registry build time depends on atom size.
  integer, parameter :: NT = 401, NNE = 161
  !--- H I plus every stage of the supported metal elements.  Only ions
  !--- with an n-level file allocate their (T, n_e) tables.
  integer, parameter :: MAXION = 1 + MAX_METAL_ELEMENTS*MAX_ION_STAGES
  real(kind=wp), parameter :: T_LO  =  2.0_wp, DT_L  = 0.01_wp
  real(kind=wp), parameter :: NE_LO = -6.0_wp, DNE_L = 0.1_wp

  !--- Absolute cooling arrays are allocated only for ions actually tabulated.
  type cool_table_type
     character(len=16) :: name = ''
     integer :: stage = 0
     integer :: nlev  = 0
     logical :: density_ratio_cooling = .false.
     real(kind=wp) :: source_ceiling_ln_K = huge(1.0_wp)
     real(kind=wp), allocatable :: logT(:)
     real(kind=wp), allocatable :: lambda(:,:)
  end type cool_table_type

  logical :: inited = .false.
  real(kind=wp) :: gT(NT), gne(NNE)
  integer :: nion_reg = 0
  type(cool_table_type) :: reg(MAXION)
  real(kind=wp) :: build_time_total = 0.0_wp

contains

  !=========================================================================
  ! Fill the (log T, log n_e) grid once (lazy, on the first add).
  !=========================================================================
  subroutine nlevel_cooling_init()
    implicit none
    integer :: i, j
    if (inited) return
    do i = 1, NT
       gT(i)  = T_LO  + DT_L *real(i-1, wp)
    end do
    do j = 1, NNE
       gne(j) = NE_LO + DNE_L*real(j-1, wp)
    end do
    inited = .true.
  end subroutine nlevel_cooling_init

  !=========================================================================
  ! Register an ion's cooling table.  Loads data/atomic/nlevel_<name>_<stage>
  ! (par%atomic_dir, '_full' suffix for Fe II/III under par%fe_levels_full, exactly
  ! as lines_mod), builds the absolute cooling table, and returns its index.
  ! Returns 0 (no table; caller keeps lam_fit) when the file does not load.
  !=========================================================================
  integer function nlevel_cooling_add(name, stage) result(idx)
    implicit none
    character(len=*), intent(in) :: name
    integer,          intent(in) :: stage
    type(nlevel_atom_type) :: atom
    character(len=PATH_LEN) :: fname
    character(len=8)   :: nlsuf
    logical :: ok, limited_grid
    real(kind=wp), allocatable :: Rcoll(:,:), Rrad(:,:), lam(:,:)
    real(kind=wp) :: Tk, nek, final_T, last_logT
    integer(kind=i8b) :: c0, c1, crate
    integer :: i, j, k, ie, nt_ion

    idx = 0
    call nlevel_cooling_init()

    !--- Expanded n-level files are supplied only for Fe II and Fe III.
    nlsuf = ''
    if (par%fe_levels_full .and. trim(name) == 'fe' .and. &
        (stage == 2 .or. stage == 3)) nlsuf = '_full'
    write(fname,'(a,a,a,i0,a,a)') trim(par%atomic_dir)//'/nlevel_', &
       trim(name), '_', stage, trim(nlsuf), '.txt'
    call nlevel_load(trim(fname), atom, ok)
    if (.not. ok) return

    if (nion_reg >= MAXION) then
       if (mpar%p_rank == 0) write(*,'(a,i0,a)') &
          ' NLCOOL: registry full (MAXION = ', MAXION, '); ion left in the '// &
          'low-density limit.'
       return
    end if

    ie = nion_reg + 1
    ! A ttype=0 collision row is the version-2 TEMPERATURE_TABLE model.
    ! Its final node stores ln(T/K). Every collision row must support a
    ! temperature used to build this ion's cooling table.
    do k = 1, atom%nups
       if (atom%ttype(k) == 0) then
          reg(ie)%source_ceiling_ln_K = min(reg(ie)%source_ceiling_ln_K, &
               atom%node_x(atom%offset(k)+atom%nc(k)-1))
       end if
    end do
    limited_grid = reg(ie)%source_ceiling_ln_K < log(10.0_wp)*gT(NT)
    if (limited_grid) then
       final_T = exp(reg(ie)%source_ceiling_ln_K)
       ! exp(ln T) can round one representable value above the last source
       ! node; keep the actual final solve within the source interval.
       if (log(final_T) > reg(ie)%source_ceiling_ln_K) &
          final_T = nearest(final_T, -1.0_wp)
       last_logT = log10(final_T)
       nt_ion = count(gT < last_logT - 32.0_wp*epsilon(last_logT)) + 1
       if (nt_ion < 2) error stop 'NLCOOL: source ceiling is below the second cooling temperature'
       allocate(reg(ie)%logT(nt_ion))
       reg(ie)%logT(1:nt_ion-1) = gT(1:nt_ion-1)
       reg(ie)%logT(nt_ion) = last_logT
    else
       nt_ion = NT
       allocate(reg(ie)%logT(NT))
       reg(ie)%logT = gT
    end if

    allocate(Rcoll(atom%nlev,atom%nlev), Rrad(atom%nlev,atom%nlev), lam(nt_ion,NNE))
    call system_clock(c0, crate)
    do i = 1, nt_ion
       Tk = 10.0_wp**reg(ie)%logT(i)
       if (i == nt_ion .and. limited_grid) Tk = final_T
       call nlevel_cooling_prepT(atom, Tk, Rcoll, Rrad)
       do j = 1, NNE
          nek = 10.0_wp**gne(j)
          lam(i,j) = nlevel_cooling_solve(atom, Rcoll, Rrad, nek)
       end do
    end do
    call system_clock(c1)
    build_time_total = build_time_total + real(c1-c0, wp)/real(crate, wp)

    nion_reg = ie
    reg(ie)%name  = name
    reg(ie)%stage = stage
    reg(ie)%nlev  = atom%nlev
    reg(ie)%density_ratio_cooling = atom%density_ratio_cooling
    allocate(reg(ie)%lambda(nt_ion, NNE))
    reg(ie)%lambda(:,:) = lam(:,:)
    idx = ie

    if (mpar%p_rank == 0) write(*,'(a,a,a,i0,a,i0,a,f8.4)') &
       ' NLCOOL: ', trim(name), ' stage ', stage, &
       ' tabulated (nlev = ', atom%nlev, &
       '), Lambda/Lambda_ref(1e4 K, 1e4 cm^-3) = ', suppression_ratio(ie, 1.0e4_wp, 1.0e4_wp)
  end function nlevel_cooling_add

  !=========================================================================
  ! Runtime cooling coefficient (the D2 combiner).  idx <= 0 -> lam_fit.
  !=========================================================================
  real(kind=wp) function nlevel_cooling_apply(idx, lam_fit, T, ne) result(lam)
    implicit none
    integer,       intent(in) :: idx
    real(kind=wp), intent(in) :: lam_fit, T, ne
    real(kind=wp) :: lamref, lam_se
    if (idx <= 0) then
       lam = lam_fit
       return
    end if
    lamref = low_density_cooling(idx, T)
    lam_se = cooling_coeff(idx, T, ne)
    ! A complete or source-converged atom supplies the density dependence
    ! of the entire fitted cooling. Fit/interpolation errors cannot be
    ! interpreted as unsuppressed emission from omitted levels. Truncated
    ! models require full-source convergence before setting DENSITY_RATIO.
    ! Current convergence checks cover 1000--100000 K and ne <= 1e10 cm^-3;
    ! table evaluations outside that range do not establish source accuracy.
    if (reg(idx)%density_ratio_cooling) then
       if (lamref > 0.0_wp) then
          lam = lam_fit*(lam_se/lamref)
       else
          lam = lam_se
       end if
       return
    end if
    if (lamref <= lam_fit) then
       lam = lam_se + (lam_fit - lamref)
    else
       lam = (lam_fit/lamref)*lam_se
    end if
  end function nlevel_cooling_apply

  !=========================================================================
  ! Reconstructed cooling coefficient Lambda_SE(T, ne) from the table (the
  ! direct-solve reference is nlevel_mod::nlevel_cooling_at).
  !=========================================================================
  real(kind=wp) function nlevel_cooling_lam_se(idx, T, ne) result(lam)
    implicit none
    integer,       intent(in) :: idx
    real(kind=wp), intent(in) :: T, ne
    if (idx <= 0) then
       lam = 0.0_wp
       return
    end if
    lam = cooling_coeff(idx, T, ne)
  end function nlevel_cooling_lam_se

  !=========================================================================
  ! Cooling coefficient Lambda_SE(T, ne) [erg cm^3 s^-1] by 2-D interpolation
  ! of the stored table.  Within a log-T cell: ln(Lambda) linear in 1/T --
  ! exact for a single line (Lambda ~ exp(-Tx/T)), so the exp-suppressed tails
  ! that linear-in-log-T cannot represent are captured, INCLUDING the high-ne
  ! corner where one high-Tx line survives after the low-n_crit lines saturate
  ! (there the leading exp does not cancel out of the suppression).  Then
  ! log-linear in log ne between the two bracketing columns.  Edge clamping:
  ! log ne below the floor -> the floor (reference) column, above the ceiling
  ! -> the last column. T is clamped below the grid; source-table ceilings
  ! reject unsupported higher temperatures.
  !=========================================================================
  real(kind=wp) function cooling_coeff(idx, T, ne) result(lam)
    implicit none
    integer,       intent(in) :: idx
    real(kind=wp), intent(in) :: T, ne
    real(kind=wp) :: lne, w1, wt, wn, L0, L1
    integer :: it, je
    call cooling_temperature_cell(idx, T, it, w1, wt)
    if (.not. ieee_is_finite(ne)) &
       error stop 'NLCOOL: electron density must be positive and finite'
    if (ne <= 0.0_wp) &
       error stop 'NLCOOL: electron density must be positive and finite'
    lne = log10(ne)
    if (lne <= gne(1)) then
       lam = col_cooling(idx, it, 1, w1, wt)
    else if (lne >= gne(NNE)) then
       lam = col_cooling(idx, it, NNE, w1, wt)
    else
       je = min(int((lne - gne(1))/DNE_L) + 1, NNE-1)
       wn = (lne - gne(je))/(gne(je+1) - gne(je))
       L0 = col_cooling(idx, it, je,   w1, wt)
       L1 = col_cooling(idx, it, je+1, w1, wt)
       if (L0 > 0.0_wp .and. L1 > 0.0_wp) then
          lam = exp(log(L0) + wn*(log(L1) - log(L0)))
       else
          lam = L0*(1.0_wp - wn) + L1*wn
       end if
    end if
  end function cooling_coeff

  !=========================================================================
  ! Cooling coefficient at (T, ne-column j): shape-preserving cubic ln(Lambda)
  ! in -1/T. A constant Boltzmann excitation energy has linear ln(Lambda)
  ! in this coordinate. Linear interpolation alone exceeded the 0.5% density-
  ! ratio criterion for supplied atoms on the 401-temperature grid. Keep
  ! nonnegative linear interpolation in log T when either endpoint is zero.
  ! Positive subnormal grid values are retained directly.
  !=========================================================================
  real(kind=wp) function col_cooling(idx, it, j, w1, wt) result(val)
    implicit none
    integer,       intent(in) :: idx, it, j
    real(kind=wp), intent(in) :: w1, wt
    real(kind=wp) :: Li, Lip, h, y0, y1, m0, m1, y, w2, w3
    Li  = reg(idx)%lambda(it,j)
    Lip = reg(idx)%lambda(it+1,j)
    if (Li > 0.0_wp .and. Lip > 0.0_wp) then
       y0 = log(Li); y1 = log(Lip)
       h = 10.0_wp**(-reg(idx)%logT(it)) - 10.0_wp**(-reg(idx)%logT(it+1))
       m0 = cooling_log_inverse_temperature_slope(idx,it,j)
       m1 = cooling_log_inverse_temperature_slope(idx,it+1,j)
       w2 = w1*w1; w3 = w2*w1
       y = (2*w3-3*w2+1)*y0 + (w3-2*w2+w1)*h*m0 &
         + (-2*w3+3*w2)*y1 + (w3-w2)*h*m1
       ! The monotone Hermite polynomial is bounded by its endpoints.
       ! Enforce that bound against floating-point cancellation only.
       y = min(max(y,min(y0,y1)),max(y0,y1))
       val = exp(y)
    else
       val = Li*(1.0_wp - wt) + Lip*wt
    end if
  end function col_cooling

  ! Derivative of ln(Lambda) with respect to -1/T at an existing table node.
  ! Weighted harmonic slopes preserve shape and a linear Boltzmann exponent.
  ! A zero neighboring coefficient restricts the slope to the positive pair;
  ! its logarithm is never evaluated. No additional temperature solves occur.
  real(kind=wp) function cooling_log_inverse_temperature_slope(idx,i,j) result(m)
    integer, intent(in) :: idx,i,j
    integer :: n
    real(kind=wp) :: h0,h1,h2,d0,d1,d2,w0,w1
    n = size(reg(idx)%logT)
    m = 0.0_wp
    if (i == 1) then
       if (reg(idx)%lambda(2,j) <= 0.0_wp) return
       h0 = 10.0_wp**(-reg(idx)%logT(1))-10.0_wp**(-reg(idx)%logT(2))
       d0 = (log(reg(idx)%lambda(2,j))-log(reg(idx)%lambda(1,j)))/h0
       m = d0
       if (n == 2) return
       if (reg(idx)%lambda(3,j) <= 0.0_wp) return
       h1 = 10.0_wp**(-reg(idx)%logT(2))-10.0_wp**(-reg(idx)%logT(3))
       d1 = (log(reg(idx)%lambda(3,j))-log(reg(idx)%lambda(2,j)))/h1
       m = cooling_log_inverse_temperature_end_slope(h0,h1,d0,d1)
    else if (i == n) then
       if (reg(idx)%lambda(n-1,j) <= 0.0_wp) return
       h0 = 10.0_wp**(-reg(idx)%logT(n-1))-10.0_wp**(-reg(idx)%logT(n))
       d0 = (log(reg(idx)%lambda(n,j))-log(reg(idx)%lambda(n-1,j)))/h0
       m = d0
       if (n == 2) return
       if (reg(idx)%lambda(n-2,j) <= 0.0_wp) return
       h1 = 10.0_wp**(-reg(idx)%logT(n-2))-10.0_wp**(-reg(idx)%logT(n-1))
       d1 = (log(reg(idx)%lambda(n-1,j))-log(reg(idx)%lambda(n-2,j)))/h1
       m = cooling_log_inverse_temperature_end_slope(h0,h1,d0,d1)
    else
       h0 = 10.0_wp**(-reg(idx)%logT(i-1))-10.0_wp**(-reg(idx)%logT(i))
       h1 = 10.0_wp**(-reg(idx)%logT(i))-10.0_wp**(-reg(idx)%logT(i+1))
       if (reg(idx)%lambda(i-1,j) <= 0.0_wp) then
          if (reg(idx)%lambda(i+1,j) <= 0.0_wp) return
          d1 = (log(reg(idx)%lambda(i+1,j))-log(reg(idx)%lambda(i,j)))/h1
          m = d1
          if (i+2 <= n) then
             if (reg(idx)%lambda(i+2,j) > 0.0_wp) then
                h2 = 10.0_wp**(-reg(idx)%logT(i+1))-10.0_wp**(-reg(idx)%logT(i+2))
                d2 = (log(reg(idx)%lambda(i+2,j))-log(reg(idx)%lambda(i+1,j)))/h2
                m = cooling_log_inverse_temperature_end_slope(h1,h2,d1,d2)
             end if
          end if
          return
       end if
       if (reg(idx)%lambda(i+1,j) <= 0.0_wp) then
          d0 = (log(reg(idx)%lambda(i,j))-log(reg(idx)%lambda(i-1,j)))/h0
          m = d0
          if (i-2 >= 1) then
             if (reg(idx)%lambda(i-2,j) > 0.0_wp) then
                h2 = 10.0_wp**(-reg(idx)%logT(i-2))-10.0_wp**(-reg(idx)%logT(i-1))
                d2 = (log(reg(idx)%lambda(i-1,j))-log(reg(idx)%lambda(i-2,j)))/h2
                m = cooling_log_inverse_temperature_end_slope(h0,h2,d0,d2)
             end if
          end if
          return
       end if
       d0 = (log(reg(idx)%lambda(i,j))-log(reg(idx)%lambda(i-1,j)))/h0
       d1 = (log(reg(idx)%lambda(i+1,j))-log(reg(idx)%lambda(i,j)))/h1
       if (d0 == 0.0_wp .or. d1 == 0.0_wp) return
       if ((d0 > 0.0_wp) .neqv. (d1 > 0.0_wp)) return
       w0 = 2*h1+h0; w1 = h1+2*h0
       m = (w0+w1)/(w0/d0+w1/d1)
    end if
  end function cooling_log_inverse_temperature_slope

  real(kind=wp) function cooling_log_inverse_temperature_end_slope(h0,h1,d0,d1) result(m)
    real(kind=wp), intent(in) :: h0,h1,d0,d1
    m = ((2*h0+h1)*d0-h0*d1)/(h0+h1)
    if ((m > 0.0_wp) .neqv. (d0 > 0.0_wp)) then
       m = 0.0_wp
    else if ((d0 > 0.0_wp) .neqv. (d1 > 0.0_wp)) then
       if (abs(m) > abs(3*d0)) m = 3*d0
    end if
  end function cooling_log_inverse_temperature_end_slope

  !=========================================================================
  ! Low-density reference from the ne-floor column, using the same T rule.
  !=========================================================================
  real(kind=wp) function low_density_cooling(idx, T) result(val)
    implicit none
    integer,       intent(in) :: idx
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: w1, wt
    integer :: it
    call cooling_temperature_cell(idx, T, it, w1, wt)
    val = col_cooling(idx, it, 1, w1, wt)
  end function low_density_cooling

  !=========================================================================
  ! Locate a temperature cell on this ion's source-limited log-T grid.
  !=========================================================================
  subroutine cooling_temperature_cell(idx, T, it, w1, wt)
    implicit none
    integer, intent(in) :: idx
    real(kind=wp), intent(in) :: T
    integer, intent(out) :: it
    real(kind=wp), intent(out) :: w1, wt
    real(kind=wp) :: lt
    integer :: ntemp
    if (.not. ieee_is_finite(T)) &
       error stop 'NLCOOL: temperature must be positive and finite'
    if (T <= 0.0_wp) &
       error stop 'NLCOOL: temperature must be positive and finite'
    if (log(T) > reg(idx)%source_ceiling_ln_K) then
       write(*,'(a,a,a,i0,a,es16.8)') ' NLCOOL: ', trim(reg(idx)%name), &
            ' stage ', reg(idx)%stage, ' unsupported temperature (K): ', T
       error stop 'NLCOOL: temperature exceeds source collision table'
    end if
    ntemp = size(reg(idx)%logT)
    lt = min(max(log10(T), reg(idx)%logT(1)), reg(idx)%logT(ntemp))
    it = min(int((lt - T_LO)/DT_L) + 1, ntemp-1)
    w1 = (1.0_wp/10.0_wp**lt - 1.0_wp/10.0_wp**reg(idx)%logT(it)) &
       / (1.0_wp/10.0_wp**reg(idx)%logT(it+1) - &
          1.0_wp/10.0_wp**reg(idx)%logT(it))
    wt = min(max((lt - reg(idx)%logT(it))/ &
         (reg(idx)%logT(it+1) - reg(idx)%logT(it)), 0.0_wp), 1.0_wp)
  end subroutine cooling_temperature_cell

  !=========================================================================
  ! Cooling ratio for the setup report. A zero reference has no defined ratio.
  !=========================================================================
  real(kind=wp) function suppression_ratio(idx, T, ne) result(S)
    implicit none
    integer,       intent(in) :: idx
    real(kind=wp), intent(in) :: T, ne
    S = low_density_cooling(idx, T)
    if (S > 0.0_wp) then
       S = cooling_coeff(idx, T, ne)/S
    else
       S = 0.0_wp
    end if
  end function suppression_ratio

  !=========================================================================
  ! Rank-0 summary of the table build (ions tabulated + wall-clock).
  !=========================================================================
  subroutine nlevel_cooling_report()
    implicit none
    if (mpar%p_rank /= 0) return
    if (nion_reg == 0) then
       write(*,'(a)') ' NLCOOL: no n-level cooling tables built '// &
          '(all ions in the low-density limit).'
       return
    end if
    write(*,'(a,i0,a,f8.3,a)') ' NLCOOL: ', nion_reg, &
       ' ion cooling tables built in ', build_time_total, ' s.'
  end subroutine nlevel_cooling_report

end module nlevel_cooling_mod
