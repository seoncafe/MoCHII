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
! is a function of (T, n_e) alone and is tabulated exactly at setup from the
! n-level atom files (data/atomic/nlevel_<el>_<stage>.txt).  Lambda_SE
! decreases with n_e (the populations saturate) -- the density suppression.
! Stored for each ion:
!   Ssup(T, ne)   = Lambda_SE(T, ne) / Lambda_SE(T, ne_floor)  (2-D, in (0,1])
!   lam_se_ref(T) = Lambda_SE(T, ne_floor)                     (1-D, low-ne
!                                                               limit)
! with ne_floor = 10^-6 cm^-3 far below every critical density, so the
! reference column is the low-density limit.
!
! Runtime cooling coefficient, combined with the existing Tier-1 fit lam_fit
! passed in as an argument (so this module never uses cooling_mod -> no
! module cycle):
!   Lambda(T, ne) = lam_fit(T) * [ f_cov(T) Ssup(T, ne) + (1 - f_cov(T)) ]
!   f_cov(T)      = min(1, lam_se_ref(T) / lam_fit(T))
! f_cov is the fraction of the low-density cooling the (possibly truncated)
! n-level model accounts for; the remaining 1 - f_cov is high-excitation
! resonance-line cooling whose critical densities are far above nebular
! values, left in its low-density limit (unsuppressed).  As n_e -> 0,
! Ssup -> 1 and Lambda -> lam_fit exactly, recovering the current default.
! Ions with no n-level model get idx = 0 -> lam_fit returned unchanged.
!
! The table is built here from the SAME solver (nlevel_mod) that produces the
! emergent line diagnostics, so the thermal balance and the line list are one
! physics.  Interpolation: bilinear in (log T, log n_e) for Ssup, linear in
! log T for lam_se_ref, with edge clamping (log n_e below the floor -> 1,
! above the ceiling -> the last column, T clamped to the grid).
!---------------------------------------------------------------------------
  use define
  use nlevel_mod, only : nlevel_atom_type, nlevel_load, MAXLEV, &
                         nlevel_cooling_prepT, nlevel_cooling_solve
  implicit none
  private

  public :: nlevel_cooling_add, nlevel_cooling_apply
  public :: nlevel_cooling_lam_se, nlevel_cooling_report

  !--- (log T, log n_e) grid: log10 T = 2..6 step 0.01 (401 points); log10 n_e
  !--- = -6..10 step 0.1 (161 points).  ne_floor = 10^-6 is the reference
  !--- (low-density) column.  This is finer in T than the initial D4 spec
  !--- (0.05/0.2 dex): the gate G-a criterion (worst 0.5% vs the direct solve)
  !--- is set by the low-T/high-n_e corner where a high-Tx line and a saturated
  !--- low-n_crit line become comparable -- a sum of two exponentials in T that
  !--- one interpolation cell cannot resolve at 0.05 dex.  The build is a few
  !--- seconds for the full registry (gate G-h) and the runtime lookup cost is
  !--- independent of the grid size, so the T axis is refined to pass G-a with
  !--- margin (plan section C1: "refine the T axis if it fails").
  integer, parameter :: NT = 401, NNE = 161, MAXION = 40
  real(kind=wp), parameter :: T_LO  =  2.0_wp, DT_L  = 0.01_wp
  real(kind=wp), parameter :: NE_LO = -6.0_wp, DNE_L = 0.1_wp

  !--- Ssup/lamref are allocatable so a low_density run (no tables built)
  !--- carries none of the ~0.5 MB/ion table storage: the static reg(MAXION)
  !--- array then costs only its scalar members until nlevel_cooling_add
  !--- allocates the arrays of an ion actually tabulated.
  type cool_table_type
     character(len=16) :: name = ''
     integer :: stage = 0
     integer :: nlev  = 0
     real(kind=wp), allocatable :: Ssup(:,:)
     real(kind=wp), allocatable :: lamref(:)
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
  ! (par%atomic_dir, '_full' suffix for Fe under par%fe_levels_full, exactly
  ! as lines_mod), builds Ssup + lam_se_ref, and returns the new index.
  ! Returns 0 (no table; caller keeps lam_fit) when the file does not load.
  !=========================================================================
  integer function nlevel_cooling_add(name, stage) result(idx)
    implicit none
    character(len=*), intent(in) :: name
    integer,          intent(in) :: stage
    type(nlevel_atom_type) :: atom
    character(len=256) :: fname
    character(len=8)   :: nlsuf
    logical :: ok
    real(kind=wp) :: Rcoll(MAXLEV,MAXLEV), Rrad(MAXLEV,MAXLEV)
    real(kind=wp) :: lam(NT, NNE), Tk, nek
    integer(kind=i8b) :: c0, c1, crate
    integer :: i, j, ie

    idx = 0
    call nlevel_cooling_init()

    !--- filename (par%fe_levels_full appends '_full' for Fe; lines_mod:296).
    nlsuf = ''
    if (par%fe_levels_full .and. trim(name) == 'fe') nlsuf = '_full'
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

    call system_clock(c0, crate)
    do i = 1, NT
       Tk = 10.0_wp**gT(i)
       call nlevel_cooling_prepT(atom, Tk, Rcoll, Rrad)
       do j = 1, NNE
          nek = 10.0_wp**gne(j)
          lam(i,j) = nlevel_cooling_solve(atom, Rcoll, Rrad, nek)
       end do
    end do
    call system_clock(c1)
    build_time_total = build_time_total + real(c1-c0, wp)/real(crate, wp)

    nion_reg = nion_reg + 1
    ie = nion_reg
    reg(ie)%name  = name
    reg(ie)%stage = stage
    reg(ie)%nlev  = atom%nlev
    allocate(reg(ie)%Ssup(NT, NNE), reg(ie)%lamref(NT))
    reg(ie)%lamref(:) = lam(:,1)
    do j = 1, NNE
       reg(ie)%Ssup(:,j) = lam(:,j)/max(lam(:,1), tiny(1.0_wp))
    end do
    idx = ie

    if (mpar%p_rank == 0) write(*,'(a,a,a,i0,a,i0,a,f8.4)') &
       ' NLCOOL: ', trim(name), ' stage ', stage, &
       ' tabulated (nlev = ', atom%nlev, &
       '), Ssup(1e4 K, 1e4 cm^-3) = ', suppression_ratio(ie, 1.0e4_wp, 1.0e4_wp)
  end function nlevel_cooling_add

  !=========================================================================
  ! Runtime cooling coefficient (the D2 combiner).  idx <= 0 -> lam_fit.
  !=========================================================================
  real(kind=wp) function nlevel_cooling_apply(idx, lam_fit, T, ne) result(lam)
    implicit none
    integer,       intent(in) :: idx
    real(kind=wp), intent(in) :: lam_fit, T, ne
    real(kind=wp) :: fcov, S, lamref
    if (idx <= 0) then
       lam = lam_fit
       return
    end if
    lamref = low_density_cooling(idx, T)
    S      = cooling_coeff(idx, T, ne)/max(lamref, tiny(1.0_wp))
    fcov   = min(1.0_wp, lamref/max(lam_fit, tiny(1.0_wp)))
    lam    = lam_fit*(fcov*S + (1.0_wp - fcov))
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
  ! -> the last column; T clamped to the grid.
  !=========================================================================
  real(kind=wp) function cooling_coeff(idx, T, ne) result(lam)
    implicit none
    integer,       intent(in) :: idx
    real(kind=wp), intent(in) :: T, ne
    real(kind=wp) :: lt, lne, w1, wt, wn, L0, L1
    integer :: it, je
    lt = min(max(log10(T), gT(1)), gT(NT))
    it = min(int((lt - gT(1))/DT_L) + 1, NT-1)
    w1 = (1.0_wp/10.0_wp**lt       - 1.0_wp/10.0_wp**gT(it)) &
       / (1.0_wp/10.0_wp**gT(it+1) - 1.0_wp/10.0_wp**gT(it))
    wt = min(max((lt - gT(it))/(gT(it+1) - gT(it)), 0.0_wp), 1.0_wp)
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
  ! Cooling coefficient at (T, ne-column j): ln(Lambda) linear in 1/T between
  ! the two log-T rows (weight w1); a zero endpoint (never in practice: the
  ! lowest-Tx line does not underflow down to 100 K) falls back to linear in
  ! log T (weight wt).  Lambda_grid(i,j) = lamref(i) Ssup(i,j).
  !=========================================================================
  real(kind=wp) function col_cooling(idx, it, j, w1, wt) result(val)
    implicit none
    integer,       intent(in) :: idx, it, j
    real(kind=wp), intent(in) :: w1, wt
    real(kind=wp) :: Li, Lip
    Li  = reg(idx)%lamref(it)  *reg(idx)%Ssup(it,j)
    Lip = reg(idx)%lamref(it+1)*reg(idx)%Ssup(it+1,j)
    if (Li > 0.0_wp .and. Lip > 0.0_wp) then
       val = exp(log(Li) + w1*(log(Lip) - log(Li)))
    else
       val = Li*(1.0_wp - wt) + Lip*wt
    end if
  end function col_cooling

  !=========================================================================
  ! lam_se_ref(T): the low-density cooling limit (the ne-floor column), used
  ! for f_cov and the report; same ln-linear-in-1/T interpolation in T.
  !=========================================================================
  real(kind=wp) function low_density_cooling(idx, T) result(val)
    implicit none
    integer,       intent(in) :: idx
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lt, w1, wt
    integer :: it
    lt = min(max(log10(T), gT(1)), gT(NT))
    it = min(int((lt - gT(1))/DT_L) + 1, NT-1)
    w1 = (1.0_wp/10.0_wp**lt       - 1.0_wp/10.0_wp**gT(it)) &
       / (1.0_wp/10.0_wp**gT(it+1) - 1.0_wp/10.0_wp**gT(it))
    wt = min(max((lt - gT(it))/(gT(it+1) - gT(it)), 0.0_wp), 1.0_wp)
    val = col_cooling(idx, it, 1, w1, wt)
  end function low_density_cooling

  !=========================================================================
  ! Ssup(T, ne) = Lambda_SE(T, ne) / lam_se_ref(T), for the setup report.
  !=========================================================================
  real(kind=wp) function suppression_ratio(idx, T, ne) result(S)
    implicit none
    integer,       intent(in) :: idx
    real(kind=wp), intent(in) :: T, ne
    S = cooling_coeff(idx, T, ne)/max(low_density_cooling(idx, T), tiny(1.0_wp))
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
