module dust_temp_mod
!---------------------------------------------------------------------------
! MoCHII: equilibrium dust temperature and IR emission consuming Heat_dust
! (PLAN section 7 item 1, first consumer).
!
! The transported band (par%eion_min may be lowered into the FUV, e.g.
! 6 eV — gas cross sections vanish below their thresholds automatically,
! so only dust absorbs there) yields the grain heating rate
! Heat_dust(leaf) [erg s^-1 cm^-3].  Equilibrium MIXTURE temperature:
!     Heat_dust = (rhokap/d_cm) * femit(T_d),
!     femit(T) = 4 pi Int s_abs(nu) B_nu(T) dnu   [erg s^-1 cm^-2]
! with s_abs(nu) = (1-albedo) C_ext(nu)/C_ext(lambda_ref) from the SAME
! kext table (par%ion_dust_kext), integrated over the full tabulated
! range and precomputed on a log-T grid at setup; T_d by interpolation.
! This is the single equilibrium mixture temperature — the SEDust
! stochastic/PAH treatment (size- and material-resolved) is the next
! stage; outputs are labeled accordingly.
!
! Outputs: T_dust(leaf) in the rates file (via gas_rates_write hook) and
! '<base>_dustir.txt' — the grid-integrated modified-blackbody IR
! spectrum L_nu = sum_leaf 4 pi n_d C_abs(nu) B_nu(T_d) V, whose
! frequency integral equals sum(Heat_dust V) by construction (energy
! closure check printed).
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: dust_temp_setup, dust_temp_compute, dust_ir_write, t_dust

  integer, parameter :: NT = 61, NLAM_IR = 240
  real(kind=wp) :: tgrid(NT), femit(NT)
  real(kind=wp), allocatable :: lam_ir(:), sabs_ir(:)
  real(kind=wp), allocatable :: t_dust(:)

contains

  !=========================================================================
  subroutine dust_temp_setup()
    use mpi
    implicit none
    real(kind=wp), allocatable :: lam(:), alb(:), cext(:)
    character(len=512) :: line
    real(kind=wp) :: v(4), cref, T, nu1, nu2, bnu, snu, w
    integer :: unit, ios, n, i, k, ierr

    !--- read the kext table (lambda, albedo, -, C_ext); sort ascending.
    open(newunit=unit, file=trim(par%ion_dust_kext), status='old', iostat=ios)
    if (ios /= 0) then
       if (mpar%p_rank == 0) write(*,'(a)') &
          'ERROR: dust_temp needs par%ion_dust_kext.'
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    n = 0
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       read(line,*,iostat=ios) v
       if (ios == 0 .and. v(1) > 0.0_wp) n = n + 1
    end do
    allocate(lam(n), alb(n), cext(n))
    rewind(unit)
    i = 0
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       read(line,*,iostat=ios) v
       if (ios == 0 .and. v(1) > 0.0_wp) then
          i = i + 1
          lam(i) = v(1);  alb(i) = v(2);  cext(i) = v(4)
       end if
    end do
    close(unit)
    call sort3l(lam, alb, cext, n)
    cref = loginterp(lam, cext, n, par%lambda_ref)

    !--- IR wavelength grid (1 - 1000 um) with s_abs from the table.
    allocate(lam_ir(NLAM_IR), sabs_ir(NLAM_IR))
    do k = 1, NLAM_IR
       lam_ir(k) = exp(log(1.0_wp) + (log(1000.0_wp) - log(1.0_wp)) &
                   *real(k-1,wp)/real(NLAM_IR-1,wp))
       sabs_ir(k) = (1.0_wp - lininterp(lam, alb, n, lam_ir(k))) &
                    * loginterp(lam, cext, n, lam_ir(k)) / cref
    end do

    !--- femit(T) = 4 pi Int s_abs B_nu dnu over the FULL table range
    !--- (integrated on the table grid itself, trapezoid in nu).
    do i = 1, NT
       tgrid(i) = 10.0_wp**(0.5_wp + 3.0_wp*real(i-1,wp)/real(NT-1,wp))  ! 3.2-3200 K
       T = tgrid(i)
       femit(i) = 0.0_wp
       do k = 1, n-1
          nu1 = clight_cgs/(lam(k+1)*1.0e-4_wp)     ! ascending lambda -> nu2>nu1? careful
          nu2 = clight_cgs/(lam(k)*1.0e-4_wp)
          snu = 0.5_wp*( (1.0_wp-alb(k))*cext(k) &
                        + (1.0_wp-alb(k+1))*cext(k+1) )/cref
          bnu = 0.5_wp*(planck_nu(nu1, T) + planck_nu(nu2, T))
          femit(i) = femit(i) + fourpi*snu*bnu*(nu2 - nu1)
       end do
    end do
    deallocate(lam, alb, cext)
    if (mpar%p_rank == 0) write(*,'(a,es11.4,a)') &
       ' DUST: femit(20 K) = ', femit_of(20.0_wp), &
       ' erg/s/cm^2 (equilibrium mixture emission function ready)'
  end subroutine dust_temp_setup

  !=========================================================================
  elemental real(kind=wp) function planck_nu(nu, T) result(b)
    real(kind=wp), intent(in) :: nu, T
    real(kind=wp) :: x
    x = h_planck_cgs*nu/(kboltz_cgs*T)
    if (x < 700.0_wp) then
       b = 2.0_wp*h_planck_cgs*nu**3/clight_cgs**2/(exp(x) - 1.0_wp)
    else
       b = 0.0_wp
    end if
  end function planck_nu

  real(kind=wp) function femit_of(T) result(f)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lt, w
    integer :: i
    lt = min(max(log10(T), log10(tgrid(1))), log10(tgrid(NT)))
    i  = min(int((lt - log10(tgrid(1)))/(3.0_wp/real(NT-1,wp))) + 1, NT-1)
    w  = (lt - log10(tgrid(i)))/(log10(tgrid(i+1)) - log10(tgrid(i)))
    f  = exp(log(femit(i))*(1.0_wp - w) + log(femit(i+1))*w)
  end function femit_of

  !=========================================================================
  ! Invert Heat_dust = (rhokap/d_cm) femit(T) per leaf (log bisection).
  !=========================================================================
  subroutine dust_temp_compute(heat_dust)
    use octree_mod, only : amr_grid
    implicit none
    real(kind=wp), intent(in) :: heat_dust(:)
    real(kind=wp) :: kd, target, tlo, thi, tm
    integer :: il, it

    if (.not. allocated(t_dust)) allocate(t_dust(amr_grid%nleaf))
    if (size(t_dust) /= amr_grid%nleaf) then
       deallocate(t_dust);  allocate(t_dust(amr_grid%nleaf))
    end if
    t_dust = 0.0_wp
    do il = 1, amr_grid%nleaf
       kd = amr_grid%rhokap(il)/par%distance2cm     ! n_d C_ext(ref) [cm^-1]
       if (kd <= 0.0_wp .or. heat_dust(il) <= 0.0_wp) cycle
       target = heat_dust(il)/kd
       tlo = tgrid(1);  thi = tgrid(NT)
       if (femit_of(thi) <= target) then
          t_dust(il) = thi
          cycle
       end if
       do it = 1, 60
          tm = sqrt(tlo*thi)
          if (femit_of(tm) < target) then
             tlo = tm
          else
             thi = tm
          end if
          if (thi/tlo - 1.0_wp < 1.0e-5_wp) exit
       end do
       t_dust(il) = sqrt(tlo*thi)
    end do
  end subroutine dust_temp_compute

  !=========================================================================
  subroutine dust_ir_write(heat_dust)
    use octree_mod, only : amr_grid
    use utility,    only : get_base_name
    implicit none
    real(kind=wp), intent(in) :: heat_dust(:)
    real(kind=wp) :: Lnu(NLAM_IR), vol, kd, nu, Lir, Labs, dnu
    character(len=192) :: outname
    integer :: il, k, unit, ic

    if (mpar%p_rank /= 0) return
    Lnu = 0.0_wp;  Labs = 0.0_wp
    do il = 1, amr_grid%nleaf
       if (t_dust(il) <= 0.0_wp) cycle
       ic  = amr_grid%icell_of_leaf(il)
       vol = (2.0_wp*amr_grid%ch(ic)*par%distance2cm)**3
       kd  = amr_grid%rhokap(il)/par%distance2cm
       Labs = Labs + heat_dust(il)*vol
       do k = 1, NLAM_IR
          nu = clight_cgs/(lam_ir(k)*1.0e-4_wp)
          Lnu(k) = Lnu(k) + fourpi*kd*sabs_ir(k)*planck_nu(nu, t_dust(il))*vol
       end do
    end do
    !--- energy closure: integrate L_nu over nu (trapezoid, descending nu)
    Lir = 0.0_wp
    do k = 1, NLAM_IR-1
       dnu = clight_cgs/(lam_ir(k)*1.0e-4_wp) - clight_cgs/(lam_ir(k+1)*1.0e-4_wp)
       Lir = Lir + 0.5_wp*(Lnu(k) + Lnu(k+1))*dnu
    end do

    outname = trim(get_base_name(par%out_file))//'_dustir.txt'
    open(newunit=unit, file=trim(outname), status='replace')
    write(unit,'(a)') '# MoCHII equilibrium-mixture dust IR spectrum '// &
       '(single T_d per leaf; SEDust stochastic/PAH bands are the next stage)'
    write(unit,'(a,es14.6,a,es14.6,a,f7.3)') '# L_abs(dust) = ', Labs, &
       ' erg/s;  integral L_nu dnu = ', Lir, ';  ratio = ', Lir/max(Labs,tinest)
    write(unit,'(a)') '# lambda[um]   L_nu[erg/s/Hz]'
    do k = 1, NLAM_IR
       write(unit,'(f12.4,es15.6)') lam_ir(k), Lnu(k)
    end do
    close(unit)
    write(*,'(2a)') ' DUST: IR spectrum written to: ', trim(outname)
    write(*,'(a,f7.3,a)') ' DUST: IR/absorbed energy closure = ', &
       Lir/max(Labs,tinest), ' (1-2% from the 1-1000 um grid cut expected)'
  end subroutine dust_ir_write

  !=========================================================================
  subroutine sort3l(x, y1, y2, n)
    implicit none
    integer,       intent(in)    :: n
    real(kind=wp), intent(inout) :: x(n), y1(n), y2(n)
    real(kind=wp) :: tx, t1, t2
    integer :: i, j
    do i = 2, n
       tx = x(i);  t1 = y1(i);  t2 = y2(i)
       j = i - 1
       do while (j >= 1)
          if (x(j) <= tx) exit
          x(j+1) = x(j);  y1(j+1) = y1(j);  y2(j+1) = y2(j)
          j = j - 1
       end do
       x(j+1) = tx;  y1(j+1) = t1;  y2(j+1) = t2
    end do
  end subroutine sort3l

  real(kind=wp) function loginterp(x, y, n, x0) result(y0)
    integer,       intent(in) :: n
    real(kind=wp), intent(in) :: x(n), y(n), x0
    integer :: i
    real(kind=wp) :: w
    if (x0 <= x(1)) then
       y0 = y(1);  return
    else if (x0 >= x(n)) then
       y0 = y(n);  return
    end if
    do i = 1, n-1
       if (x0 < x(i+1)) exit
    end do
    w = log(x0/x(i))/log(x(i+1)/x(i))
    y0 = exp(log(max(y(i),tinest))*(1.0_wp-w) + log(max(y(i+1),tinest))*w)
  end function loginterp

  real(kind=wp) function lininterp(x, y, n, x0) result(y0)
    integer,       intent(in) :: n
    real(kind=wp), intent(in) :: x(n), y(n), x0
    integer :: i
    real(kind=wp) :: w
    if (x0 <= x(1)) then
       y0 = y(1);  return
    else if (x0 >= x(n)) then
       y0 = y(n);  return
    end if
    do i = 1, n-1
       if (x0 < x(i+1)) exit
    end do
    w = (x0 - x(i))/(x(i+1) - x(i))
    y0 = y(i)*(1.0_wp-w) + y(i+1)*w
  end function lininterp

end module dust_temp_mod
