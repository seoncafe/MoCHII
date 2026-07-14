module gas_opacity_mod
!---------------------------------------------------------------------------
! MoCHII: ionizing-band absorption coefficients.
!
! kap_ion(inu, leaf) is the absorption coefficient per CODE LENGTH in
! ionizing bin inu:
!   kap = nH * [ x_HI sigma_HI(E) +
!                y_He (x_HeI sigma_HeI(E) + x_HeII sigma_HeII(E)) ]
!         * distance2cm   (+ grey dust rhokap when par%ion_add_dust)
! with y_He = par%He_abund = n_He/n_H.  The grey s_ext rescale of the dust
! SED band does NOT apply here: opacity varies with both leaf and bin, so
! the band carries the full (nnu, nleaf) array.
!
! With gas_niter = 0 the array is filled once from the fixed gas state; the
! iteration refills it from the updated ion fractions between transport
! passes.
!---------------------------------------------------------------------------
  use define
  use octree_mod,    only : amr_grid
  use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, gas_nleaf
  use photo_xsec,    only : sigma_HI, sigma_HeI, sigma_HeII
  use ion_band_mod,  only : ion_e
  use memory_mod,    only : create_shared_mem
  implicit none
  private

  public :: gas_opacity_setup, gas_opacity_fill, kap_ion
  public :: ion_dust_sabs, ion_dust_ssca, ion_dust_g

  real(kind=wp), pointer :: kap_ion(:,:) => null()   ! (nnu_ion, nleaf)

  !--- ionizing-band dust ABSORPTION relative to the reference extinction:
  !--- ion_dust_sabs(b) = (1 - albedo(E_b)) C_ext(E_b) / C_ext(lambda_ref),
  !--- read from par%ion_dust_kext (MoCafe kext table format: lambda[um],
  !--- albedo, <cos>, C_ext/H, ...; D03 and astrodust layouts share these
  !--- first four columns).  The leaf grey rhokap (reference extinction per
  !--- code length) times this ratio gives the band absorption.  EUV dust
  !--- SCATTERING (albedo ~ 0.2-0.45) is deferred: the analytic edge walk
  !--- has no interaction sampling yet, so the scattered fraction is
  !--- treated as unabsorbed (dropped from the dust budget, not
  !--- redistributed) UNLESS par%ion_dust_scatter, in which case
  !--- ion_dust_ssca(b) = albedo C_ext/C_ext(ref) joins the extinction
  !--- and interactions are sampled (HG asymmetry ion_dust_g(b) from the
  !--- table <cos> column).
  real(kind=wp), allocatable :: ion_dust_sabs(:)
  real(kind=wp), allocatable :: ion_dust_ssca(:)
  real(kind=wp), allocatable :: ion_dust_g(:)

contains

  !=========================================================================
  subroutine gas_opacity_setup()
    use mpi
    implicit none
    if (par%ion_add_dust) call ion_dust_read()
    call create_shared_mem(kap_ion, [par%nnu_ion, gas_nleaf])
    call gas_opacity_fill()
    if (mpar%p_rank == 0) write(*,'(a,i4,a,i12,a)') &
       ' ION: opacity array kap_ion(', par%nnu_ion, ' bins, ', gas_nleaf, ' leaves) filled'
  end subroutine gas_opacity_setup

  !=========================================================================
  subroutine ion_dust_read()
    use mpi
    implicit none
    real(kind=wp), allocatable :: lam(:), alb(:), cext(:), gcos(:)
    character(len=512) :: line
    real(kind=wp) :: v(4), lref, cref, aref, lb, w
    integer :: unit, ios, n, i, b, ierr

    if (allocated(ion_dust_sabs)) return
    allocate(ion_dust_sabs(par%nnu_ion), ion_dust_ssca(par%nnu_ion), &
             ion_dust_g(par%nnu_ion))
    open(newunit=unit, file=trim(par%ion_dust_kext), status='old', iostat=ios)
    if (ios /= 0) then
       if (mpar%p_rank == 0) write(*,'(2a)') &
          'ERROR: ion_add_dust needs par%ion_dust_kext; cannot open ', &
          trim(par%ion_dust_kext)
       call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
    end if
    !--- count numeric rows, then read (lambda, albedo, -, C_ext)
    n = 0
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       read(line,*,iostat=ios) v
       if (ios == 0 .and. v(1) > 0.0_wp) n = n + 1
    end do
    allocate(lam(n), alb(n), cext(n), gcos(n))
    rewind(unit)
    i = 0
    do
       read(unit,'(a)',iostat=ios) line
       if (ios /= 0) exit
       read(line,*,iostat=ios) v
       if (ios == 0 .and. v(1) > 0.0_wp) then
          i = i + 1
          lam(i) = v(1);  alb(i) = v(2);  gcos(i) = v(3);  cext(i) = v(4)
       end if
    end do
    close(unit)
    !--- sort ascending in lambda (the D03 table is NOT monotonic: it has
    !--- an out-of-order seam, so a simple reversal is not enough).
    call sort4(lam, alb, gcos, cext, n)
    cref = interp_log(lam, cext, n, par%lambda_ref)
    do b = 1, par%nnu_ion
       lb = 1.23984_wp/ion_e_of(b)          ! um
       ion_dust_sabs(b) = (1.0_wp - interp_lin(lam, alb, n, lb)) &
                          * interp_log(lam, cext, n, lb) / cref
       ion_dust_ssca(b) = interp_lin(lam, alb, n, lb) &
                          * interp_log(lam, cext, n, lb) / cref
       ion_dust_g(b)    = interp_lin(lam, gcos, n, lb)
    end do
    if (mpar%p_rank == 0) then
       write(*,'(2a)')       ' ION: dust kext table = ', trim(par%ion_dust_kext)
       write(*,'(a,es11.4)') ' ION: C_ext/H at lambda_ref     = ', cref
       write(*,'(a,2f8.3)')  ' ION: sabs at band edges        = ', &
          ion_dust_sabs(1), ion_dust_sabs(par%nnu_ion)
       if (par%ion_dust_scatter) write(*,'(a,2f8.3,a,2f7.3)') &
          ' ION: ssca at band edges        = ', &
          ion_dust_ssca(1), ion_dust_ssca(par%nnu_ion), &
          ',  g = ', ion_dust_g(1), ion_dust_g(par%nnu_ion)
    end if
    deallocate(lam, alb, cext, gcos)
  end subroutine ion_dust_read

  real(kind=wp) function ion_e_of(b) result(e)
    use ion_band_mod, only : ion_e
    integer, intent(in) :: b
    e = ion_e(b)
  end function ion_e_of

  subroutine sort4(x, y1, y2, y3, n)
    implicit none
    integer,       intent(in)    :: n
    real(kind=wp), intent(inout) :: x(n), y1(n), y2(n), y3(n)
    real(kind=wp) :: tx, t1, t2, t3
    integer :: i, j
    do i = 2, n                       ! insertion sort (n ~ 1e3, setup only)
       tx = x(i);  t1 = y1(i);  t2 = y2(i);  t3 = y3(i)
       j = i - 1
       do while (j >= 1)
          if (x(j) <= tx) exit
          x(j+1) = x(j);  y1(j+1) = y1(j);  y2(j+1) = y2(j);  y3(j+1) = y3(j)
          j = j - 1
       end do
       x(j+1) = tx;  y1(j+1) = t1;  y2(j+1) = t2;  y3(j+1) = t3
    end do
  end subroutine sort4

  real(kind=wp) function interp_log(x, y, n, x0) result(y0)
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
  end function interp_log

  real(kind=wp) function interp_lin(x, y, n, x0) result(y0)
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
  end function interp_lin

  !=========================================================================
  ! (Re)fill kap_ion from the current gas state.  Called at setup and after
  ! every equilibrium update (opacity feedback: the stellar tally must be
  ! recomputed each iteration).
  !=========================================================================
  subroutine gas_opacity_fill()
    use mpi
    use physics_amr_mod, only : laursen09_ndust
    use species_mod,     only : n_elements, species_opacity_add
    implicit none
    integer  :: il, inu, ierr
    real(kind=wp) :: sHI, sHeI, sHeII, kap

    if (mpar%h_rank == 0) then
       !--- dust density tied to the COMPUTED
       !--- ionization state, refreshed every iteration; the PAH share
       !--- f_pah carries its own ionized-gas survival f_ion_pah
       !--- (default 0 = PAHs destroyed in ionized gas).
       if (trim(par%dust_model) == 'laursen09_live') then
          block
            real(kind=wp) :: xh, fac
            do il = 1, gas_nleaf
               xh  = gas_xHI(il)
               fac = (1.0_wp - par%f_pah)*(xh + par%f_ion_dust*(1.0_wp - xh)) &
                     + par%f_pah*(xh + par%f_ion_pah*(1.0_wp - xh))
               amr_grid%rhokap(il) = (par%Z_global/max(par%Z_ref,1.0e-30_wp)) &
                  * gas_nH(il)*fac * par%cext_dust * par%DGR * par%distance2cm
            end do
          end block
       end if
       do inu = 1, par%nnu_ion
          sHI   = sigma_HI(ion_e(inu))
          sHeI  = sigma_HeI(ion_e(inu))
          sHeII = sigma_HeII(ion_e(inu))
          do il = 1, gas_nleaf
             kap = gas_nH(il) * ( gas_xHI(il)*sHI &
                   + par%He_abund*(gas_xHeI(il)*sHeI + gas_xHeII(il)*sHeII) ) &
                   * par%distance2cm
             if (par%ion_add_dust) then
                kap = kap + amr_grid%rhokap(il)*ion_dust_sabs(inu)
                !--- scattering joins the EXTINCTION when interactions
                !--- are sampled; otherwise scattered photons pass.
                if (par%ion_dust_scatter) &
                   kap = kap + amr_grid%rhokap(il)*ion_dust_ssca(inu)
             end if
             kap_ion(inu, il) = kap
          end do
       end do
       !--- metal photoionization absorption (the only GAS opacity in the
       !--- FUV bins; negligible next to H/He above 13.6 eV).
       if (par%use_metals .and. par%ion_metal_abs .and. n_elements > 0) &
          call species_opacity_add(kap_ion, par%nnu_ion, gas_nleaf)
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)
  end subroutine gas_opacity_fill

end module gas_opacity_mod
