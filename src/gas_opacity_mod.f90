module gas_opacity_mod
!---------------------------------------------------------------------------
! MoCHII: ionizing-band dust optics and the leaf state the packet opacity is
! built from.
!
! A packet's absorption coefficient per CODE LENGTH is evaluated at its own
! photon energy by ion_packet_opacity (ion_packet_mod):
!   kap = nH * [ x_HI sigma_HI(E) +
!                y_He (x_HeI sigma_HeI(E) + x_HeII sigma_HeII(E)) + metals ]
!         * distance2cm   (+ grey dust rhokap when par%ion_add_dust)
! with y_He = par%He_abund = n_He/n_H.  This module owns the two leaf-resident
! ingredients of that expression: the dust extinction table (read once, then
! evaluated at the exact energy through ion_dust_*_at), and the grey grain
! extinction density rhokap, which follows the live ionization state.
! gas_opacity_fill refreshes them between transport passes; the grey s_ext
! rescale of the dust SED band does NOT apply here.
!---------------------------------------------------------------------------
  use define
  use octree_mod,    only : amr_grid
  use gas_state_mod, only : gas_nH, gas_xHI, gas_nleaf
  use ion_band_mod,  only : nnu_band
  implicit none
  private

  public :: gas_opacity_setup, gas_opacity_fill
  public :: ion_dust_sabs, ion_dust_ssca, ion_dust_g
  public :: ion_dust_sabs_at, ion_dust_ssca_at, ion_dust_g_at

  !--- ionizing-band dust ABSORPTION relative to the reference extinction:
  !--- ion_dust_sabs(b) = (1 - albedo(E_b)) C_ext(E_b) / C_ext(lambda_ref).
  !--- The (lambda[um], albedo, <cos>, C_ext/H) optics come from one of two
  !--- sources: par%ion_dust_kext (a kext file: MoCafe table format, D03 and
  !--- astrodust layouts share these first four columns) when it is set, which
  !--- OVERRIDES the model; otherwise the SEDust grain model itself
  !--- (grain_model_mod), so extinction/absorption and the reemitted IR refer
  !--- to the same grains.  The leaf grey rhokap (reference extinction per
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

  !--- The sorted extinction table is retained (read once) so the dust cross
  !--- sections can be evaluated at a photon's exact energy, not only at the bin
  !--- centers (ion_dust_*_at below).  The bin arrays above are unchanged: they
  !--- are still filled by the inline loop in ion_dust_read.
  real(kind=wp), allocatable :: dust_tab_lam(:), dust_tab_alb(:)
  real(kind=wp), allocatable :: dust_tab_gcos(:), dust_tab_cext(:)
  integer       :: dust_tab_n    = 0
  real(kind=wp) :: dust_tab_cref = 0.0_wp

contains

  !=========================================================================
  subroutine gas_opacity_setup()
    use mpi
    implicit none
    if (par%ion_add_dust) call ion_dust_read()
    call gas_opacity_fill()
  end subroutine gas_opacity_setup

  !=========================================================================
  subroutine ion_dust_read()
    use mpi
    use grain_model_mod, only : grain_extinction_table
    implicit none
    real(kind=wp), allocatable :: lam(:), alb(:), cext(:), gcos(:)
    character(len=512) :: line
    real(kind=wp) :: v(4), lref, cref, aref, lb, w, lam_short, lam_long
    integer :: unit, ios, n, i, b, ierr
    logical :: from_model

    if (allocated(ion_dust_sabs)) return
    allocate(ion_dust_sabs(nnu_band), ion_dust_ssca(nnu_band), &
             ion_dust_g(nnu_band))

    !--- optics source: a kext file OVERRIDES the model (backward compatible);
    !--- an empty par%ion_dust_kext takes the (lambda, albedo, <cos>, C_ext/H)
    !--- from the SEDust grain model so extinction and emission share grains.
    from_model = (len_trim(par%ion_dust_kext) == 0)
    if (.not. from_model) then
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
    else
       call grain_extinction_table(lam, alb, gcos, cext, n)
    end if
    !--- sort ascending in lambda (the D03 file is NOT monotonic: it has an
    !--- out-of-order seam, so a simple reversal is not enough; the model grid
    !--- is already ascending, so this is a no-op there).
    call sort4(lam, alb, gcos, cext, n)

    !--- the model grid must span lambda_ref AND the whole transported band,
    !--- because the interpolators clamp at the grid edges: a grid that stops
    !--- at 0.0912 um (e.g. Zubko) would silently freeze the EUV optics.
    if (from_model) then
       lam_short = 1.23984_wp/par%eion_max
       if (par%add_fuv) then
          lam_long = 1.23984_wp/par%efuv_min
       else
          lam_long = 1.23984_wp/par%eion_min
       end if
       if (lam(1) > lam_short*1.0001_wp .or. lam(n) < lam_long*0.9999_wp .or. &
           par%lambda_ref < lam(1) .or. par%lambda_ref > lam(n)) then
          if (mpar%p_rank == 0) then
             write(*,'(3a)') 'ERROR: par%dust_model_sed = ''', &
                trim(par%dust_model_sed), ''' does not cover the ionizing band;'
             write(*,'(a,es10.3,a,es10.3,a)') '       model grid [um] = [', &
                lam(1), ', ', lam(n), ']'
             write(*,'(a,es10.3,a,es10.3,a,es10.3)') &
                '       band + ref need [um] = [', lam_short, ', ', lam_long, &
                '] and ', par%lambda_ref
             write(*,'(a)') '       set par%ion_dust_kext to a table that covers it.'
          end if
          call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
       end if
    end if

    cref = interp_log(lam, cext, n, par%lambda_ref)
    !--- retain the sorted table for exact-energy evaluation (ion_dust_*_at);
    !--- the bin loop below is left untouched so ion_dust_*(b) stay identical.
    dust_tab_n = n;  dust_tab_cref = cref
    if (allocated(dust_tab_lam)) deallocate(dust_tab_lam, dust_tab_alb, &
       dust_tab_gcos, dust_tab_cext)
    allocate(dust_tab_lam(n), dust_tab_alb(n), dust_tab_gcos(n), dust_tab_cext(n))
    dust_tab_lam  = lam(1:n);   dust_tab_alb  = alb(1:n)
    dust_tab_gcos = gcos(1:n);  dust_tab_cext = cext(1:n)
    do b = 1, nnu_band
       lb = 1.23984_wp/ion_e_of(b)          ! um
       ion_dust_sabs(b) = (1.0_wp - interp_lin(lam, alb, n, lb)) &
                          * interp_log(lam, cext, n, lb) / cref
       ion_dust_ssca(b) = interp_lin(lam, alb, n, lb) &
                          * interp_log(lam, cext, n, lb) / cref
       ion_dust_g(b)    = interp_lin(lam, gcos, n, lb)
    end do
    !--- S0 invariant: the bin arrays are exactly the *_at evaluators at the bin
    !--- centers (same retained table, same expression).  Silent on success, so
    !--- the run output is unchanged; guards against future drift of either side.
    do b = 1, nnu_band
       if (abs(ion_dust_sabs(b) - ion_dust_sabs_at(ion_e_of(b))) > &
              1.0e-12_wp*max(abs(ion_dust_sabs(b)), tiny(1.0_wp)) .or. &
           abs(ion_dust_ssca(b) - ion_dust_ssca_at(ion_e_of(b))) > &
              1.0e-12_wp*max(abs(ion_dust_ssca(b)), tiny(1.0_wp)) .or. &
           abs(ion_dust_g(b)    - ion_dust_g_at(ion_e_of(b)))    > &
              1.0e-12_wp*max(abs(ion_dust_g(b)),    tiny(1.0_wp))) &
          error stop 'ion_dust_*_at disagrees with the bin array'
    end do
    if (mpar%p_rank == 0) then
       if (from_model) then
          write(*,'(3a,i0,a)') ' ION: dust optics = ', trim(par%dust_model_sed), &
             ' (SEDust grain model, ', n, ' wavelengths)'
          write(*,'(a,es10.3,a,es10.3,a,f7.1,a)') &
             ' ION: model grid [um]           = [', lam(1), ', ', lam(n), &
             '] (to ', 1.23984_wp/lam(1), ' eV)'
       else
          write(*,'(2a)')       ' ION: dust kext table = ', trim(par%ion_dust_kext)
       end if
       write(*,'(a,es11.4)') ' ION: C_ext/H at lambda_ref     = ', cref
       write(*,'(a,2f8.3)')  ' ION: sabs at band edges        = ', &
          ion_dust_sabs(1), ion_dust_sabs(nnu_band)
       if (par%ion_dust_scatter) write(*,'(a,2f8.3,a,2f7.3)') &
          ' ION: ssca at band edges        = ', &
          ion_dust_ssca(1), ion_dust_ssca(nnu_band), &
          ',  g = ', ion_dust_g(1), ion_dust_g(nnu_band)
    end if
    deallocate(lam, alb, cext, gcos)
  end subroutine ion_dust_read

  real(kind=wp) function ion_e_of(b) result(e)
    use ion_band_mod, only : ion_e
    integer, intent(in) :: b
    e = ion_e(b)
  end function ion_e_of

  !--- Dust cross sections at a photon's exact energy E [eV] from the retained
  !--- extinction table, mirroring the inline bin-loop expressions in
  !--- ion_dust_read: ion_dust_*(b) equals ion_dust_*_at(ion_e(b)) by
  !--- construction.  Not yet wired into transport (isolated S0 addition).
  real(kind=wp) function ion_dust_sabs_at(E) result(s)
    real(kind=wp), intent(in) :: E
    real(kind=wp) :: lb
    lb = 1.23984_wp/E
    s = (1.0_wp - interp_lin(dust_tab_lam, dust_tab_alb, dust_tab_n, lb)) &
        * interp_log(dust_tab_lam, dust_tab_cext, dust_tab_n, lb) / dust_tab_cref
  end function ion_dust_sabs_at

  real(kind=wp) function ion_dust_ssca_at(E) result(s)
    real(kind=wp), intent(in) :: E
    real(kind=wp) :: lb
    lb = 1.23984_wp/E
    s = interp_lin(dust_tab_lam, dust_tab_alb, dust_tab_n, lb) &
        * interp_log(dust_tab_lam, dust_tab_cext, dust_tab_n, lb) / dust_tab_cref
  end function ion_dust_ssca_at

  real(kind=wp) function ion_dust_g_at(E) result(g)
    real(kind=wp), intent(in) :: E
    real(kind=wp) :: lb
    lb = 1.23984_wp/E
    g = interp_lin(dust_tab_lam, dust_tab_gcos, dust_tab_n, lb)
  end function ion_dust_g_at

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
  ! (Re)build the leaf-resident ingredients of the packet opacity from the
  ! current gas state: the metal stage fractions and, for a live dust model,
  ! the grey grain extinction density.  Called at setup and after every
  ! equilibrium update (opacity feedback: the stellar tally must be recomputed
  ! each iteration).
  !=========================================================================
  subroutine gas_opacity_fill()
    use mpi
    use physics_amr_mod, only : laursen09_ndust, dust_opac_norm
    use species_mod,     only : n_elements, species_stage_cache_refresh
    implicit none
    integer  :: il, ierr

    if (par%use_metals .and. n_elements > 0) call species_stage_cache_refresh()
    if (mpar%h_rank == 0) then
       !--- dust density tied to the COMPUTED
       !--- ionization state, refreshed every iteration; the PAH share
       !--- f_pah carries its own ionized-gas survival f_ion_pah
       !--- (default 0 = PAHs destroyed in ionized gas).  dust_opac_norm is
       !--- the scale solved once from par%taumax/par%tauhomo at grid setup
       !--- (1 without a target); reapplying it here is what makes the
       !--- target the optical depth of the INITIAL state rather than a
       !--- value discarded at the first refill.  It is never re-solved, so
       !--- the realized tau then follows the grain survival factor freely.
       if (trim(par%dust_model) == 'laursen09_live') then
          block
            real(kind=wp) :: xh, fac
            do il = 1, gas_nleaf
               xh  = gas_xHI(il)
               fac = (1.0_wp - par%f_pah)*(xh + par%f_ion_dust*(1.0_wp - xh)) &
                     + par%f_pah*(xh + par%f_ion_pah*(1.0_wp - xh))
               amr_grid%rhokap(il) = (par%Z_global/max(par%Z_ref,1.0e-30_wp)) &
                  * gas_nH(il)*fac * par%cext_dust * par%DGR * par%distance2cm &
                  * dust_opac_norm
            end do
          end block
       end if
    end if
    call MPI_BARRIER(mpar%hostcomm, ierr)
  end subroutine gas_opacity_fill

end module gas_opacity_mod
