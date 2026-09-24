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
! with y_He = par%abund_He = n_He/n_H.  This module owns the two leaf-resident
! ingredients of that expression: the dust extinction table (read once, then
! evaluated at the exact energy through ion_dust_*_at), and the grey grain
! extinction density rhokap, which follows the live ionization state.
! gas_opacity_fill refreshes them between transport passes; the grey s_ext
! rescale of the dust SED band does NOT apply here.
!---------------------------------------------------------------------------
  use define
  use octree_mod,    only : amr_grid
  use gas_state_mod, only : gas_nH, gas_xHI, gas_nleaf
  implicit none
  private

  public :: gas_opacity_setup, gas_opacity_fill, dust_reference_extinction
  public :: ion_dust_sabs_at, ion_dust_ssca_at, ion_dust_g_at

  !--- Ionizing-band dust optics relative to the reference extinction, at a
  !--- packet's exact energy E:
  !---   ion_dust_sabs_at(E) = (1 - albedo(E)) C_ext(E) / C_ext(lambda_ref),
  !---   ion_dust_ssca_at(E) = albedo(E) C_ext(E) / C_ext(lambda_ref),
  !---   ion_dust_g_at(E)    = <cos>(E), the Henyey-Greenstein asymmetry.
  !--- The (lambda[um], albedo, <cos>, C_ext/H) optics come from one of two
  !--- sources: par%ion_dust_kext (a kext file: MoCafe table format, D03 and
  !--- astrodust layouts share these first four columns) when it is set, which
  !--- OVERRIDES the model; otherwise the SEDust grain model itself
  !--- (grain_model_mod), so extinction/absorption and the reemitted IR refer
  !--- to the same grains.  The leaf grey rhokap (reference extinction per
  !--- code length) times these ratios gives the band opacity.  With
  !--- par%ion_dust_scatter the scattering joins the extinction and
  !--- interactions are sampled; without it the scattered fraction passes
  !--- (absorption-only treatment: dropped from the dust budget, not
  !--- redistributed).
  !--- The sorted extinction table is retained (read once) and evaluated at
  !--- each photon's energy.
  real(kind=wp), allocatable :: dust_tab_lam(:), dust_tab_alb(:)
  real(kind=wp), allocatable :: dust_tab_gcos(:), dust_tab_cext(:)
  integer       :: dust_tab_n    = 0
  real(kind=wp) :: dust_tab_cref = 0.0_wp

  !--- how far par%cext_dust, when a run states one outright, may sit from the
  !--- table's own C_ext/H at par%lambda_ref before the run says so.  What this
  !--- has to catch is a value left over from OTHER grains, which is O(1) --
  !--- astrodust to MRN is 1.68x -- so 1% is loose by two orders and still
  !--- unmissable.  Tighter would cry wolf on a deliberately rounded
  !--- normalization: the 19 inputs on the WD01/D03 kext file pin 4.868e-22,
  !--- the quoted Milky Way C_ext(V) per H, where that file interpolated at
  !--- 0.55 um gives 4.8606e-22 -- 0.15% apart, and the rounded number is the
  !--- one they mean.
  real(kind=wp), parameter :: CEXT_DUST_TOL = 1.0e-2_wp

contains

  !=========================================================================
  subroutine gas_opacity_setup()
    use mpi
    implicit none
    if (par%ion_add_dust) call ion_dust_read()
    call gas_opacity_fill()
  end subroutine gas_opacity_setup

  !=========================================================================
  !=========================================================================
  !--- The extinction of the grains this run carries, at par%lambda_ref, and
  !--- the (lambda, albedo, <cos>, C_ext/H) table it is interpolated from.
  !---
  !--- This is what par%cext_dust IS.  Every leaf's grey opacity is built as
  !--- n_dust C_ext(lambda_ref) per code length (grid_mod_amr), so that number
  !--- is a PROPERTY OF THE NAMED GRAINS, not an independent choice: it is the
  !--- same C_ext/H the band ratios below are divided by.  It used to be typed
  !--- into the namelist beside par%dust_model, which gave one quantity two
  !--- sources -- and the wrong one is silent, because a stale value does not
  !--- fail, it just scales every leaf's dust optical depth.  So par%cext_dust
  !--- < 0, the default, means "the grains this run carries", and a value given
  !--- outright is still taken as given, with both numbers printed so that a
  !--- deliberate renormalization stays visible and an accidental one does not.
  !---
  !--- Called BEFORE the grid is built (main.f90), because grid_create_amr
  !--- multiplies par%cext_dust into every leaf's rhokap.  Loading the table
  !--- depends on nothing but par% and the files it names, so that order is
  !--- free; the band tabulation ion_dust_read does with the same table needs
  !--- the band, which needs the grid, and stays where it is.
  !---
  !--- A run with no dust in the band loads no optics and has no reference
  !--- extinction to speak of: par%cext_dust < 0 becomes 0 there, and the grey
  !--- rhokap it would otherwise scale is read by nothing (the dust heating,
  !--- the SED and the packet opacity are all under par%ion_add_dust).
  subroutine dust_reference_extinction()
    use mpi
    use grain_model_mod, only : grain_extinction_table, dmodel
    implicit none
    real(kind=wp), allocatable :: lam(:), alb(:), cext(:), gcos(:)
    character(len=512) :: line
    real(kind=wp) :: v(4), cref, lam_short, lam_long, given
    integer :: unit, ios, n, i, ierr
    logical :: from_model

    if (dust_tab_n > 0) return

    if (.not. par%ion_add_dust) then
       if (par%cext_dust < 0.0_wp) par%cext_dust = 0.0_wp
       return
    end if

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
    !--- because the interpolators clamp at the grid edges: an unextended
    !--- astrodust or DL07 grid stops at the T-matrix Q table's 0.0912 um end
    !--- and would silently freeze the EUV optics at that value.  (Zubko never
    !--- does: its own DustEM optics already reach 0.001 um.)
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
             write(*,'(3a)') 'ERROR: par%dust_model = ''', &
                trim(par%dust_model), ''' does not cover the ionizing band;'
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
    !--- retain the sorted table: the bin arrays and the exact-energy
    !--- evaluators (ion_dust_*_at) are both built on it, so they cannot come
    !--- from two different reads.
    dust_tab_n = n;  dust_tab_cref = cref
    if (allocated(dust_tab_lam)) deallocate(dust_tab_lam, dust_tab_alb, &
       dust_tab_gcos, dust_tab_cext)
    allocate(dust_tab_lam(n), dust_tab_alb(n), dust_tab_gcos(n), dust_tab_cext(n))
    dust_tab_lam  = lam(1:n);   dust_tab_alb  = alb(1:n)
    dust_tab_gcos = gcos(1:n);  dust_tab_cext = cext(1:n)
    deallocate(lam, alb, cext, gcos)

    !--- resolve par%cext_dust against the table just loaded.
    given = par%cext_dust
    if (given < 0.0_wp) par%cext_dust = cref

    if (mpar%p_rank == 0) then
       if (from_model) then
          write(*,'(3a,i0,a)') ' ION: dust optics = ', trim(par%dust_model), &
             ' (SEDust grain model, ', n, ' wavelengths)'
          write(*,'(a,es10.3,a,es10.3,a,f7.1,a)') &
             ' ION: model grid [um]           = [', dust_tab_lam(1), ', ', &
             dust_tab_lam(n), '] (to ', 1.23984_wp/dust_tab_lam(1), ' eV)'
          !--- the model name does not identify the optics: the curve above is
          !--- the size-integrated table the builder loaded, and which table it
          !--- was decides the extinction this run transports.
          write(*,'(2a)') ' ION: extinction table          = ', &
             trim(dmodel%kext_path)
       else
          write(*,'(2a)')       ' ION: dust kext table = ', trim(par%ion_dust_kext)
       end if
       write(*,'(a,es11.4)') ' ION: C_ext/H at lambda_ref     = ', cref
       if (given < 0.0_wp) then
          write(*,'(a)') ' ION: cext_dust                 = that value (par%cext_dust < 0)'
       else
          write(*,'(a,es11.4,a,f9.5,a)') &
             ' ION: cext_dust                 = ', given, ' as given (x', &
             given/max(cref, tiny(1.0_wp)), ' of the table)'
          if (abs(given/max(cref, tiny(1.0_wp)) - 1.0_wp) > CEXT_DUST_TOL) then
             write(*,'(a)') ' ION: WARNING: par%cext_dust is not this table''s '// &
                'C_ext/H at lambda_ref, so the'
             write(*,'(a)') ' ION:          grid dust optical depth is scaled by '// &
                'the ratio above.  Leave'
             write(*,'(a)') ' ION:          par%cext_dust negative to take the '// &
                'grains'' own value.'
          end if
       end if
    end if
  end subroutine dust_reference_extinction

  !=========================================================================
  !--- Load the extinction table the band optics are evaluated from, and
  !--- report them at the two edges of the transported band.
  subroutine ion_dust_read()
    use ion_band_mod, only : ion_eedge
    implicit none
    logical, save :: done = .false.
    real(kind=wp) :: elo, ehi

    if (done) return
    done = .true.
    call dust_reference_extinction()
    elo = ion_eedge(1);  ehi = ion_eedge(size(ion_eedge))
    if (mpar%p_rank == 0) then
       write(*,'(a,2f8.3)')  ' ION: sabs at band edges        = ', &
          ion_dust_sabs_at(elo), ion_dust_sabs_at(ehi)
       if (par%ion_dust_scatter) write(*,'(a,2f8.3,a,2f7.3)') &
          ' ION: ssca at band edges        = ', &
          ion_dust_ssca_at(elo), ion_dust_ssca_at(ehi), &
          ',  g = ', ion_dust_g_at(elo), ion_dust_g_at(ehi)
    end if
  end subroutine ion_dust_read

  !--- Dust cross sections at a photon's exact energy E [eV] from the retained
  !--- extinction table: log-log in C_ext, linear in albedo and <cos>, over
  !--- wavelength.  Called for every packet (ion_packet_mod) and for the
  !--- infrared emission bands (dust_emis_transport_mod).
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
       if (trim(par%dust_density_law) == 'laursen09_live') then
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
