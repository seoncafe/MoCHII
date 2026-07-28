module ion_packet_mod
!---------------------------------------------------------------------------
! Packet-local cross-section cache and grouped dynamic-opacity shadow.
!
! The cache is built once on entry to ionizing transport.  During the current
! grouped phase photon%energy_eV equals ion_e(inu), but every dynamic term is
! evaluated through the packet cache and the leaf state.  Transport still
! returns and uses the authoritative precomputed kap_ion value; the dynamic
! value is validation-only until the continuous-energy cutover.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: build_ion_packet_physics, ion_packet_opacity
  public :: packet_opacity_shadow_reset, packet_opacity_shadow_reduce
  public :: packet_opacity_shadow_report

  real(kind=wp) :: opacity_max_abs = 0.0_wp
  real(kind=wp) :: opacity_scale = 0.0_wp
  integer(kind=int64) :: opacity_checks = 0_int64

contains

  subroutine build_ion_packet_physics(photon)
    use ion_band_mod,   only : ion_e
    use photo_xsec,     only : sigma_HI, sigma_HeI, sigma_HeII
    use gas_opacity_mod,only : ion_dust_sabs, ion_dust_ssca, ion_dust_g
    implicit none
    type(photon_type), intent(inout) :: photon
    real(kind=wp) :: energy

    energy = photon%energy_eV
    photon%ionphys%energy_eV = energy
    photon%ionphys%sigma_HI = sigma_HI(energy)
    photon%ionphys%sigma_HeI = sigma_HeI(energy)
    photon%ionphys%sigma_HeII = sigma_HeII(energy)
    photon%ionphys%n_metal = 0
    photon%ionphys%sigma_metal = 0.0_wp
    if (par%use_metals) then
       block
         use species_mod, only : species_packet_cross_sections
         call species_packet_cross_sections(energy, &
              photon%ionphys%sigma_metal, photon%ionphys%n_metal)
       end block
    end if
    photon%ionphys%dust_abs = 0.0_wp
    photon%ionphys%dust_sca = 0.0_wp
    photon%ionphys%dust_g = 0.0_wp
    if (par%ion_add_dust) then
       !--- Phase-16.4 grouped shadow: the table has so far been retained only
       !--- on the band grid, so these are the values at ion_e(inu).  The
       !--- continuous sampler phase will retain/interpolate the source table
       !--- at energy_eV instead.
       photon%ionphys%dust_abs = ion_dust_sabs(photon%inu)
       photon%ionphys%dust_sca = ion_dust_ssca(photon%inu)
       photon%ionphys%dust_g = ion_dust_g(photon%inu)
    end if

    if (trim(par%ion_energy_mode) == 'grouped') then
       if (energy /= ion_e(photon%inu)) then
          write(*,'(a,i0,2(a,es24.16))') &
             'ERROR: grouped packet energy mismatch in bin ', photon%inu, &
             ': packet=', energy, ', center=', ion_e(photon%inu)
          error stop 1
       end if
    end if
  end subroutine build_ion_packet_physics

  real(kind=wp) function ion_packet_opacity(photon, il) result(kap)
    use octree_mod,    only : amr_grid
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII
    use gas_opacity_mod, only : kap_ion
    implicit none
    type(photon_type), intent(in) :: photon
    integer,           intent(in) :: il
    real(kind=wp) :: dynamic, metal

    kap = kap_ion(photon%inu,il)
    if (trim(par%ion_energy_mode) /= 'continuous' .and. &
        .not. par%ion_shadow_rates) return

    metal = 0.0_wp
    if (par%use_metals .and. par%ion_metal_abs) then
       block
         use species_mod, only : species_cached_opacity
         metal = species_cached_opacity(photon%ionphys%sigma_metal, &
                 photon%ionphys%n_metal, il)
       end block
    end if
    dynamic = gas_nH(il)*( &
         gas_xHI(il)*photon%ionphys%sigma_HI &
       + par%He_abund*(gas_xHeI(il)*photon%ionphys%sigma_HeI &
       + gas_xHeII(il)*photon%ionphys%sigma_HeII) + metal) &
       * par%distance2cm
    if (par%ion_add_dust) then
       dynamic = dynamic + amr_grid%rhokap(il)*photon%ionphys%dust_abs
       if (par%ion_dust_scatter) &
          dynamic = dynamic + amr_grid%rhokap(il)*photon%ionphys%dust_sca
    end if
    if (trim(par%ion_energy_mode) == 'continuous') then
       kap = dynamic
       return
    end if
    opacity_max_abs = max(opacity_max_abs, abs(dynamic-kap))
    opacity_scale = max(opacity_scale, abs(kap))
    opacity_checks = opacity_checks + 1_int64
  end function ion_packet_opacity

  subroutine packet_opacity_shadow_reset()
    opacity_max_abs = 0.0_wp
    opacity_scale = 0.0_wp
    opacity_checks = 0_int64
  end subroutine packet_opacity_shadow_reset

  subroutine packet_opacity_shadow_reduce()
    use mpi
    implicit none
    real(kind=wp) :: maxima(2)
    integer(kind=int64) :: count
    integer :: ierr
    if (.not. par%ion_shadow_rates) return
    maxima = [opacity_max_abs, opacity_scale]
    call MPI_ALLREDUCE(MPI_IN_PLACE, maxima, 2, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    count = opacity_checks
    call MPI_ALLREDUCE(count, opacity_checks, 1, MPI_INTEGER8, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    opacity_max_abs = maxima(1)
    opacity_scale = maxima(2)
  end subroutine packet_opacity_shadow_reduce

  subroutine packet_opacity_shadow_report(rtol)
    use mpi
    implicit none
    real(kind=wp), intent(in) :: rtol
    real(kind=wp) :: rel
    integer :: ierr
    logical :: passed
    if (.not. par%ion_shadow_rates) return
    if (opacity_scale > 0.0_wp) then
       rel = opacity_max_abs/opacity_scale
    else
       rel = opacity_max_abs
    end if
    passed = opacity_checks > 0 .and. rel <= rtol
    if (mpar%p_rank == 0) then
       write(*,'(a,i0,a,es11.3)') ' ION: dynamic-opacity shadow checks=', &
          opacity_checks, ', relative error=', rel
       write(*,'(a,a)') ' ION: grouped dynamic-opacity shadow: ', &
                        merge('PASS', 'FAIL', passed)
    end if
    if (.not. passed) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  end subroutine packet_opacity_shadow_report

end module ion_packet_mod
