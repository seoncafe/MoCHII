module ion_packet_mod
!---------------------------------------------------------------------------
! Packet-local cross sections and the absorption coefficient at a packet's
! exact photon energy.
!
! The cache is built once on entry to ionizing transport: every cross section
! (H I, He I, He II, the active metal stages, dust absorption/scattering and
! the HG asymmetry) is evaluated at photon%energy_eV, the energy the source or
! diffuse sampler actually drew.  ion_packet_opacity then combines that cache
! with the leaf state, so the extinction a packet sees is the physical
! sigma(E) n, not a bin representative.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: build_ion_packet_physics, ion_packet_opacity

contains

  subroutine build_ion_packet_physics(photon)
    use photo_xsec,     only : sigma_HI, sigma_HeI, sigma_HeII
    use gas_opacity_mod,only : ion_dust_sabs_at, ion_dust_ssca_at, ion_dust_g_at
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
       !--- Dust cross sections are evaluated at each packet's exact energy,
       !--- from the same retained extinction table the band-edge diagnostics
       !--- are tabulated from (gas_opacity_mod).
       photon%ionphys%dust_abs = ion_dust_sabs_at(energy)
       photon%ionphys%dust_sca = ion_dust_ssca_at(energy)
       photon%ionphys%dust_g = ion_dust_g_at(energy)
    end if
  end subroutine build_ion_packet_physics

  real(kind=wp) function ion_packet_opacity(photon, il) result(kap)
    !--- Absorption coefficient per CODE LENGTH seen by this packet in leaf il:
    !---   kap = nH [ x_HI sigma_HI(E) + y_He (x_HeI sigma_HeI(E)
    !---              + x_HeII sigma_HeII(E)) + metals ] * distance2cm
    !--- plus the leaf grey dust extinction times the packet's dust absorption
    !--- (and scattering, when interactions are sampled).  Every cross section
    !--- comes from the packet cache at photon%energy_eV.
    use octree_mod,    only : amr_grid
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII
    implicit none
    type(photon_type), intent(in) :: photon
    integer,           intent(in) :: il
    real(kind=wp) :: metal

    metal = 0.0_wp
    if (par%use_metals .and. par%ion_metal_abs) then
       block
         use species_mod, only : species_cached_opacity
         metal = species_cached_opacity(photon%ionphys%sigma_metal, &
                 photon%ionphys%n_metal, il)
       end block
    end if
    kap = gas_nH(il)*( &
         gas_xHI(il)*photon%ionphys%sigma_HI &
       + par%He_abund*(gas_xHeI(il)*photon%ionphys%sigma_HeI &
       + gas_xHeII(il)*photon%ionphys%sigma_HeII) + metal) &
       * par%distance2cm
    if (par%ion_add_dust) then
       kap = kap + amr_grid%rhokap(il)*photon%ionphys%dust_abs
       if (par%ion_dust_scatter) &
          kap = kap + amr_grid%rhokap(il)*photon%ionphys%dust_sca
    end if
  end function ion_packet_opacity

end module ion_packet_mod
