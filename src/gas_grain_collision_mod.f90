module gas_grain_collision_mod
!---------------------------------------------------------------------------
! MoCHII: energy exchanged in collisions of neutral atoms with grains
! (par%gas_grain_energy_exchange).
!
! Each neutral atom of species i that strikes a grain of radius a brings in a
! mean kinetic energy 2 k T (the Maxwellian flux average) and leaves with
! 2 k [T - alpha_i (T - T_d)], alpha_i the thermal accommodation coefficient.
! Summed over the size bins of the run's grains (grain_size_bins_mod)
! the gas cooling rate is
!   Lambda_gd = Sum_bins n_gr pi a^2 Sum_i n_i (8 k T / pi m_i)^(1/2)
!               alpha_i 2 k (T - T_d)                      [erg cm^-3 s^-1],
! Hollenbach & McKee (1979, ApJS 41, 555; HM79) eq. 3.25 with the species
! written out; it is negative (gas heating) where T_d > T.
!
! Accommodation, HM79 eqs. 3.23-3.24 (read):
!   H:  alpha = (1 - alpha0) exp{-[2 k (T_d + T)/D]^(1/2)} + alpha0,
!       D = 500 k (the adsorption energy), alpha0 = 2 m_H m_gr/(m_H + m_gr)^2,
!       the high-temperature hard-sphere limit;
!   He: the hard-sphere value alpha0 with m_He.
! m_gr, the mass of the surface atom struck: 12 amu for carbonaceous grains;
! 172/7 = 24.6 amu for silicates and astrodust, the mean atomic mass of the
! MgFeSiO4 stoichiometry Weingartner & Draine (2001, ApJ 548, 296, sec. 2.4)
! adopt.
!
! Only neutral atoms are here.  Electrons and ions stick to the grain (the
! ions are neutralized there); the energy they bring in, with the Coulomb
! focusing of the charged grain, and the 2 k T_d the neutralized ions take
! back are counted with the grain charge in grain_photoelectric_mod.  H2 is
! not a species of MoCHII.
!
! T_d: the leaf's equilibrium mixture temperature (dust_temp_mod), applied to
! every size.  Grains small enough to be stochastically heated spend most of
! the time below the equilibrium temperature, so (T - T_d) is underestimated
! for them by (T_eq - <T_d>); in ionized gas (T ~ 1e4 K, T_d ~ 1e2 K) that is
! below 1%, in a cold neutral zone (T ~ 1e2 K) it can reach ~10%.
!
! Default on with dust: gas and grains at different temperatures exchange
! energy whenever atoms strike grains, which happens wherever both exist.
!---------------------------------------------------------------------------
  use define
  use grain_size_bins_mod, only : gbin, n_gbin, CARBONACEOUS
  implicit none
  private

  public :: gas_grain_neutral_cooling, thermal_accommodation_H, &
            thermal_accommodation_He

  real(kind=wp), parameter :: AMU = 1.66053907e-24_wp
  real(kind=wp), parameter :: M_H  = 1.00794_wp, M_HE = 4.002602_wp      ! [amu]
  real(kind=wp), parameter :: M_SURFACE_CARB = 12.011_wp
  real(kind=wp), parameter :: M_SURFACE_SIL  = 172.24_wp/7.0_wp
  real(kind=wp), parameter :: D_ADS_K = 500.0_wp

contains

  !--- HM79 eq. 3.23-3.24 for H on a surface of atomic mass mgr [amu]
  pure real(kind=wp) function thermal_accommodation_H(T, Td, mgr) result(alpha)
    real(kind=wp), intent(in) :: T, Td, mgr
    real(kind=wp) :: a0
    a0 = 2.0_wp*M_H*mgr/(M_H + mgr)**2
    alpha = (1.0_wp - a0)*exp(-sqrt(2.0_wp*(Td + T)/D_ADS_K)) + a0
  end function thermal_accommodation_H

  !--- HM79: the hard-sphere value for He
  pure real(kind=wp) function thermal_accommodation_He(mgr) result(alpha)
    real(kind=wp), intent(in) :: mgr
    alpha = 2.0_wp*M_HE*mgr/(M_HE + mgr)**2
  end function thermal_accommodation_He

  !=========================================================================
  ! Lambda_gd [erg cm^-3 s^-1] of H0 (nHI) and He0 (nHeI) on the run's
  ! grains; dens_other / dens_pah = n_H times the leaf's abundance scale of
  ! the other / PAH bins (grain_size_bins_mod::grain_leaf_scale).
  !=========================================================================
  real(kind=wp) function gas_grain_neutral_cooling(T, Td, nHI, nHeI, dens_other, &
                                                   dens_pah) result(cool)
    real(kind=wp), intent(in) :: T, Td, nHI, nHeI, dens_other, dens_pah
    real(kind=wp) :: vH, vHe, flux_c, flux_s, area_c, area_s, d, ac, as
    integer :: ib

    cool = 0.0_wp
    if (n_gbin == 0) return
    area_c = 0.0_wp;  area_s = 0.0_wp
    do ib = 1, n_gbin
       d = merge(dens_pah, dens_other, gbin(ib)%pah)
       if (gbin(ib)%material == CARBONACEOUS) then
          area_c = area_c + d*gbin(ib)%n_per_H*pi*gbin(ib)%a_cm**2
       else
          area_s = area_s + d*gbin(ib)%n_per_H*pi*gbin(ib)%a_cm**2
       end if
    end do
    vH  = sqrt(8.0_wp*kboltz_cgs*T/(pi*M_H*AMU))
    vHe = sqrt(8.0_wp*kboltz_cgs*T/(pi*M_HE*AMU))
    ac = thermal_accommodation_H(T, Td, M_SURFACE_CARB)
    as = thermal_accommodation_H(T, Td, M_SURFACE_SIL)
    flux_c = nHI*vH*ac + nHeI*vHe*thermal_accommodation_He(M_SURFACE_CARB)
    flux_s = nHI*vH*as + nHeI*vHe*thermal_accommodation_He(M_SURFACE_SIL)
    cool = (area_c*flux_c + area_s*flux_s)*2.0_wp*kboltz_cgs*(T - Td)
  end function gas_grain_neutral_cooling

end module gas_grain_collision_mod
