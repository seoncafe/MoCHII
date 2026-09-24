module photoelectron_degradation_mod
!---------------------------------------------------------------------------
! MoCHII: how the excess energy E0 = h nu - E_th of a fast photoelectron is
! degraded in a partly ionized H/He gas -- the share that ends as heat and
! the share that ionizes each neutral target (par%use_sec_ion).
!
! Two published determinations, each for what it provides:
!
! * Shull & van Steenberg (1985, ApJ 298, 268; SvS85), Table 2 and their
!   eqs. (1)-(2), read from the published scan: the ENERGY BUDGET of the
!   secondary-ionization channels against the ionized fraction,
!     phi_HI(x)  = 0.3908 (1 - x^0.4092)^1.7592,
!     phi_HeI(x) = 0.0554 (1 - x^0.4614)^1.6660,
!   the fraction of E0 spent ionizing H I and He I.  The form vanishes at
!   x = 1.  The table footnote reads "These formulae are only accurate for
!   E0 > 100 eV"; their gas has n(He)/n(H) = 0.1 with equal H and He
!   ionization fractions, x = n(H+)/n(H_tot).
!
! * Dalgarno, Yan & Liu (1999, ApJS 125, 237; DYL99), read from the
!   published PDF: the dependence on E0 over 30-1000 eV that SvS85 do not
!   resolve, through the mean energy per ion pair, their eq. (13)
!     W(x, E0) = W0 (1 + C x^a)          (Table 4, "H and He Gas Mixture"),
!   and the heating efficiency, their eq. (14)
!     eta(x, E0) = 1 + (eta0 - 1)/(1 + C x^a)   (Table 7, "H, He" mixture),
!   both stated for x <= 0.1.
!
! Combination (the form EXHALE uses, electron_energy_degradation.f90):
!   f_HI (x,E0) = phi_HI (x) W_H+ (x, 1 keV)/W_H+ (x, E0),
!   f_HeI(x,E0) = phi_HeI(x) W_He+(x, 1 keV)/W_He+(x, E0),
! so the 1 keV limit is exactly SvS85 and a slower photoelectron makes fewer
! ion pairs for its energy (W_He+ is 487 eV at 1 keV but 6740 eV at 30 eV:
! 30 eV is barely above the 24.59 eV He threshold).  The budget f_ion =
! f_HI + f_HeI is divided among the targets actually present, on SvS85's
! reference composition n(He I)/n(H I) = 0.1:
!   w_tot = 0.1 f_HI n_HI + f_HeI n_HeI,
!   c_HI  = f_ion 0.1 f_HI / w_tot,   c_HeI = f_ion f_HeI / w_tot   [cm^3],
! so n_HI c_HI + n_HeI c_HeI = f_ion at every composition, and at the
! reference one c_HI n_HI = f_HI, c_HeI n_HeI = f_HeI exactly.  The heat
! share is
!   f_heat(x,E0) = 1 - [1 - eta(x,E0)] (1 - x),  clamped to [0, 1 - f_ion],
! where the (1 - x) closure (EXHALE's, not either paper's) sends the
! unthermalized share to zero in a fully ionized gas, which eq. (14) alone
! does not (it leaves (1 - eta0)/(1 + C) = 10% at 1 keV and x = 1).  The
! rest, 1 - f_heat - f_ion, is collisional excitation, taken to leave as
! line radiation.
!
! x is n_e/(n_H + n_He): the Coulomb losses that compete with ionization
! are to free electrons, so the metal electrons count; it equals SvS85's
! n(H+)/n(H_tot) under their equal-ionization assumption.  Clamped to [0,1].
!
! VALIDITY AND WHAT IS APPROXIMATE.
!  * E0 below 30 eV (the bottom of DYL99) thermalizes fully: no tabulation
!    covers 13.6-30 eV, where a photoelectron can still ionize H, so the
!    heat there is an upper bound and the partition steps at 30 eV.  Above
!    1 keV the coefficients are held at their 1 keV values.
!  * DYL99's fits are stated for x <= 0.1 and SvS85's for 1e-4 < x < 1;
!    above 0.1 the E0 dependence and eq. (14) are extrapolated, as in
!    EXHALE (which measured eq. (14) with the closure against SvS85's own
!    3 keV Monte Carlo over x > 0.1 to within 5.6%).
!  * Both calculations have n(He)/n(H) = 0.1; the redistribution corrects
!    the branching at other compositions, not the total budget.
!  * Metal photoelectrons are not partitioned (their energy is all heat) and
!    metals are not targets: they are trace absorbers.
!---------------------------------------------------------------------------
  use define, only : wp
  implicit none
  private

  public :: N_EDEP, E_SECONDARY, photoelectron_energy_weights
  public :: photoelectron_partition, ionized_fraction_svs
  public :: svs85_fion_HI, svs85_fion_HeI

  !--- DYL99 energy nodes [eV] and the threshold below which E0 is all heat.
  integer, parameter :: N_EDEP = 6
  real(kind=wp), parameter :: E_NODE(N_EDEP) = &
       [30.0_wp, 50.0_wp, 100.0_wp, 200.0_wp, 500.0_wp, 1000.0_wp]
  real(kind=wp), parameter :: E_SECONDARY = 30.0_wp
  !--- SvS85 reference neutral ratio n(He I)/n(H I)
  real(kind=wp), parameter :: HEH_SVS = 0.1_wp

  !--- DYL99 Table 4, H and He gas mixture: (W0 [eV], a, C) at each node.
  real(kind=wp), parameter :: W_HP(3,N_EDEP) = reshape([ &
       60.0_wp, 1.06_wp,  221.0_wp, &
       48.9_wp, 1.01_wp,   87.9_wp, &
       42.6_wp, 0.964_wp,  41.0_wp, &
       40.3_wp, 0.928_wp,  24.5_wp, &
       39.8_wp, 0.890_wp,  15.5_wp, &
       39.8_wp, 0.866_wp,  12.2_wp], [3,N_EDEP])
  real(kind=wp), parameter :: W_HEP(3,N_EDEP) = reshape([ &
       6740.0_wp, 1.08_wp,  522.0_wp, &
       2150.0_wp, 1.07_wp,   98.4_wp, &
       1030.0_wp, 1.03_wp,   40.9_wp, &
        700.0_wp, 1.02_wp,   24.1_wp, &
        540.0_wp, 1.00_wp,   15.6_wp, &
        487.0_wp, 0.994_wp,  12.5_wp], [3,N_EDEP])
  !--- DYL99 Table 7, H, He two-gas mixture: (eta0, a, C) at each node.
  real(kind=wp), parameter :: ETA_HHE(3,N_EDEP) = reshape([ &
       0.259_wp, 0.984_wp, 181.0_wp, &
       0.186_wp, 0.834_wp,  53.4_wp, &
       0.148_wp, 0.781_wp,  25.8_wp, &
       0.132_wp, 0.743_wp,  15.8_wp, &
       0.122_wp, 0.703_wp,  10.1_wp, &
       0.117_wp, 0.678_wp,   7.95_wp], [3,N_EDEP])

contains

  !=========================================================================
  ! Interpolation weights of a photoelectron energy E0 [eV] on the DYL99
  ! nodes: linear in log10(E0), clamped at 1 keV.  The weights sum to one and
  ! do not depend on the gas state, so a path estimator can accumulate the
  ! deposited energy node by node and the state-dependent partition is
  ! applied afterwards exactly.  E0 <= E_SECONDARY returns all zeros (the
  ! energy is heat, outside the partition).
  !=========================================================================
  pure subroutine photoelectron_energy_weights(E0, w)
    real(kind=wp), intent(in)  :: E0
    real(kind=wp), intent(out) :: w(N_EDEP)
    real(kind=wp) :: le, t
    integer :: k
    w = 0.0_wp
    if (E0 <= E_SECONDARY) return
    if (E0 >= E_NODE(N_EDEP)) then
       w(N_EDEP) = 1.0_wp
       return
    end if
    le = log10(E0)
    k = 1
    do while (k < N_EDEP - 1 .and. E0 > E_NODE(k+1))
       k = k + 1
    end do
    t = (le - log10(E_NODE(k)))/(log10(E_NODE(k+1)) - log10(E_NODE(k)))
    w(k) = 1.0_wp - t
    w(k+1) = t
  end subroutine photoelectron_energy_weights

  !=========================================================================
  ! x = n_e/(n_H + n_He), clamped to [0, 1].
  !=========================================================================
  pure real(kind=wp) function ionized_fraction_svs(ne, nH, nHe) result(x)
    real(kind=wp), intent(in) :: ne, nH, nHe
    if (nH + nHe > 0.0_wp) then
       x = min(max(ne/(nH + nHe), 0.0_wp), 1.0_wp)
    else
       x = 1.0_wp
    end if
  end function ionized_fraction_svs

  elemental real(kind=wp) function svs85_fion_HI(x) result(f)
    real(kind=wp), intent(in) :: x
    f = 0.3908_wp*(1.0_wp - x**0.4092_wp)**1.7592_wp
  end function svs85_fion_HI

  elemental real(kind=wp) function svs85_fion_HeI(x) result(f)
    real(kind=wp), intent(in) :: x
    f = 0.0554_wp*(1.0_wp - x**0.4614_wp)**1.6660_wp
  end function svs85_fion_HeI

  !=========================================================================
  ! The partition at each node for the state (x, n_HI, n_HeI):
  !   f_heat(k) heat share of E0,
  !   c_HI(k), c_HeI(k) [cm^3] secondary-ionization energy share per target
  !   particle (n_HI c_HI + n_HeI c_HeI = f_ion(k)).
  !=========================================================================
  pure subroutine photoelectron_partition(x, nHI, nHeI, f_heat, c_HI, c_HeI)
    real(kind=wp), intent(in)  :: x, nHI, nHeI
    real(kind=wp), intent(out) :: f_heat(N_EDEP), c_HI(N_EDEP), c_HeI(N_EDEP)
    real(kind=wp) :: xc, phiH, phiHe, WtopH, WtopHe, fH, fHe, fion, wtot, eta
    integer :: k
    xc     = min(max(x, 0.0_wp), 1.0_wp)
    phiH   = svs85_fion_HI(xc)
    phiHe  = svs85_fion_HeI(xc)
    WtopH  = W_HP(1,N_EDEP)*(1.0_wp + W_HP(3,N_EDEP)*xc**W_HP(2,N_EDEP))
    WtopHe = W_HEP(1,N_EDEP)*(1.0_wp + W_HEP(3,N_EDEP)*xc**W_HEP(2,N_EDEP))
    do k = 1, N_EDEP
       fH   = phiH *WtopH /(W_HP(1,k) *(1.0_wp + W_HP(3,k) *xc**W_HP(2,k)))
       fHe  = phiHe*WtopHe/(W_HEP(1,k)*(1.0_wp + W_HEP(3,k)*xc**W_HEP(2,k)))
       fion = fH + fHe
       wtot = HEH_SVS*fH*nHI + fHe*nHeI
       if (wtot > 0.0_wp) then
          c_HI(k)  = fion*HEH_SVS*fH/wtot
          c_HeI(k) = fion*fHe/wtot
       else
          c_HI(k)  = 0.0_wp
          c_HeI(k) = 0.0_wp
          fion     = 0.0_wp          ! no neutral target: nothing ionized
       end if
       eta = 1.0_wp + (ETA_HHE(1,k) - 1.0_wp) &
                      /(1.0_wp + ETA_HHE(3,k)*xc**ETA_HHE(2,k))
       f_heat(k) = min(max(1.0_wp - (1.0_wp - eta)*(1.0_wp - xc), 0.0_wp), &
                       max(1.0_wp - fion, 0.0_wp))
    end do
  end subroutine photoelectron_partition

end module photoelectron_degradation_mod
