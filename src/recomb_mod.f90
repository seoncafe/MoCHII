module recomb_mod
!---------------------------------------------------------------------------
! MoCHII: recombination and collisional-ionization rate coefficients (G1).
!
! Recombination — Hui & Gnedin (1997, MNRAS 292, 27) analytic fits, case A
! and case B, for H II, He II, He III, as functions of T [K]; accuracy
! better than ~0.7% over 1-1e9 K.  He II dielectronic recombination is
! negligible at HII-region temperatures and is deferred to G2, where the
! fitting pipeline supplies Badnell (2006) RR + DR coefficients for the
! metal cascade anyway (docs/PLAN.md section 2.3).  At T = 1e4 K:
! alpha_B(H II) = 2.594e-13 cm^3 s^-1 (canonical 2.59e-13).
!
! Collisional ionization — Voronov (1997, ADNDT 65, 1) fits for H I, He I,
! He II:  k(T) = A (1 + P sqrt(U)) U^K exp(-U) / (X + U),  U = dE/T_eV.
! Negligible next to photoionization at ~1e4 K but kept in the balance
! (docs/PLAN.md section 2.5); G2 reuses them.
!---------------------------------------------------------------------------
  use define, only : wp, kboltz_cgs, ev2erg
  implicit none
  private

  public :: alphaA_HII, alphaB_HII, alphaA_HeII, alphaB_HeII
  public :: alphaA_HeIII, alphaB_HeIII
  public :: ci_HI, ci_HeI, ci_HeII

  !--- ionization-threshold temperatures T_TR = dE/k [K] (Hui & Gnedin).
  real(kind=wp), parameter :: T_TR_HI   = 157807.0_wp
  real(kind=wp), parameter :: T_TR_HeI  = 285335.0_wp
  real(kind=wp), parameter :: T_TR_HeII = 631515.0_wp

contains

  !=========================================================================
  ! Hui & Gnedin (1997) recombination fits; lambda = 2 T_TR / T.
  !=========================================================================
  elemental real(kind=wp) function alphaA_HII(T) result(a)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*T_TR_HI/T
    a = 1.269e-13_wp * lam**1.503_wp / (1.0_wp + (lam/0.522_wp)**0.470_wp)**1.923_wp
  end function alphaA_HII

  elemental real(kind=wp) function alphaB_HII(T) result(a)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*T_TR_HI/T
    a = 2.753e-14_wp * lam**1.500_wp / (1.0_wp + (lam/2.740_wp)**0.407_wp)**2.242_wp
  end function alphaB_HII

  elemental real(kind=wp) function alphaA_HeII(T) result(a)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*T_TR_HeI/T
    a = 3.000e-14_wp * lam**0.654_wp
  end function alphaA_HeII

  elemental real(kind=wp) function alphaB_HeII(T) result(a)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*T_TR_HeI/T
    a = 1.260e-14_wp * lam**0.750_wp
  end function alphaB_HeII

  !--- He III is hydrogenic (Z=2): the H fit evaluated on lambda(T_TR_HeII),
  !--- scaled by Z (alpha_Z(T) = Z alpha_H(T/Z^2)).
  elemental real(kind=wp) function alphaA_HeIII(T) result(a)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*T_TR_HeII/T
    a = 2.0_wp * 1.269e-13_wp * lam**1.503_wp / (1.0_wp + (lam/0.522_wp)**0.470_wp)**1.923_wp
  end function alphaA_HeIII

  elemental real(kind=wp) function alphaB_HeIII(T) result(a)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*T_TR_HeII/T
    a = 2.0_wp * 2.753e-14_wp * lam**1.500_wp / (1.0_wp + (lam/2.740_wp)**0.407_wp)**2.242_wp
  end function alphaB_HeIII

  !=========================================================================
  ! Voronov (1997) collisional-ionization rate coefficients [cm^3 s^-1].
  !=========================================================================
  elemental real(kind=wp) function voronov(T, dE_eV, P, A, X, Kexp) result(rate)
    real(kind=wp), intent(in) :: T, dE_eV, P, A, X, Kexp
    real(kind=wp) :: U
    U = dE_eV * ev2erg / (kboltz_cgs * T)
    rate = A * (1.0_wp + P*sqrt(U)) * U**Kexp * exp(-U) / (X + U)
  end function voronov

  elemental real(kind=wp) function ci_HI(T) result(k)
    real(kind=wp), intent(in) :: T
    k = voronov(T, 13.6_wp, 0.0_wp, 2.91e-8_wp, 0.232_wp, 0.39_wp)
  end function ci_HI

  elemental real(kind=wp) function ci_HeI(T) result(k)
    real(kind=wp), intent(in) :: T
    k = voronov(T, 24.6_wp, 0.0_wp, 1.75e-8_wp, 0.180_wp, 0.35_wp)
  end function ci_HeI

  elemental real(kind=wp) function ci_HeII(T) result(k)
    real(kind=wp), intent(in) :: T
    k = voronov(T, 54.4_wp, 1.0_wp, 2.05e-9_wp, 0.265_wp, 0.25_wp)
  end function ci_HeII

end module recomb_mod
