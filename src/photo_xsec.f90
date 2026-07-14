module photo_xsec
!---------------------------------------------------------------------------
! MoCHII: photoionization cross sections.
!
! The one-electron ions H I (Z=1) and He II (Z=2) use the EXACT ground-state
! hydrogenic cross section (sigma_hydrogenic); the analytic VFKY96 fit is kept
! for He I (two electrons) and every metal ion (registry).
!
! Exact hydrogenic (Osterbrock & Ferland 2006, Eq. 2.4; same form as EXHALE
! cross_sec.f90):
!   eps   = sqrt(E/E_th - 1)
!   sigma = sigma_0 (E_th/E)^4 exp(4 - 4 arctan(eps)/eps) / (1 - exp(-2 pi/eps))
! with sigma_0 = 6.30e-18/Z^2 cm^2 the threshold value.  E_th is the true
! ionization potential from define.f90 (eth_HI = 13.598, eth_HeII = 54.416 eV),
! so the edge stays consistent with the band grid / ion_align_edges.
!
! VFKY96: analytic single-shell fit from Verner, Ferland, Korista & Yakovlev
! 1996, ApJ 465, 487, Eqs. (1)-(4).  Parameters are the PH2 table of the
! authors' phfit2.f (verified against ~/RT_Codes/Verner_Fortran_sub/phfit2.f
! on 2026-07-11):
!   He I  (2,2): E_0=13.61, sigma_0=949.2, y_a=1.469, P=3.188,
!                y_w=2.039, y_0=0.4434, y_1=2.136
! Threshold eth_HeI = 24.587 eV.  sigma_0 is in Mb; results returned in cm^2.
! The same fit form serves every metal ion.
!
! Cross-checked against EXHALE cross_sec.f90 and its comparison script
! docs/compare_photoion_cross_sections.py.
!---------------------------------------------------------------------------
  use define, only : wp, pi, eth_HI, eth_HeI, eth_HeII
  implicit none
  private

  public :: sigma_vfky96, sigma_vy95, sigma_hydrogenic, &
            sigma_HI, sigma_HeI, sigma_HeII

  real(kind=wp), parameter :: Mb2cm2 = 1.0e-18_wp

contains

  !=========================================================================
  ! Full VFKY96 Eq. (1) with the (y_0, y_1) offset/asymptote parameters:
  !   x = E/E_0 - y_0,  z = sqrt(x^2 + y_1^2),  Q = 5.5 - 0.5*P
  !   sigma = sigma_0 * ((x-1)^2 + y_w^2) * z^(-Q) * (1 + sqrt(z/y_a))^(-P)
  ! E and E_th in eV, sigma_0 in Mb; returns cm^2 (0 below threshold).
  !=========================================================================
  elemental real(kind=wp) function sigma_vfky96(E, E_th, E_0, s_0, y_a, P, &
                                                y_w, y_0, y_1) result(sigma)
    real(kind=wp), intent(in) :: E, E_th, E_0, s_0, y_a, P, y_w, y_0, y_1
    real(kind=wp) :: x, z, Q, Fy
    if (E < E_th) then
       sigma = 0.0_wp
       return
    end if
    x  = E/E_0 - y_0
    z  = sqrt(x*x + y_1*y_1)
    Q  = 5.5_wp - 0.5_wp*P
    Fy = ((x - 1.0_wp)**2 + y_w**2) * z**(-Q) * (1.0_wp + sqrt(z/y_a))**(-P)
    sigma = s_0 * Fy * Mb2cm2
  end function sigma_vfky96

  !=========================================================================
  ! Verner & Yakovlev (1995) subshell fit — the outer-shell cross section
  ! for the elements ABSENT from the VFKY96 outer-shell table (P, Cl, K;
  ! the same hybrid as MOCASSIN phFitEl):
  !   y = E/E_0,  Q = 5.5 + l - 0.5 P
  !   sigma = sigma_0 ((y-1)^2 + y_w^2) y^(-Q) (1 + sqrt(y/y_a))^(-P)
  ! with l the subshell orbital quantum number.
  !=========================================================================
  elemental real(kind=wp) function sigma_vy95(E, E_th, E_0, s_0, y_a, P, &
                                              y_w, l) result(sigma)
    real(kind=wp), intent(in) :: E, E_th, E_0, s_0, y_a, P, y_w, l
    real(kind=wp) :: y, Q, Fy
    if (E < E_th) then
       sigma = 0.0_wp
       return
    end if
    y  = E/E_0
    Q  = 5.5_wp + l - 0.5_wp*P
    Fy = ((y - 1.0_wp)**2 + y_w**2) * y**(-Q) &
         * (1.0_wp + sqrt(y/y_a))**(-P)
    sigma = s_0 * Fy * Mb2cm2
  end function sigma_vy95

  !=========================================================================
  ! Exact ground-state hydrogenic photoionization cross section for a
  ! one-electron ion of nuclear charge Z (Osterbrock & Ferland 2006, Eq. 2.4):
  !   eps   = sqrt(E/E_th - 1)
  !   sigma = sigma_0 (E_th/E)^4 exp(4 - 4 arctan(eps)/eps) / (1 - exp(-2 pi/eps))
  ! sigma_0 = 6.30e-18/Z^2 cm^2 is the threshold value; E, E_th in eV.
  ! Returns 0 below threshold and the analytic sigma_0 limit exactly at E=E_th.
  !=========================================================================
  elemental real(kind=wp) function sigma_hydrogenic(E, E_th, Z) result(sigma)
    real(kind=wp), intent(in) :: E, E_th, Z
    real(kind=wp) :: eps, s0
    real(kind=wp), parameter :: sig0 = 6.30e-18_wp     ! Z=1 threshold [cm^2]
    s0 = sig0 / (Z*Z)
    if (E < E_th) then
       sigma = 0.0_wp
    else if (E > E_th) then
       eps   = sqrt(E/E_th - 1.0_wp)
       sigma = s0 * (E_th/E)**4 &
                  * exp(4.0_wp - 4.0_wp*atan(eps)/eps) &
                  / (1.0_wp - exp(-2.0_wp*pi/eps))
    else
       sigma = s0                                       ! E == E_th (eps -> 0)
    end if
  end function sigma_hydrogenic

  !=========================================================================
  elemental real(kind=wp) function sigma_HI(E) result(sigma)
    real(kind=wp), intent(in) :: E    ! photon energy [eV]
    sigma = sigma_hydrogenic(E, eth_HI, 1.0_wp)
  end function sigma_HI

  !=========================================================================
  elemental real(kind=wp) function sigma_HeI(E) result(sigma)
    real(kind=wp), intent(in) :: E
    sigma = sigma_vfky96(E, eth_HeI, 1.361e1_wp, 9.492e2_wp, 1.469_wp, &
                         3.188_wp, 2.039_wp, 4.434e-1_wp, 2.136_wp)
  end function sigma_HeI

  !=========================================================================
  elemental real(kind=wp) function sigma_HeII(E) result(sigma)
    real(kind=wp), intent(in) :: E
    sigma = sigma_hydrogenic(E, eth_HeII, 2.0_wp)
  end function sigma_HeII

end module photo_xsec
