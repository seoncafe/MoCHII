module photo_xsec
!---------------------------------------------------------------------------
! MoCHII: photoionization cross sections.
!
! Analytic single-shell fits from Verner, Ferland, Korista & Yakovlev 1996,
! ApJ 465, 487 (VFKY96), Eqs. (1)-(4).  Parameters are the PH2 table of the
! authors' phfit2.f (verified against ~/RT_Codes/Verner_Fortran_sub/phfit2.f
! on 2026-07-11):
!   H I   (1,1): E_0=0.4298, sigma_0=5.475e4, y_a=32.88, P=2.963, 0,0,0
!   He II (2,1): E_0=1.720,  sigma_0=1.369e4, y_a=32.88, P=2.963, 0,0,0
!   He I  (2,2): E_0=13.61,  sigma_0=949.2,   y_a=1.469, P=3.188,
!                y_w=2.039, y_0=0.4434, y_1=2.136
! Thresholds from the PH1 table (13.60, 54.42, 24.59 eV); define.f90 carries
! them as eth_HI/eth_HeI/eth_HeII.  sigma_0 is in Mb; results returned in
! cm^2.  The same fit form serves every metal ion (registry).
!
! Cross-checked against EXHALE cross_sec.f90 and its comparison script
! docs/compare_photoion_cross_sections.py.
!---------------------------------------------------------------------------
  use define, only : wp, eth_HI, eth_HeI, eth_HeII
  implicit none
  private

  public :: sigma_vfky96, sigma_vy95, sigma_HI, sigma_HeI, sigma_HeII

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
  elemental real(kind=wp) function sigma_HI(E) result(sigma)
    real(kind=wp), intent(in) :: E    ! photon energy [eV]
    sigma = sigma_vfky96(E, eth_HI, 4.298e-1_wp, 5.475e4_wp, 3.288e1_wp, &
                         2.963_wp, 0.0_wp, 0.0_wp, 0.0_wp)
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
    sigma = sigma_vfky96(E, eth_HeII, 1.720_wp, 1.369e4_wp, 3.288e1_wp, &
                         2.963_wp, 0.0_wp, 0.0_wp, 0.0_wp)
  end function sigma_HeII

end module photo_xsec
