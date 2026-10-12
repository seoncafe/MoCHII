module free_free_mod
!---------------------------------------------------------------------------
! MoCHII: thermal free-free (bremsstrahlung) emission of an ion of net
! charge Z in a Maxwellian electron gas -- the constants and the Gaunt
! factor every free-free term of the code uses: the cooling (cooling_mod,
! species_mod::metal_freefree), the written nebular continuum and the
! transported free-free packets (recomb_continuum_mod).
!
! EMISSIVITY.  4 pi j_nu = FF_SPECTRAL Z^2 g_ff(nu, T, Z) T^(-1/2)
! exp(-h nu/kT) n_e n_i [erg s^-1 cm^-3 Hz^-1], with
!     FF_SPECTRAL = (2^5 pi e^6 / 3 m_e c^3) (2 pi / 3 k m_e)^(1/2)
!                 = 6.841593e-38 erg cm^3 s^-1 Hz^-1 K^(1/2)
! (Rybicki & Lightman 1979 eq. 5.14a; CODATA 2018 constants), and its
! frequency integral, FF_TOTAL Z^2 gbar_ff(T, Z) T^(1/2) n_e n_i with
! FF_TOTAL = FF_SPECTRAL k/h = 1.425557e-27 erg cm^3 s^-1 K^(-1/2).
! The two used to be separate rounded numbers: 6.8391e-38 (MOCASSIN's
! fb_ff, 3.7e-4 below the constant) in the written continuum and 1.42554e-27
! in the cooling, so the spectrum the gas radiated and the energy it lost
! differed by 3.7e-4 of the free-free term.  One constant now gives both.
!
! GAUNT FACTOR (free_free_gaunt).  The ported Hummer (1988) fit (gaunt.f90)
! inside its stated range, 1e-4 <= u = h nu/kT <= 10^1.5,
! 1e-3 <= gamma^2 = Z^2 Ry/kT <= 1e3 and T >= 10^2.5 K; the van Hoof et al.
! (2014) table (gaunt_vh14_mod) outside it, and everywhere with
! par%gaunt_vh14.  This is the rule cooling_mod already applied to the
! thermal average; the spectral coefficient of the written continuum called
! Hummer at every u, where above u = 10^1.5 (h nu > 27 eV at 1e4 K, i.e.
! across the whole ionizing band) the fit returns Ferguson's asymptote
! min(1.05, sqrt(12 gamma^2/u)) rather than a computed Gaunt factor.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: free_free_gaunt
  public :: FF_SPECTRAL, FF_TOTAL

  real(kind=wp), parameter :: FF_SPECTRAL = 6.841593e-38_wp
  real(kind=wp), parameter :: FF_TOTAL    = FF_SPECTRAL*kboltz_cgs/h_planck_cgs
  !--- Rydberg energy [eV] used for gamma^2 and for the Hummer frequency
  !--- argument, as in cooling_mod's thermal average
  real(kind=wp), parameter :: RYD_EV = 13.605693_wp

  interface
     subroutine getGauntFF(z, log10Te, xlf, g, iflag)
        implicit none
        integer, intent(out) :: iflag
        real, intent(in)  :: log10Te, z
        real, dimension(:), intent(in)  :: xlf
        real, dimension(size(xlf)), intent(out) :: g
     end subroutine getGauntFF
  end interface

contains

  !=========================================================================
  ! Free-free Gaunt factor g_ff at net charge Zc, temperature T [K] and
  ! u = h nu / kT (see the module header for the rule).
  !=========================================================================
  real(kind=wp) function free_free_gaunt(Zc, T, u) result(g)
    use gaunt_vh14_mod, only : gaunt_vh14_setup, gauntff_vh14
    use utility,        only : fatal_error
    implicit none
    real(kind=wp), intent(in) :: Zc, T, u
    real(kind=wp) :: ryd_over_kt, gamma2
    real :: xlf(1), gf(1)
    integer :: iflag
    ryd_over_kt = RYD_EV*ev2erg/(kboltz_cgs*T)
    gamma2 = Zc*Zc*ryd_over_kt
    if (par%gaunt_vh14 .or. T < 10.0_wp**2.5_wp .or. gamma2 < 1.0e-3_wp .or. &
        gamma2 > 1.0e3_wp .or. u < 1.0e-4_wp .or. u > 10.0_wp**1.5_wp) then
       if (log10(gamma2) < -6.0_wp .or. log10(gamma2) > 10.0_wp .or. &
           log10(u) < -16.0_wp .or. log10(u) > 13.0_wp) &
          call fatal_error('free-free Gaunt arguments outside van Hoof table')
       call gaunt_vh14_setup()
       g = gauntff_vh14(Zc, T, u/ryd_over_kt)
    else
       xlf(1) = real(log10(u/ryd_over_kt))
       call getGauntFF(real(log10(Zc)), real(log10(T)), xlf, gf, iflag)
       g = real(gf(1), wp)
    end if
  end function free_free_gaunt

end module free_free_mod
