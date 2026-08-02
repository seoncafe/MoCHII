! MoCHII: copied from MoCafe_v2.00/src/physics_amr_mod.f90 (2026-07-11)
module physics_amr_mod
!---------------------------------------------------------------------------
! Dust physics for the AMR grid (dust-only port of LaRT physics_amr_mod.f90).
!
! Only the metallicity->dust relation (Laursen+09) is kept.  For dust_density_law =
! 'laursen09' an explicit xHI column in the file is required.
!---------------------------------------------------------------------------
  use define, only: wp
  implicit none
  private

  public :: laursen09_ndust, dust_opac_norm

  !--- Scale factor that put the INITIAL dust opacity on the requested
  !--- par%taumax (radial pole) or par%tauhomo (volume average).  It is
  !--- solved once, at grid setup, and then held fixed: for dust_density_law =
  !--- 'laursen09_live' every later opacity refill multiplies by it, so the
  !--- target sets the optical depth of the initial configuration and the
  !--- subsequent evolution comes from the ionization-dependent grain
  !--- survival alone.  It is deliberately never re-solved -- pinning the
  !--- realized tau to the target every iteration would defeat the live
  !--- dust model.  Without a target it stays 1.  Kept here rather than in
  !--- par because par%taumax/par%tauhomo are overwritten with the realized
  !--- values once the grid is built, so the factor cannot be recovered
  !--- from them afterward.
  real(wp), save :: dust_opac_norm = 1.0_wp

contains

  !=========================================================================
  ! Dust pseudo-number density from metallicity (Laursen+09).
  !   ndust = (Z / Z_ref) * (nH*xHI + f_ion * nH*(1-xHI))
  ! Multiplied by cext_dust * distance2cm in the caller to get the opacity.
  !=========================================================================
  elemental function laursen09_ndust(nH, xHI, Z, Z_ref, f_ion) result(ndust)
    real(wp), intent(in) :: nH, xHI, Z, Z_ref, f_ion
    real(wp) :: ndust, nHI, nHII
    nHI   = nH * xHI
    nHII  = nH * (1.0_wp - xHI)
    ndust = (Z / max(Z_ref, 1.0e-30_wp)) * (nHI + f_ion * nHII)
  end function laursen09_ndust

end module physics_amr_mod
