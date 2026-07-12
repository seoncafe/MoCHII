! MoCHII: copied from MoCafe_v2.00/src/physics_amr_mod.f90 (2026-07-11)
module physics_amr_mod
!---------------------------------------------------------------------------
! Dust physics for the AMR grid (dust-only port of LaRT physics_amr_mod.f90).
!
! Only the metallicity->dust relation (Laursen+09) is kept; the CIE neutral
! fraction (T->xHI) and the Case-B Lyman-alpha emissivity are deferred with
! temperature (AMR_CLUMPS_PLAN.md Part A / Section 6).  For dust_model =
! 'laursen09' MoCafe therefore requires an explicit xHI column in the file.
!---------------------------------------------------------------------------
  use define, only: wp
  implicit none
  private

  public :: laursen09_ndust

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
