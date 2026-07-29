module hei_cascade_mod
!---------------------------------------------------------------------------
! MoCHII: branch weights of the He I recombination cascade to the ground
! state, for the diffuse field's fourth channel.
!
! Every case-B He II recombination reaches 1^1S through exactly one of the
! n = 2 terms, so what the channel needs are the three BRANCH FRACTIONS; the
! absolute rate is already set by alpha_B(He II) in diffuse_build.  The
! fractions come from one internally consistent set of effective
! recombination coefficients,
!
!     alpha(2^3S) = 2.10e-13 (T/1e4)^-0.381
!     alpha(2^1S) = 2.06e-14 (T/1e4)^-0.451
!     alpha(2^1P) = 4.17e-14 (T/1e4)^-0.695     [cm^3 s^-1]
!
! from Benjamin, Skillman & Smits (1999) Table 1 with J. S. Mathis, as used
! in K. Wood's ionize (probset.f).  Taking all three from one source matters
! more than the individual values: a fraction assembled from coefficients
! with unrelated temperature dependences is not a fraction of anything.
!
! Normalizing by their SUM rather than against an external alpha_B is Wood's
! own choice, made for the reason its source states -- "different
! T-dependences in fitting formulae".  The three fits account for 97% of
! MoCHII's alpha_B(He II) at 10^4 K and 90% at 8000 K, and forcing them onto
! alpha_B would dump that mismatch into whichever branch was solved for last.
! The shortfall is reported at setup instead.
!
! At 8000 K the fractions are 0.762 : 0.076 : 0.162, against the low-density
! statistical weights 0.750 : 0.083 : 0.167 the code used before.  The
! correction is modest and temperature dependent, not a reordering.
!
! An earlier version of this module read nebcont_mod's 6.23e-14 as
! alpha_eff(2^1S) and closed 2^1P against alpha_B, which gave 0.746 : 0.224 :
! 0.031 and an accompanying argument that case-B trapping inverts the singlet
! split.  That was wrong: 2.06e-14 + 4.17e-14 = 6.23e-14 exactly, so the
! coefficient is the singlet TOTAL, not the 2^1S branch.
!
! NOTE ON THE 2^3S COEFFICIENT.  hei_metastable_mod carries a different fit,
! 2.10e-13 (T/1e4)^-0.778 (Oklopcic & Hirata 2018): the same prefactor with a
! steeper exponent, 9% higher at 8000 K.  That one is an ABSOLUTE production
! rate for the metastable population; this one is a member of a set whose
! ratios are the branch fractions.  They serve different purposes and are
! deliberately not shared, because mixing sources would break the internal
! consistency the fractions depend on.
!---------------------------------------------------------------------------
  use define
  implicit none
  private

  public :: hei_cascade_setup, hei_branch_fractions, hei_584_ionizes_H

  !--- BSS99 Table 1 + JSM, as tabulated in Wood's ionize/probset.f:
  !--- prefactor [cm^3 s^-1] and exponent of (T/1e4).
  real(kind=wp), parameter :: A_2t3S(2) = [2.10e-13_wp, -0.381_wp]
  real(kind=wp), parameter :: A_2s1S(2) = [2.06e-14_wp, -0.451_wp]
  real(kind=wp), parameter :: A_2p1P(2) = [4.17e-14_wp, -0.695_wp]
  !--- the fits are nebular; clamp rather than extrapolate
  real(kind=wp), parameter :: T_LO = 3.0e3_wp, T_HI = 2.0e4_wp

contains

  real(kind=wp) function alpha_fit(T, cp) result(a)
    real(kind=wp), intent(in) :: T, cp(2)
    a = cp(1)*(min(max(T, T_LO), T_HI)/1.0e4_wp)**cp(2)
  end function alpha_fit

  !=========================================================================
  ! Branch fractions of the case-B He I cascade at T.  They sum to one by
  ! construction, so the channel emits exactly the recombinations that
  ! alpha_B(He II) feeds it.
  !=========================================================================
  subroutine hei_branch_fractions(T, f3s, f1s, f1p)
    real(kind=wp), intent(in)  :: T
    real(kind=wp), intent(out) :: f3s, f1s, f1p
    real(kind=wp) :: a3, a1s, a1p, s
    a3  = alpha_fit(T, A_2t3S)
    a1s = alpha_fit(T, A_2s1S)
    a1p = alpha_fit(T, A_2p1P)
    s   = a3 + a1s + a1p
    f3s = a3/s;  f1s = a1s/s;  f1p = a1p/s
  end subroutine hei_branch_fractions

  !=========================================================================
  ! Fraction of 2^1P decays whose 584 A photon ends up photoionizing a
  ! hydrogen atom, the rest being collisionally transferred 2^1P -> 2^1S and
  ! re-emitted as the two-photon continuum.
  !
  ! He I 584 A is a resonance line: the medium is optically thick to it in
  ! He^0, so the photon scatters many times before anything irreversible
  ! happens to it.  Two irreversible outcomes compete during that walk --
  ! absorption by H^0 (21.22 > 13.60 eV, so it photoionizes and deposits
  ! 7.62 eV of heat) and collisional transfer to 2^1S in the He atom that
  ! last absorbed it.  Which wins is set by the H^0 and He^0 opacities the
  ! trapped photon sees, hence this escape-probability closure
  !
  !     p = [1 + 77 x(He^0)/(sqrt(T) x(H^0))]^-1 ,
  !
  ! from Wood, Mathis & Ercolano (2004) sect. 3.3.  The same expression is in
  ! K. Wood's ionize (mcionize.f, h0find.f, written as
  ! 1/(1 + 0.77/(sqrt(T)/100) x(He^0)/x(H^0))) and in CMacIonize
  ! (PhysicalDiffuseReemissionHandler, written as
  ! sqrt(T) x(H^0)/(sqrt(T) x(H^0) + 77 x(He^0))); the three forms are
  ! algebraically identical.
  !
  ! DOMAIN OF VALIDITY.  The closure assumes the 584 A line is optically
  ! thick, static and dust free.  Dust competes for the trapped photon and is
  ! not in it: at the Lexington benchmark dust content that is a small
  ! correction, but at higher dust-to-gas it is not, and the closure will
  ! then overestimate both outcomes it does describe.  Bulk motions that
  ! shift the line out of resonance are likewise absent.
  !
  ! Limits: x(H^0) -> 0 gives p -> 0, nothing left to ionize, everything
  ! converts; x(He^0) -> 0 gives p -> 1, nothing left to convert.
  !=========================================================================
  real(kind=wp) function hei_584_ionizes_H(T, xHI, xHeI) result(p)
    real(kind=wp), intent(in) :: T, xHI, xHeI
    p = 1.0_wp/(1.0_wp + 77.0_wp/sqrt(max(T, tinest)) &
                *max(xHeI, 0.0_wp)/max(xHI, tinest))
  end function hei_584_ionizes_H

  !=========================================================================
  ! Report the fractions and how much of alpha_B(He II) the three fits
  ! account for.  The shortfall is a statement about two independent atomic
  ! data sets and is printed, not absorbed.
  !=========================================================================
  subroutine hei_cascade_setup()
    use recomb_mod, only : alphaB_HeII
    implicit none
    real(kind=wp), parameter :: TCHK(4) = &
       [5.0e3_wp, 8.0e3_wp, 1.0e4_wp, 1.5e4_wp]
    real(kind=wp) :: T, f3s, f1s, f1p, s
    integer :: i

    if (mpar%p_rank /= 0) return
    write(*,'(a)') ' HECAS: He I cascade branch fractions (BSS99 Table 1 ' // &
       '+ JSM, one consistent set)'
    write(*,'(a)') '        T [K]      2^3S      2^1S      2^1P   ' // &
       'sum(alpha)/alphaB'
    do i = 1, 4
       T = TCHK(i)
       call hei_branch_fractions(T, f3s, f1s, f1p)
       s = alpha_fit(T, A_2t3S) + alpha_fit(T, A_2s1S) + alpha_fit(T, A_2p1P)
       write(*,'(a,f9.0,3f10.5,f14.5)') '       ', T, f3s, f1s, f1p, &
          s/alphaB_HeII(T)
    end do
    write(*,'(a)') '        (the low-density statistical weights would be' // &
       '   0.75000   0.08333   0.16667)'
  end subroutine hei_cascade_setup

end module hei_cascade_mod
