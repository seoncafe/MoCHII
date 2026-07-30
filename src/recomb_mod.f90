module recomb_mod
!---------------------------------------------------------------------------
! MoCHII: recombination and collisional-ionization rate coefficients.
!
! Recombination — selected by par%recomb_model, for H II, He II, He III as
! functions of T [K].
!   'badnell_milne' (default): total radiative recombination alpha_A from the
!     Badnell (2023) RR fit; the ground-level (n=1) direct-capture alpha_1
!     from the Milne relation applied to the transport photoionization cross
!     sections (see alpha1_* below); case B alpha_B = alpha_A - alpha_1.
!     He II -> He I includes the Badnell three-term DR sum (negligible below
!     ~5e4 K); H II and He III have no DR.  The name carries both sources
!     because the two coefficients come from two determinations, and that
!     alpha_B is a difference rather than a fit is what the name records.
!     Only the FIT FORM of alpha_1 is still Mao & Kaastra's, hence rr_mao.
!   'hui_gnedin': Hui & Gnedin (1997, MNRAS 292, 27) case-A/B fits, where one
!     source supplies both alpha_A and alpha_B directly.  Retained so the
!     recorded gates reproduce.
!
! Collisional ionization — Voronov (1997, ADNDT 65, 1) fits for H I, He I,
! He II:  k(T) = A (1 + P sqrt(U)) U^K exp(-U) / (X + U),  U = dE/T_eV.
! Negligible next to photoionization at ~1e4 K but kept in the balance.
!---------------------------------------------------------------------------
  use define, only : wp, kboltz_cgs, ev2erg, nan64, par
  implicit none
  private

  public :: alphaA_HII, alphaB_HII, alphaA_HeII, alphaB_HeII
  public :: alphaA_HeIII, alphaB_HeIII
  public :: ci_HI, ci_HeI, ci_HeII, ci_dere_ratio

  !--- ionization-threshold temperatures T_TR = dE/k [K] (Hui & Gnedin).
  real(kind=wp), parameter :: T_TR_HI   = 157807.0_wp
  real(kind=wp), parameter :: T_TR_HeI  = 285335.0_wp
  real(kind=wp), parameter :: T_TR_HeII = 631515.0_wp

contains

  !=========================================================================
  ! Public recombination coefficients: dispatch on par%recomb_model, which
  ! read_input has already restricted to the two names named here.  Each model
  ! is matched by name and there is no catch-all branch: a third model added to
  ! that whitelist but not to these select constructs returns NaN, where the
  ! former 'anything that is not hui_gnedin' would have served it the
  ! badnell_milne rate and reported success.
  !=========================================================================
  elemental real(kind=wp) function alphaA_HII(T) result(a)
    real(kind=wp), intent(in) :: T
    select case (trim(par%recomb_model))
    case ('badnell_milne')
       a = rr_badnell(T, 8.318e-11_wp, 0.7472_wp, 2.965_wp, 7.001e5_wp, 0.0_wp, 0.0_wp)
    case ('hui_gnedin')
       a = hg_alphaA_HII(T)
    case default
       a = nan64
    end select
  end function alphaA_HII

  elemental real(kind=wp) function alphaB_HII(T) result(a)
    real(kind=wp), intent(in) :: T
    select case (trim(par%recomb_model))
    case ('badnell_milne')
       a = alphaA_HII(T) - alpha1_HII(T)
    case ('hui_gnedin')
       a = hg_alphaB_HII(T)
    case default
       a = nan64
    end select
  end function alphaB_HII

  elemental real(kind=wp) function alphaA_HeII(T) result(a)
    real(kind=wp), intent(in) :: T
    select case (trim(par%recomb_model))
    case ('badnell_milne')
       a = rr_badnell(T, 5.235e-11_wp, 0.6988_wp, 7.301_wp, 4.475e6_wp, 0.0829_wp, 1.682e5_wp) &
           + dr_HeII_badnell(T)
    case ('hui_gnedin')
       a = hg_alphaA_HeII(T)
    case default
       a = nan64
    end select
  end function alphaA_HeII

  elemental real(kind=wp) function alphaB_HeII(T) result(a)
    real(kind=wp), intent(in) :: T
    select case (trim(par%recomb_model))
    case ('badnell_milne')
       a = alphaA_HeII(T) - alpha1_HeII(T)
    case ('hui_gnedin')
       a = hg_alphaB_HeII(T)
    case default
       a = nan64
    end select
  end function alphaB_HeII

  elemental real(kind=wp) function alphaA_HeIII(T) result(a)
    real(kind=wp), intent(in) :: T
    select case (trim(par%recomb_model))
    case ('badnell_milne')
       a = rr_badnell(T, 1.818e-10_wp, 0.7492_wp, 1.017e1_wp, 2.786e6_wp, 0.0_wp, 0.0_wp)
    case ('hui_gnedin')
       a = hg_alphaA_HeIII(T)
    case default
       a = nan64
    end select
  end function alphaA_HeIII

  elemental real(kind=wp) function alphaB_HeIII(T) result(a)
    real(kind=wp), intent(in) :: T
    select case (trim(par%recomb_model))
    case ('badnell_milne')
       a = alphaA_HeIII(T) - alpha1_HeIII(T)
    case ('hui_gnedin')
       a = hg_alphaB_HeIII(T)
    case default
       a = nan64
    end select
  end function alphaB_HeIII

  !=========================================================================
  ! Badnell (2023) radiative recombination fit and the He II -> He I
  ! dielectronic three-term sum.
  !=========================================================================
  elemental real(kind=wp) function rr_badnell(T, A, B, T0, T1, Cc, T2) result(rr)
    real(kind=wp), intent(in) :: T, A, B, T0, T1, Cc, T2
    real(kind=wp) :: bp, tt
    bp = B + Cc*exp(-T2/T)
    tt = sqrt(T/T0)
    rr = A / ( tt * (1.0_wp+tt)**(1.0_wp-bp) * (1.0_wp + sqrt(T/T1))**(1.0_wp+bp) )
  end function rr_badnell

  elemental real(kind=wp) function dr_HeII_badnell(T) result(dr)
    real(kind=wp), intent(in) :: T
    dr = T**(-1.5_wp) * ( 1.417e-3_wp*exp(-4.633e5_wp/T) &
                        + 2.235e-4_wp*exp(-5.532e5_wp/T) &
                        - 2.185e-5_wp*exp(-8.887e5_wp/T) )
  end function dr_HeII_badnell

  !=========================================================================
  ! Fit FORM of Mao & Kaastra (2016, A&A 587, A84) for a level radiative-
  ! recombination coefficient.  T [K] is converted to eV; a0 is in units of
  ! 1e-10 cm^3 s^-1.  The form is theirs; the coefficients the alpha1_* below
  ! pass in are NOT — they come from the Milne relation (see there).
  !=========================================================================
  elemental real(kind=wp) function rr_mao(T, a0, b0, c0, a1, b1, a2, b2) result(rr)
    real(kind=wp), intent(in) :: T, a0, b0, c0, a1, b1, a2, b2
    real(kind=wp) :: te
    te = T * (kboltz_cgs/ev2erg)
    rr = a0*1.0e-10_wp * te**(-b0 - c0*log(te)) &
         * (1.0_wp + a2*te**(-b2)) / (1.0_wp + a1*te**(-b1))
  end function rr_mao

  !=========================================================================
  ! n=1 ground-level direct-capture rate alpha_1(T), the term that turns the
  ! Badnell total alpha_A into case B.
  !
  ! The coefficients are fits to the Milne relation evaluated with the SAME
  ! ground-state photoionization cross sections the transport uses — exact
  ! hydrogenic for H I(1s) and He II(1s), VFKY96 for He I(1^1S), all from
  ! photo_xsec —
  !
  !   alpha_1(T) = (g_n/2 g_+) (2/sqrt(pi)) (kT)^(-3/2) sqrt(2/m_e)/(m_e c^2)
  !                Int E^2 sigma_pi(E) exp[-(E-E_th)/kT] dE ,
  !
  ! the same integral milne_recomb_spectrum_mod samples the diffuse ground
  ! continua from.  Generated by tools/fitting/fit_alpha1_milne.py; max
  ! |fit/Milne - 1| over 1e3-1e5 K is 0.019% (H I), 0.001% (He I), 0.007%
  ! (He II).
  !
  ! WHY NOT the Mao & Kaastra (2016) level-resolved n=1 fit these coefficients
  ! replace: alpha_B = alpha_A - alpha_1 with alpha_A from Badnell and alpha_1
  ! from Mao subtracts two unrelated determinations, so alpha_B inherits their
  ! disagreement.  For He I that disagreement is 2.4% of alpha_1.  Detailed
  ! balance fixes alpha_1 exactly from sigma_pi, so taking it from the Milne
  ! relation makes the opacity, the diffuse ground-continuum spectrum and the
  ! ground recombination rate one consistent set instead of three inputs.
  !
  ! Cloudy c25.00 derives its H-like/He-like recombination the same way, level
  ! by level from its own cross sections (source/iso_radiative_recomb.cpp;
  ! ground level = row 0 of data/h_iso_recomb.dat, data/he_iso_recomb.dat).
  ! alpha_1 at 1e4 K [cm^3 s^-1], these fits vs Cloudy vs the Mao coefficients:
  !   H I    1.5791e-13   1.5840e-13   1.5873e-13
  !   He I   1.6086e-13   1.6032e-13   1.5650e-13
  !   He II  6.5117e-13   6.5145e-13   6.5366e-13
  ! i.e. <= 0.4% from Cloudy on all three, where Mao was 2.4% low on He I.
  ! Gate: tests/milne/check_alpha1.f90.
  !=========================================================================
  elemental real(kind=wp) function alpha1_HII(T) result(a)
    real(kind=wp), intent(in) :: T
    a = rr_mao(T, 1.87007e-2_wp, 1.29077_wp, -2.17465e-3_wp, 1.29189e1_wp, &
               8.05060e-1_wp, 8.44399e-2_wp, 5.93940e-4_wp)
  end function alpha1_HII

  elemental real(kind=wp) function alpha1_HeIII(T) result(a)
    real(kind=wp), intent(in) :: T
    a = rr_mao(T, 2.61525e-1_wp, 1.39527_wp, -3.22079e-4_wp, 6.40708e1_wp, &
               8.97582e-1_wp, 5.01551e-1_wp, 1.21533e-3_wp)
  end function alpha1_HeIII

  elemental real(kind=wp) function alpha1_HeII(T) result(a)
    real(kind=wp), intent(in) :: T
    a = rr_mao(T, 3.26862e-3_wp, 6.97361e-1_wp, 3.14907e-4_wp, 5.76341e1_wp, &
               1.18534_wp, 2.58417e1_wp, 9.90134e-1_wp)
  end function alpha1_HeII

  !=========================================================================
  ! Hui & Gnedin (1997) recombination fits; lambda = 2 T_TR / T.
  !=========================================================================
  elemental real(kind=wp) function hg_alphaA_HII(T) result(a)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*T_TR_HI/T
    a = 1.269e-13_wp * lam**1.503_wp / (1.0_wp + (lam/0.522_wp)**0.470_wp)**1.923_wp
  end function hg_alphaA_HII

  elemental real(kind=wp) function hg_alphaB_HII(T) result(a)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*T_TR_HI/T
    a = 2.753e-14_wp * lam**1.500_wp / (1.0_wp + (lam/2.740_wp)**0.407_wp)**2.242_wp
  end function hg_alphaB_HII

  elemental real(kind=wp) function hg_alphaA_HeII(T) result(a)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*T_TR_HeI/T
    a = 3.000e-14_wp * lam**0.654_wp
  end function hg_alphaA_HeII

  elemental real(kind=wp) function hg_alphaB_HeII(T) result(a)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*T_TR_HeI/T
    a = 1.260e-14_wp * lam**0.750_wp
  end function hg_alphaB_HeII

  !--- He III is hydrogenic (Z=2): the H fit evaluated on lambda(T_TR_HeII),
  !--- scaled by Z (alpha_Z(T) = Z alpha_H(T/Z^2)).
  elemental real(kind=wp) function hg_alphaA_HeIII(T) result(a)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*T_TR_HeII/T
    a = 2.0_wp * 1.269e-13_wp * lam**1.503_wp / (1.0_wp + (lam/0.522_wp)**0.470_wp)**1.923_wp
  end function hg_alphaA_HeIII

  elemental real(kind=wp) function hg_alphaB_HeIII(T) result(a)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: lam
    lam = 2.0_wp*T_TR_HeII/T
    a = 2.0_wp * 2.753e-14_wp * lam**1.500_wp / (1.0_wp + (lam/2.740_wp)**0.407_wp)**2.242_wp
  end function hg_alphaB_HeIII

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
    k = voronov(T, 13.6_wp, 0.0_wp, 2.91e-8_wp, 0.232_wp, 0.39_wp) &
        * ci_dere_ratio(1, 1)
  end function ci_HI

  elemental real(kind=wp) function ci_HeI(T) result(k)
    real(kind=wp), intent(in) :: T
    k = voronov(T, 24.6_wp, 0.0_wp, 1.75e-8_wp, 0.180_wp, 0.35_wp) &
        * ci_dere_ratio(2, 1)
  end function ci_HeI

  elemental real(kind=wp) function ci_HeII(T) result(k)
    real(kind=wp), intent(in) :: T
    k = voronov(T, 54.4_wp, 1.0_wp, 2.05e-9_wp, 0.265_wp, 0.25_wp) &
        * ci_dere_ratio(2, 2)
  end function ci_HeII

  !=========================================================================
  ! Dere (2007) / Voronov ('Dima') collisional-ionization ratio.  Cloudy
  ! c23.01's HYBRID model multiplies the Voronov rate by this constant per
  ! (element, stage) factor (source/atmdat_adfa.cpp DereRatio, "evaluated
  ! where the ion is abundant").  stage = MoCHII transition index
  ! (1 = X0 -> X+).  Returns 1 unless par%ci_model = 'dere_hybrid'.
  !=========================================================================
  pure real(kind=wp) function ci_dere_ratio(Z, stage) result(r)
    integer, intent(in) :: Z, stage
    real(kind=wp) :: tab(4)
    integer :: n
    r = 1.0_wp
    if (trim(par%ci_model) /= 'dere_hybrid') return
    tab = 1.0_wp;  n = 0
    select case (Z)
    case (1);  tab(1)   = 0.9063_wp;                                    n = 1  ! H
    case (2);  tab(1:2) = [1.0389_wp, 1.0686_wp];                       n = 2  ! He
    case (6);  tab(1:3) = [1.0499_wp, 0.913_wp, 1.0377_wp];             n = 3  ! C
    case (7);  tab(1:2) = [1.0421_wp, 1.1966_wp];                       n = 2  ! N
    case (8);  tab(1:2) = [1.041_wp,  1.1181_wp];                       n = 2  ! O
    case (10); tab(1:2) = [0.8089_wp, 1.1395_wp];                       n = 2  ! Ne
    case (12); tab(1:2) = [0.3793_wp, 0.9857_wp];                       n = 2  ! Mg
    case (14); tab(1:4) = [0.7328_wp, 0.8798_wp, 0.4492_wp, 0.8221_wp]; n = 4  ! Si
    case (16); tab(1:3) = [1.3572_wp, 0.8925_wp, 0.8119_wp];            n = 3  ! S
    case (17); tab(1:4) = [0.5412_wp, 0.8428_wp, 0.9237_wp, 0.819_wp];  n = 4  ! Cl
    case (18); tab(1:3) = [0.9242_wp, 0.8644_wp, 0.9752_wp];            n = 3  ! Ar
    case (20); tab(1:4) = [0.7652_wp, 1.1668_wp, 1.0422_wp, 0.8705_wp]; n = 4  ! Ca
    case (26); tab(1:3) = [0.9904_wp, 1.0568_wp, 1.824_wp];             n = 3  ! Fe
    end select
    if (stage >= 1 .and. stage <= n) r = tab(stage)
  end function ci_dere_ratio

end module recomb_mod
