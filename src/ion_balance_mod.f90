module ion_balance_mod
!---------------------------------------------------------------------------
! MoCHII: H/He photoionization equilibrium at fixed T_e.
!
! Per leaf, given the rate integrals Gamma_i from the ionizing-band J tally
! and T_e = par%te_fixed, solve
!     n_HI  (Gamma_HI  + n_e C_HI )  = alpha_[AB](H II)   n_e n_HII
!     n_HeI (Gamma_HeI + n_e C_HeI)  = alpha_[AB](He II)  n_e n_HeII
!     n_HeII(Gamma_HeII+ n_e C_HeII) = alpha_[AB](He III) n_e n_HeIII
! with the closure  n_e = nH [ x_HII + y_He (x_HeII + 2 x_HeIII) ]
! by damped fixed-point iteration on n_e.
! Case A/B via par%case_ab ('B' = on-the-spot, the default).
! With par%helium_charge_exchange (default on) the He-H charge-exchange
! reactions join both balances (cx_HeII_H and its companions below).
!
! The converged fractions are written to the shared gas-state arrays with
! under-relaxation par%ion_relax (1 = none; < 1 damps oscillation at sharp
! I-fronts, a known failure mode).
! Only the node-local rank 0 runs the solve (the inputs are ALLREDUCEd
! tallies, so every node's rank 0 is identical); it writes the shared
! arrays directly and broadcasts max|delta x_HII| to the other node ranks,
! which skip that work entirely.
!---------------------------------------------------------------------------
  use define
  use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, gas_ne, &
                            gas_Te, gas_nleaf
  use gas_rates_mod, only : gamma_HI, gamma_HeI, gamma_HeII, &
                            secondary_ionization
  use photoelectron_degradation_mod, only : ionized_fraction_svs
  use recomb_mod
  implicit none
  private

  public :: gas_equilibrium_update, solve_ion_cell, charge_neutrality_ne
  public :: helium_change_measures
  public :: hydrogen_balance_residual, hydrogen_balance_residual_grid
  public :: cx_HeII_H, cx_HeI_Hp, cx_HeIII_H, cx_HeII_Hp

  !--- ionization-potential differences for the detailed-balance reverses
  !--- [eV] (NIST ASD: H 13.598434, He 24.587387, He+ 54.417765):
  !--- He + H+ -> He+ + H and He+ + H+ -> He++ + H
  real(kind=wp), parameter :: DE_HEI_H  = 24.587387_wp - 13.598434_wp
  real(kind=wp), parameter :: DE_HEII_H = 54.417765_wp - 13.598434_wp

contains

  !=========================================================================
  ! Single-cell H/He equilibrium at temperature T: damped fixed point on
  ! n_e.  x* are in/out (start = current state).
  ! Shared by gas_equilibrium_update (fixed Te) and thermal_mod (Te solve).
  !=========================================================================
  ! With par%metal_ne and the optional leaf index il, electrons from the
  ! metal cascade join the n_e closure (species_ne; beyond the I-front
  ! they are the only electrons, e.g. n_e ~ n_CII in the PDR zone).
  !
  ! With par%charge_exchange = 'two_way' the reactions
  !     X^i + H^+ -> X^(i+1) + H^0        (rate R_cx_in  per H^+)
  !     X^(i+1) + H^0 -> X^i + H^+        (rate R_cx_out per H^0)
  ! enter hydrogen's balance as well as the metal cascade's, so the charge
  ! they move is conserved in the reacting pair:
  !     n_HI (Gamma_H + C_H n_e + R_cx_out) = (alpha_H n_e + R_cx_in) n_HII.
  ! Eliminating n_HII = n_H (1 - x_HI) keeps the closed form; BOTH terms
  ! appear in the denominator, R_cx_in for the same reason alpha_H n_e does.
  !
  ! With par%helium_charge_exchange the He-H reactions (rate coefficients
  ! k21 = cx_HeII_H, k12 = cx_HeI_Hp, k32 = cx_HeIII_H, k23 = cx_HeII_Hp)
  ! add, for each H,
  !     R_he_out = n_He (x_HeII k21 + x_HeIII k32)    (H0 -> H+)
  !     R_he_in  = n_He (x_HeI  k12 + x_HeII  k23)    (H+ -> H0)
  ! to hydrogen's balance, and n_HI k21 / n_HII k12 (k32 / k23 for the next
  ! stage) to helium's recombination / ionization, so the charge each
  ! reaction moves is counted on both sides.  The He fractions in R_he are
  ! those of the previous iterate (the lag nem carries); because an exchange
  ! moves charge without changing n_e, the loop then also requires x_HI and
  ! x_HeII to have stopped moving before it exits.
  subroutine solve_ion_cell(gH, gHe1, gHe2, nH, T, caseA, &
                            xHI, xHeI, xHeII, ne, il)
    use species_mod, only : species_ne_prepare, species_metal_sums, n_elements
    real(kind=wp), intent(in)    :: gH, gHe1, gHe2, nH, T
    logical,       intent(in)    :: caseA
    real(kind=wp), intent(inout) :: xHI, xHeI, xHeII
    real(kind=wp), intent(out)   :: ne
    integer,       intent(in), optional :: il
    real(kind=wp) :: aH, aHe2, aHe3, cH, cHe1, cHe2, yHe
    real(kind=wp) :: ne_old, xHeIII, r1, r2, denH
    real(kind=wp) :: nem, r_cx_in, r_cx_out
    real(kind=wp) :: k21, k12, k32, k23, rhe_in, rhe_out, xHI_old, xHeII_old
    real(kind=wp) :: gHt, gHe1t, dgHI, dgHeI, dh1, dh2, dh3
    logical :: with_metal_ne, with_metal_cx, with_he_cx, with_sec, settled
    integer :: it

    !--- the two feedbacks are independent: counting the metal electrons is
    !--- par%metal_ne, conserving charge in the charge-exchange pair is
    !--- par%charge_exchange.  Either one needs the (leaf, T) rate cache.
    with_metal_ne = .false.
    with_metal_cx = .false.
    if (present(il) .and. par%use_metals .and. n_elements > 0) then
       with_metal_ne = par%metal_ne
       with_metal_cx = trim(par%charge_exchange) == 'two_way'
    end if
    !--- T is fixed inside the fixed point: evaluate the metal rate
    !--- coefficients once (species_metal_sums reduces each iteration
    !--- below to multiply-adds).
    if (with_metal_ne .or. with_metal_cx) call species_ne_prepare(il, T)

    yHe = par%abund_He
    if (caseA) then
       aH = alphaA_HII(T);  aHe2 = alphaA_HeII(T);  aHe3 = alphaA_HeIII(T)
    else
       aH = alphaB_HII(T);  aHe2 = alphaB_HeII(T);  aHe3 = alphaB_HeIII(T)
    end if
    cH = ci_HI(T);  cHe1 = ci_HeI(T);  cHe2 = ci_HeII(T)
    with_he_cx = par%helium_charge_exchange .and. yHe > 0.0_wp
    if (with_he_cx) then
       k21 = cx_HeII_H(T);   k12 = cx_HeI_Hp(T)
       k32 = cx_HeIII_H(T);  k23 = cx_HeII_Hp(T)
    end if
    rhe_in = 0.0_wp;  rhe_out = 0.0_wp
    !--- secondary ionization by fast photoelectrons (par%use_sec_ion) at the
    !--- current iterate: x, n_HI and n_HeI of this very solve, not of the
    !--- previous transport pass.
    with_sec = par%use_sec_ion .and. present(il)
    gHt = gH;  gHe1t = gHe1

    xHeIII = max(0.0_wp, 1.0_wp - xHeI - xHeII)
    ne = nH * ((1.0_wp - xHI) + yHe*(xHeII + 2.0_wp*xHeIII))
    ne = max(ne, 1.0e-12_wp*nH)
    nem = 0.0_wp;  r_cx_in = 0.0_wp;  r_cx_out = 0.0_wp
    do it = 1, 200
       ne_old = ne
       xHI_old = xHI;  xHeII_old = xHeII
       if (with_sec) then
          call secondary_ionization(il, ionized_fraction_svs(ne, nH, yHe*nH), &
               nH*xHI, nH*yHe*xHeI, nH*yHe*xHeII, dgHI, dgHeI, dh1, dh2, dh3)
          gHt = gH + dgHI;  gHe1t = gHe1 + dgHeI
       end if
       if (with_he_cx) then
          rhe_out = yHe*nH*(xHeII*k21 + xHeIII*k32)
          rhe_in  = yHe*nH*(xHeI*k12  + xHeII*k23)
       end if
       !--- the stage densities the two charge-exchange sums need depend on
       !--- n_HI and n_HII through the charge-exchange terms of the cascade
       !--- itself, so they are taken from the previous iterate of this fixed
       !--- point — the same lag the electron sum carries.
       !--- This lag is what nem carries.  It used to bind x_HI as well, and
       !--- there it was the whole difficulty: substituting into hydrogen's
       !--- balance with the charge-exchange terms lagged gives a map of slope
       !--- R_in/(R_in + alpha_H n_e), which tends to 1 in partially neutral
       !--- gas, so the exit test below was met only after 124 passes at
       !--- x_HII = 0.18, 239 at 0.10 and 1285 at 0.021 against the 200-pass
       !--- cap, and such a cell left this loop with x_HII high by ~8%.
       !--- hydrogen_neutral_fraction removes that by solving the balance at
       !--- this n_e rather than substituting into it: one call now reaches the
       !--- converged state to 3.6e-11 at every x_HII down to 0.02.  What
       !--- remains lagged is only nem, whose own map is a long way from
       !--- marginal (md/ION_BALANCE_CONVERGENCE_PLAN.md, S1 to S3).
       if (with_metal_ne .or. with_metal_cx) then
          call species_metal_sums(nH, ne_old, nH*xHI, nH*(1.0_wp - xHI), &
                                  nem, r_cx_in, r_cx_out)
          if (.not. with_metal_cx) then
             r_cx_in  = 0.0_wp
             r_cx_out = 0.0_wp
          end if
       end if
       if (with_metal_cx .and. with_he_cx) then
          !--- the He terms are constants of this root find (lagged)
          xHI = hydrogen_neutral_fraction(nH, ne, aH*ne + rhe_in, &
                                          gHt + (cH + aH)*ne + rhe_in + rhe_out)
       else if (with_metal_cx) then
          !--- solve the balance at this n_e instead of substituting into it:
          !--- the charge-exchange terms make the substitution map marginally
          !--- stable in partially neutral gas (see hydrogen_neutral_fraction).
          xHI = hydrogen_neutral_fraction(nH, ne, aH*ne, gHt + (cH + aH)*ne)
       else if (with_he_cx) then
          denH = gHt + (cH + aH)*ne + r_cx_in + r_cx_out + rhe_in + rhe_out
          if (denH > 0.0_wp) then
             xHI = (aH*ne + r_cx_in + rhe_in) / denH
          else
             xHI = 1.0_wp
          end if
       else
          !--- kept verbatim, including the two terms that are zero on this
          !--- branch, so the generated arithmetic is unchanged and
          !--- 'metal_only' stays bit-identical to the pre-M1 code: under
          !--- -fp-model fast=2 even dropping an exact +0.0 lets the compiler
          !--- reassociate and move x_HeII by one ulp.
          denH = gHt + (cH + aH)*ne + r_cx_in + r_cx_out
          if (denH > 0.0_wp) then
             xHI = (aH*ne + r_cx_in) / denH
          else
             xHI = 1.0_wp
          end if
       end if
       !--- r1 = He+/He0, r2 = He++/He+.  Cap at 1e150 (physical values stay
       !--- below ~1e15 given the ne floor) so r1*r2 cannot overflow to Inf
       !--- and produce a 0*Inf = NaN in xHeII below.
       if (with_he_cx) then
          r1 = min((gHe1t + cHe1*ne + k12*nH*(1.0_wp - xHI)) &
                   / (aHe2*ne + k21*nH*xHI), 1.0e150_wp)
          r2 = min((gHe2 + cHe2*ne + k23*nH*(1.0_wp - xHI)) &
                   / (aHe3*ne + k32*nH*xHI), 1.0e150_wp)
       else
          r1 = min((gHe1t + cHe1*ne) / (aHe2*ne), 1.0e150_wp)
          r2 = min((gHe2 + cHe2*ne) / (aHe3*ne), 1.0e150_wp)
       end if
       xHeI   = 1.0_wp / (1.0_wp + r1 + r1*r2)
       xHeII  = xHeI * r1
       xHeIII = xHeII * r2
       ne = nH * ((1.0_wp - xHI) + yHe*(xHeII + 2.0_wp*xHeIII))
       if (with_metal_ne) ne = ne + nem
       ne = max(0.5_wp*(ne + ne_old), 1.0e-12_wp*nH)
       settled = .true.
       if (with_he_cx .or. with_sec) settled = &
          abs(xHI - xHI_old) <= 1.0e-10_wp*min(xHI, 1.0_wp - xHI) + 1.0e-15_wp &
          .and. abs(xHeII - xHeII_old) <= 1.0e-10_wp*xHeII + 1.0e-15_wp
       if (abs(ne - ne_old) <= 1.0e-10_wp*ne .and. settled) exit
    end do
  end subroutine solve_ion_cell

  !=========================================================================
  ! He-H charge-exchange rate coefficients [cm^3 s^-1]
  ! (par%helium_charge_exchange).  The Kingdon & Ferland (1996, ApJS 106,
  ! 205; KF96) Table 1 fits are their eq. (7),
  !   k = a 1e-9 t4^b [1 + c exp(d t4)],  t4 = T/1e4.
  !
  ! cx_HeII_H: He+ + H -> He + H+, the sum of two channels that the common
  ! calculation (Zygelman, Dalgarno, Kimura & Lane 1989, PRA 40, 2340) keeps
  ! apart:
  !  - radiative, He+ + H -> He + H+ + photon, through the A 1Sigma+ ->
  !    X 1Sigma+ transition: 1.25e-15 (T/300)^0.25 (Glover & Jappsen 2007,
  !    ApJ 666, 1, Table 1 R26, citing Zygelman et al. 1989).  That fit is
  !    one quarter of Zygelman's "Direct" column (5.36, 4.21, 4.41, 4.50,
  !    4.83, 5.99e-15 at 1, 10, 100, 200, 400, 1000 K; the fit is 0.86x the
  !    quarter at 100 K and 1.13x at 1000 K), the quarter being the
  !    probability of approach in the singlet,
  !    (2S+1)(2-delta_0Lambda) / [(2L_a+1)(2S_a+1)(2L_b+1)(2S_b+1)] = 1/4
  !    for He+(2S) + H(2S) (Courtney et al. 2021, ApJ 919, 70, eqs. 4-6).
  !    Zygelman tabulates only to 1000 K, so above that the T^0.25 is an
  !    extrapolation; Courtney et al. recompute the process and find the
  !    smooth background in agreement with Zygelman.
  !  - non-radiative: KF96 He+ row, a = 7.47e-6, b = 2.06, c = 9.93,
  !    d = -3.89, fitted over 6e3-1e5 K to Zygelman's calculation.  Below
  !    6000 K the fit is evaluated as written: it falls as t4^2.06, as a
  !    curve-crossing channel must at low collision energy, whereas holding
  !    the 6000 K value would overstate it (10x at 1000 K).  Above 1e5 K the
  !    fit grows without bound and T is held at 1e5 K.
  !   At 8000 K the two give 2.84e-15 and 6.8e-15 cm^3 s^-1.
  ! cx_HeI_Hp: He + H+ -> He+ + H, the detailed-balance reverse of the
  ! non-radiative channel alone (the radiative channel has no collisional
  ! reverse):
  !   k = [g(He+) g(H) / (g(He) g(H+))] exp(-dE/kT) k_nr = 4 exp(-dE/kT) k_nr,
  ! dE = 10.989 eV; the reduced-mass factor (mu/mu')^(3/2) of the two pairs
  ! is 1 + 6e-4 and is dropped.  One process, so one calculation and its
  ! reverse.
  !   The direct calculation of this endothermic direction by Kimura, Lane,
  ! Dalgarno & Dixson (1993, ApJ 405, 801; Glover & Jappsen 2007 R27 fit it)
  ! is NOT used, for a reason read from the published paper (p. 802): its
  ! Table 1 gives charge-transfer cross sections from 20 eV only, above a
  ! 10.99 eV threshold, and the Table 3 rates it integrates to imply an
  ! apparent activation energy of 4.33 eV between 8000 and 1e4 K (15.2 eV
  ! between 6000 and 8000 K).  For an endothermic reaction the Tolman
  ! activation energy cannot fall below about E_th - (3/2) kT = 9.7 eV
  ! there, so the table is not reliable at nebular temperatures.  It was
  ! adopted for a few hours on 2026-09-24 (35x this rate at 1e4 K, 0.16x at
  ! 3e4 K); the Lexington benchmarks moved by 0.02 K or less either way.
  ! cx_HeIII_H: He++ + H -> He+ + H+, KF96 He++ row, a = 1e-5, b = c = d = 0
  ! (1e3-1e7 K; West, Lane & Cohen 1982).  cx_HeII_Hp is its detailed-
  ! balance reverse: the statistical-weight ratio g(He++) g(H) /
  ! (g(He+) g(H+)) is 1 and dE = 40.819 eV.
  !=========================================================================
  real(kind=wp) function cx_HeII_H(T) result(k)
    real(kind=wp), intent(in) :: T
    k = 1.25e-15_wp*(T/300.0_wp)**0.25_wp + nonradiative_HeII_H(T)
  end function cx_HeII_H

  real(kind=wp) function cx_HeI_Hp(T) result(k)
    real(kind=wp), intent(in) :: T
    k = 4.0_wp*exp(-DE_HEI_H*ev2erg/(kboltz_cgs*T))*nonradiative_HeII_H(T)
  end function cx_HeI_Hp

  real(kind=wp) function cx_HeIII_H(T) result(k)
    real(kind=wp), intent(in) :: T
    k = 1.0e-14_wp + 0.0_wp*T
  end function cx_HeIII_H

  real(kind=wp) function cx_HeII_Hp(T) result(k)
    real(kind=wp), intent(in) :: T
    k = exp(-DE_HEII_H*ev2erg/(kboltz_cgs*T))*cx_HeIII_H(T)
  end function cx_HeII_Hp

  real(kind=wp) function nonradiative_HeII_H(T) result(k)
    real(kind=wp), intent(in) :: T
    real(kind=wp) :: t4
    t4 = min(T, 1.0e5_wp)/1.0e4_wp
    k = 7.47e-15_wp*t4**2.06_wp*(1.0_wp + 9.93_wp*exp(-3.89_wp*t4))
  end function nonradiative_HeII_H

  !=========================================================================
  ! Neutral hydrogen fraction that balances ionization against recombination
  ! at fixed n_e, WITH the charge-exchange terms, whose rates depend on x_HI
  ! themselves through the metal cascade.  Writing that balance as a residual,
  !
  !     R(x) = x [B + R_in(x) + R_out(x)] - [A + R_in(x)]
  !     A = alpha_H n_e ,  B = Gamma_H + (C_H + alpha_H) n_e
  !
  ! the root is bracketed in [0,1] BY CONSTRUCTION, with no clamp needed:
  ! R(0) = -(A + R_in(0)) < 0 because A > 0, and R(1) = B + R_out(1) > 0
  ! because B > 0.  Illinois (modified false position) therefore keeps the
  ! bracket at every step while converging superlinearly, so this cannot
  ! return an unphysical fraction — the property plain substitution has and a
  ! bare secant would forfeit.
  !
  ! Why a root find rather than substituting.  Where the near-resonant
  ! O0 + H+ pair locks oxygen to hydrogen the two charge-exchange terms cancel
  ! identically in this balance, so the substitution map has slope
  ! R_in/(R_in + alpha_H n_e), which tends to 1 in partially neutral gas: it
  ! converges, but takes 1285 passes at x_HII = 0.021 against a 200-pass cap.
  ! Every trial x here re-evaluates the whole cascade at that x, so unlike an
  ! extrapolator it never leaves the state off the slow manifold — which is why
  ! Aitken on n_e made the pass count worse and this does not
  ! (md/ION_BALANCE_CONVERGENCE_PLAN.md, S1 to S3).
  !=========================================================================
  real(kind=wp) function hydrogen_neutral_fraction(nH, ne, A, B) result(x)
    use species_mod, only : species_metal_sums
    real(kind=wp), intent(in) :: nH, ne, A, B
    real(kind=wp) :: xa, xb, Ra, Rb, Rx, nem, rin, rout, scale
    integer :: k

    xa = 0.0_wp
    call species_metal_sums(nH, ne, 0.0_wp, nH, nem, rin, rout)
    Ra = -(A + rin)
    xb = 1.0_wp
    call species_metal_sums(nH, ne, nH, 0.0_wp, nem, rin, rout)
    Rb = B + rout
    scale = A + B + abs(Ra) + abs(Rb)

    x = xb
    do k = 1, 80
       if (Rb == Ra) exit
       x = xb - Rb*(xb - xa)/(Rb - Ra)
       !--- keep the trial strictly inside the bracket
       x = min(max(x, min(xa, xb)), max(xa, xb))
       call species_metal_sums(nH, ne, nH*x, nH*(1.0_wp - x), nem, rin, rout)
       Rx = x*(B + rin + rout) - (A + rin)
       if (abs(Rx) <= 1.0e-15_wp*scale) exit
       if (abs(xb - xa) <= 1.0e-15_wp) exit
       if (Rx*Rb < 0.0_wp) then
          xa = xb;  Ra = Rb
       else
          Ra = 0.5_wp*Ra          ! Illinois: halve the stale endpoint
       end if
       xb = x;  Rb = Rx
    end do
  end function hydrogen_neutral_fraction

  !=========================================================================
  ! Electron density of a state whose H/He fractions are already fixed, from
  ! charge neutrality including the metal cascade (par%metal_ne):
  !     n_e = n_e(H,He) + sum_el n_H A_el sum_i (i-1) f_i(n_e)
  ! The metal stage fractions f_i depend on n_e themselves (recombination
  ! goes as n_e alpha), so this is a scalar fixed point, solved with the
  ! same damped iteration and tolerance solve_ion_cell uses for its own
  ! n_e closure.  A single evaluation of the metal term at the H/He-only
  ! electron density is NOT the solution: beyond the ionization front the
  ! H/He term vanishes, so that argument is the n_e -> 0 limit in which
  ! nothing can recombine and every low-IP metal comes out fully ionized.
  !
  ! Its purpose is the WRITTEN state: the under-relaxed H/He fractions are
  ! not the ones solve_ion_cell converged with, and the stage fractions
  ! written alongside (species_write) are evaluated at the written n_e and
  ! T_e — so n_e must satisfy charge neutrality at exactly those fractions
  ! or the output describes two different gases.
  !
  ! ne_seed is the electron density the cell solve converged to; starting
  ! there keeps the iteration on the branch the physics solve selected.
  ! T is the temperature the state is WRITTEN at (par%te_fixed, or the
  ! newly solved T_e), which is not in general the last trial temperature
  ! of the thermal bisection, so the (leaf, T) rate cache is refilled here
  ! rather than inherited.
  !=========================================================================
  real(kind=wp) function charge_neutrality_ne(il, T, nH, ne_HHe, xHI, ne_seed) &
                                             result(ne)
    use species_mod, only : species_ne_prepare, species_ne_cached
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, ne_HHe, xHI, ne_seed
    real(kind=wp) :: ne_old, nHI, nHII
    integer :: it

    call species_ne_prepare(il, T)
    nHI  = nH*xHI
    nHII = nH*(1.0_wp - xHI)
    ne   = max(ne_seed, 1.0e-12_wp*nH)
    do it = 1, 200
       ne_old = ne
       ne = ne_HHe + species_ne_cached(nH, ne_old, nHI, nHII)
       ne = max(0.5_wp*(ne + ne_old), 1.0e-12_wp*nH)
       if (abs(ne - ne_old) <= 1.0e-10_wp*ne) exit
    end do
  end function charge_neutrality_ne

  !=========================================================================
  ! Relative residual of hydrogen's ionization balance on the WRITTEN state
  ! of leaf il, including both charge-exchange terms:
  !
  !   eps = | n_HI (Gamma_H + C_H n_e + R_cx_out) - (alpha_H n_e + R_cx_in) n_HII |
  !         / [ n_HI (Gamma_H + C_H n_e + R_cx_out) ]
  !
  ! Inside solve_ion_cell the two sides of every charge-exchange pair are
  ! formed from the same cached coefficients and the same stage densities, so
  ! the pair balances by construction there.  What this measures is the state
  ! actually written out, which differs from the solved one because the H/He
  ! fractions are under-relaxed toward it (par%ion_relax) while the metal
  ! stage fractions written alongside are evaluated at the written n_e and T.
  !
  ! With par%charge_exchange = 'metal_only' the metal R_cx terms are dropped
  ! here as well, and the He-H terms follow par%helium_charge_exchange, so
  ! the residual always measures the equation the run actually solves.
  !
  ! Domain of validity: the state carries x_HI, so n_HII = n_H (1 - x_HI) has
  ! the relative precision eps/(1 - x_HI) and the residual cannot resolve
  ! anything below that.  In a fully neutral leaf 1 - x_HI rounds to zero and
  ! eps saturates at 1 whatever the physics does; such leaves are excluded
  ! from the reported maximum and mean by XHII_FLOOR below rather than being
  ! allowed to set them.
  !=========================================================================
  real(kind=wp) function hydrogen_balance_residual(il, T, nH, xHI, ne) result(eps)
    use species_mod, only : species_ne_prepare, species_metal_sums, n_elements
    integer,       intent(in) :: il
    real(kind=wp), intent(in) :: T, nH, xHI, ne
    real(kind=wp) :: aH, cH, nHI, nHII, nem, r_cx_in, r_cx_out, loss, gain
    real(kind=wp) :: xHeIII, dgHI, dgHeI, dh1, dh2, dh3

    r_cx_in  = 0.0_wp
    r_cx_out = 0.0_wp
    if (par%use_metals .and. n_elements > 0 .and. &
        trim(par%charge_exchange) == 'two_way') then
       call species_ne_prepare(il, T)
       call species_metal_sums(nH, ne, nH*xHI, nH*(1.0_wp - xHI), &
                               nem, r_cx_in, r_cx_out)
    end if
    if (trim(par%case_ab) == 'A') then
       aH = alphaA_HII(T)
    else
       aH = alphaB_HII(T)
    end if
    cH   = ci_HI(T)
    nHI  = nH*xHI
    nHII = nH*(1.0_wp - xHI)
    !--- He-H charge exchange on the written He fractions
    if (par%helium_charge_exchange .and. par%abund_He > 0.0_wp) then
       xHeIII = max(0.0_wp, 1.0_wp - gas_xHeI(il) - gas_xHeII(il))
       r_cx_out = r_cx_out + par%abund_He*nH*(gas_xHeII(il)*cx_HeII_H(T) &
                                              + xHeIII*cx_HeIII_H(T))
       r_cx_in  = r_cx_in  + par%abund_He*nH*(gas_xHeI(il)*cx_HeI_Hp(T) &
                                              + gas_xHeII(il)*cx_HeII_Hp(T))
    end if
    !--- secondary ionization at the written state (par%use_sec_ion)
    dgHI = 0.0_wp
    if (par%use_sec_ion) call secondary_ionization(il, &
         ionized_fraction_svs(ne, nH, par%abund_He*nH), nHI, &
         par%abund_He*nH*gas_xHeI(il), par%abund_He*nH*gas_xHeII(il), &
         dgHI, dgHeI, dh1, dh2, dh3)
    loss = nHI*(gamma_HI(il) + dgHI + cH*ne + r_cx_out)
    gain = (aH*ne + r_cx_in)*nHII
    eps  = abs(loss - gain) / max(loss, tinest)
  end function hydrogen_balance_residual

  !=========================================================================
  ! The same residual over the grid: its maximum over leaves and its
  ! volume-weighted mean, over the leaves where it resolves anything
  ! (x_HII above XHII_FLOOR, so n_HII = n_H (1 - x_HI) keeps at least six
  ! significant digits).  Diagnostic only — it is printed, not acted on.
  !=========================================================================
  subroutine hydrogen_balance_residual_grid(eps_max, eps_vol)
    use mpi
    use octree_mod, only : leaf_half
    implicit none
    real(kind=wp), intent(out) :: eps_max, eps_vol
    real(kind=wp), parameter :: XHII_FLOOR = 1.0e-10_wp
    real(kind=wp) :: eps, vol, sum_ev, sum_v
    integer :: il, ierr

    eps_max = 0.0_wp
    eps_vol = 0.0_wp
    !--- the gas state is node-shared and every node's rank 0 holds the same
    !--- values, so one rank per node evaluates and broadcasts.
    if (mpar%h_rank == 0) then
       sum_ev = 0.0_wp
       sum_v  = 0.0_wp
       do il = 1, gas_nleaf
          if (gas_nH(il) <= 0.0_wp) cycle
          if (1.0_wp - gas_xHI(il) < XHII_FLOOR) cycle
          eps = hydrogen_balance_residual(il, gas_Te(il), gas_nH(il), &
                                          gas_xHI(il), gas_ne(il))
          eps_max = max(eps_max, eps)
          vol = (2.0_wp*leaf_half(il))**3
          sum_ev = sum_ev + vol*eps
          sum_v  = sum_v  + vol
       end do
       eps_vol = sum_ev / max(sum_v, tinest)
    end if
    call MPI_BCAST(eps_max, 1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
    call MPI_BCAST(eps_vol, 1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
  end subroutine hydrogen_balance_residual_grid

  !=========================================================================
  ! Returns both convergence measures (par%conv_crit picks which gates the
  ! iteration): max_dx = max leaf |delta x_HII| (stalls at the front-cell
  ! Monte Carlo noise floor); dx_vol = sum V |delta x_HII| / sum V x_HII,
  ! the volume-weighted L1 change of the ionized volume — front jitter
  ! occupies little volume, so this measure converges when the global
  ! solution has.  dhe_vol(2): the same measure for helium, reported and
  ! not gating (helium_change_measures).
  subroutine gas_equilibrium_update(max_dx, dx_vol, dhe_vol)
    use mpi
    use octree_mod, only : amr_grid, leaf_half
    implicit none
    real(kind=wp), intent(out) :: max_dx, dx_vol, dhe_vol(2)

    real(kind=wp), allocatable :: xHI_new(:), xHeI_new(:), xHeII_new(:), ne_new(:)
    real(kind=wp) :: T, yHe, w, vol, sum_dxv, sum_xv, sum_he(4)
    real(kind=wp) :: nH, ne, xHI, xHeI, xHeII, xHeIII
    logical :: caseA
    integer :: il, ierr

    T   = par%te_fixed
    yHe = par%abund_He
    w   = par%ion_relax
    caseA = trim(par%case_ab) == 'A'

    !--- only the node-local rank 0 solves.  jt_ion is ALLREDUCEd over
    !--- MPI_COMM_WORLD, so every node's rank 0 sees the identical tally and
    !--- produces identical values; the other node ranks skip the (expensive)
    !--- solve for each leaf and receive the state through shared memory + the
    !--- convergence metrics by broadcast.
    if (mpar%h_rank == 0) then
    allocate(xHI_new(gas_nleaf), xHeI_new(gas_nleaf), xHeII_new(gas_nleaf), &
             ne_new(gas_nleaf))

    max_dx  = 0.0_wp
    sum_dxv = 0.0_wp
    sum_xv  = 0.0_wp
    sum_he  = 0.0_wp
    do il = 1, gas_nleaf
       nH = gas_nH(il)
       if (nH <= 0.0_wp) then
          xHI_new(il) = gas_xHI(il);  xHeI_new(il) = gas_xHeI(il)
          xHeII_new(il) = gas_xHeII(il);  ne_new(il) = 0.0_wp
          cycle
       end if

       xHI = gas_xHI(il);  xHeI = gas_xHeI(il);  xHeII = gas_xHeII(il)
       call solve_ion_cell(gamma_HI(il), gamma_HeI(il), gamma_HeII(il), &
                           nH, T, caseA, xHI, xHeI, xHeII, ne, il)

       !--- under-relaxation toward the solved state.
       xHI_new(il)   = (1.0_wp - w)*gas_xHI(il)   + w*xHI
       xHeI_new(il)  = (1.0_wp - w)*gas_xHeI(il)  + w*xHeI
       xHeII_new(il) = (1.0_wp - w)*gas_xHeII(il) + w*xHeII
       xHeIII = max(0.0_wp, 1.0_wp - xHeI_new(il) - xHeII_new(il))
       ne_new(il) = nH * ((1.0_wp - xHI_new(il)) &
                    + yHe*(xHeII_new(il) + 2.0_wp*xHeIII))
       !--- metal electrons in the written state (par%metal_ne), solved
       !--- self-consistently: the stage fractions depend on n_e, so the
       !--- closure has to be iterated at the under-relaxed H/He state, not
       !--- evaluated once at the H/He-only electron density.
       if (par%metal_ne .and. par%use_metals) then
          block
            use species_mod, only : n_elements
            if (n_elements > 0) ne_new(il) = &
               charge_neutrality_ne(il, T, nH, ne_new(il), xHI_new(il), ne)
          end block
       end if
       max_dx = max(max_dx, abs(xHI_new(il) - gas_xHI(il)))
       vol = (2.0_wp*leaf_half(il))**3
       sum_dxv = sum_dxv + vol*abs(xHI_new(il) - gas_xHI(il))
       sum_xv  = sum_xv  + vol*(1.0_wp - xHI_new(il))
       sum_he(1) = sum_he(1) + vol*abs(xHeI_new(il) - gas_xHeI(il))
       sum_he(2) = sum_he(2) + vol*(1.0_wp - xHeI_new(il))
       sum_he(3) = sum_he(3) + vol*abs(xHeIII &
                 - max(0.0_wp, 1.0_wp - gas_xHeI(il) - gas_xHeII(il)))
       sum_he(4) = sum_he(4) + vol*xHeIII
    end do
    dx_vol = sum_dxv / max(sum_xv, tinest)
    dhe_vol = helium_change_measures(sum_he)

    !--- write the solved state straight into the node-shared arrays.
    do il = 1, gas_nleaf
       gas_xHI(il)   = xHI_new(il)
       gas_xHeI(il)  = xHeI_new(il)
       gas_xHeII(il) = xHeII_new(il)
       gas_ne(il)    = ne_new(il)
    end do
    deallocate(xHI_new, xHeI_new, xHeII_new, ne_new)
    end if   ! h_rank == 0

    !--- broadcast the metrics to the node's other ranks and barrier so the
    !--- shared-memory writes are visible before this returns.
    call MPI_BCAST(max_dx, 1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
    call MPI_BCAST(dx_vol, 1, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
    call MPI_BCAST(dhe_vol, 2, MPI_DOUBLE_PRECISION, 0, mpar%hostcomm, ierr)
    call MPI_BARRIER(mpar%hostcomm, ierr)
  end subroutine gas_equilibrium_update

  !=========================================================================
  ! Volume-integrated change of the helium state over one iteration, the
  ! helium counterpart of dx_vol, from the sums
  !   s = [sum V |d x_HeI|, sum V (1 - x_HeI), sum V |d x_HeIII|,
  !        sum V x_HeIII]
  ! (x at the new state):
  !   d(1) = s1/s2, the He0 <-> He+ front, normalized on ionized helium so
  !          that He+ <-> He2+ exchange does not enter it;
  !   d(2) = s3/s4, the He+ <-> He2+ front, set to -1 (not measured) where
  !          He2+ is less than HE3_SHARE of the ionized helium.  Below that
  !          its relative change is Monte Carlo noise on a trace species:
  !          in the HII40 benchmark He2+ fills 1.2e-4 of the volume and this
  !          ratio sits at 5-9e-2 per iteration for 80 iterations, so it
  !          would never pass a stopping test while carrying nothing the
  !          model depends on.
  !=========================================================================
  pure function helium_change_measures(s) result(d)
    real(kind=wp), intent(in) :: s(4)
    real(kind=wp) :: d(2)
    real(kind=wp), parameter :: HE3_SHARE = 1.0e-2_wp
    d(1) = s(1)/max(s(2), tinest)
    if (s(4) > HE3_SHARE*s(2)) then
       d(2) = s(3)/s(4)
    else
       d(2) = -1.0_wp
    end if
  end function helium_change_measures

end module ion_balance_mod
