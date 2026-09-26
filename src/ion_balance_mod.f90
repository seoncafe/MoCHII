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
! alpha_B when the ground-recombination photons are absorbed on the spot
! (par%recombination_continuum = 'on_the_spot', the default), alpha_A
! when they are transported or escape.
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
                            gas_Te, gas_nleaf, gas_state_publish, &
                            gas_leaf_owner_first, gas_leaf_owner_stride
  use gas_rates_mod, only : gamma_HI, gamma_HeI, gamma_HeII, &
                            secondary_ionization
  use photoelectron_degradation_mod, only : ionized_fraction_svs
  use recomb_mod
  implicit none
  private

  public :: gas_equilibrium_update, solve_ion_cell, charge_neutrality_ne
  public :: helium_change_measures, signed_change_measures
  public :: hydrogen_balance_residual, hydrogen_balance_residual_grid
  public :: cx_HeII_H, cx_HeI_Hp, cx_HeIII_H, cx_HeII_Hp
  public :: helium_charge_exchange_heating, nonradiative_HeII_H

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
  !
  ! With par%hei_metastable_coupling (and il present) the He I 2^3S
  ! population n_3 = n_e [alpha_3 n_HeII + (q_13 + q_feed) n_HeI] / D
  ! (hei_metastable_mod) enters both balances at the current iterate:
  !  - helium: ionization of 2^3S by electrons (C_3) and photons (Phi_3)
  !    returns the atom to He+, adding (n_e C_3 + Phi_3) n_3 to He0 -> He+.
  !    n_3 is linear in n_HeI and n_HeII, so this is solved JOINTLY in the
  !    closed form, not lagged: with u = (n_e C_3 + Phi_3)/D, the share of
  !    2^3S atoms that end ionized,
  !      n_HeII/n_HeI = [Gamma + C n_e + u (q_13 + q_feed) n_e (+ CX)]
  !                     / [(alpha - u alpha_3) n_e (+ CX)],
  !    i.e. a share u of the triplet recombinations is undone and a share u
  !    of the collisional triplet excitations becomes a two-step ionization.
  !    alpha - u alpha_3 > 0 always (alpha_3 = f_3S alpha_B < alpha_B <= alpha).
  !    D is evaluated at this iterate's n_e and n_HI.
  !  - hydrogen: Penning + associative ionization He(2^3S) + H0 adds n_3 Q_31
  !    to the H0 ionization rate; n_3 there uses the previous iterate's
  !    helium fractions (the lag rhe_* already carry), so the exit test then
  !    also requires x_HI and x_HeII to have stopped moving.
  subroutine solve_ion_cell(gH, gHe1, gHe2, nH, T, caseA, &
                            xHI, xHeI, xHeII, ne, il, f_h_only)
    use species_mod, only : species_ne_prepare, species_metal_sums, n_elements
    use hei_metastable_mod, only : hei_triplet_rates_type, hei_triplet_rates, &
                                   hei_triplet_destruction, &
                                   hei_triplet_population, hei_phi3_rate
    use hei_recombination_photons_mod, only : hei_ground_nodes_type, &
                                   hei_ground_nodes, hei_on_the_spot_sources
    use heii_recombination_photons_mod, only : heii_ots_nodes_type, &
                                   heii_ots_nodes, heii_on_the_spot_sources
    real(kind=wp), intent(in)    :: gH, gHe1, gHe2, nH, T
    logical,       intent(in)    :: caseA
    real(kind=wp), intent(inout) :: xHI, xHeI, xHeII
    real(kind=wp), intent(out)   :: ne
    integer,       intent(in), optional :: il
    !--- fraction of the on-the-spot photons only H0 can absorb (the He I
    !--- excited cascade, the He II cascade below the He0 edge) that H0 did
    !--- absorb: 1 unless they outnumber hydrogen's recombinations
    real(kind=wp), intent(out), optional :: f_h_only
    real(kind=wp) :: aH, aHe2, aHe3, cH, cHe1, cHe2, yHe
    real(kind=wp) :: ne_old, xHeIII, r1, r2, denH
    real(kind=wp) :: nem, r_cx_in, r_cx_out
    real(kind=wp) :: k21, k12, k32, k23, rhe_in, rhe_out, xHI_old, xHeII_old
    real(kind=wp) :: gHt, gHe1t, dgHI, dgHeI, dh1, dh2, dh3
    logical :: with_metal_ne, with_metal_cx, with_he_cx, with_sec, settled
    logical :: with_meta, with_ots, with_heii_ots
    type(hei_triplet_rates_type) :: hr
    type(hei_ground_nodes_type)  :: gnodes
    type(heii_ots_nodes_type)    :: hnodes
    real(kind=wp) :: aHe3B, yH3, yHe3, sH3, sHe3, num2, den2, sH3o, sfix, afix
    real(kind=wp) :: aHb, bHb, s0fix, sabs, supply
    real(kind=wp) :: aHe2B, a1He, y_ots, s_gnd, s_exc
    real(kind=wp) :: phi3, r_pen, u3, num1, den1, gHe1b
    logical :: joint
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
    !--- He I 2^3S coupling: the rate set at this T once, Phi_3 of the leaf
    with_meta = par%hei_metastable_coupling .and. present(il) .and. yHe > 0.0_wp
    if (with_meta) then
       hr   = hei_triplet_rates(T)
       phi3 = hei_phi3_rate(il)
    end if
    r_pen = 0.0_wp
    !--- case-B on-the-spot He I recombination photons (par%hei_on_the_spot,
    !--- hei_recombination_photons_mod): the excited cascade and the share y
    !--- of the ground captures that H0 absorbs ionize hydrogen; helium keeps
    !--- alpha_B + y alpha_1.  The share y falls with n(H0), so the ground
    !--- captures join hydrogen's photoionization rate per H0 at the previous
    !--- iterate.  The excited photons lie below the He0 edge and have no
    !--- other absorber: their number does not fall with n(H0), so they are
    !--- an absolute source sfix (per n_H) in hydrogen's balance, which
    !--- x_HI = 0 bounds where they outnumber hydrogen's recombinations
    !--- (heii_recombination_photons_mod, "Absorption limit").
    with_ots = par%hei_on_the_spot .and. (.not. caseA) .and. present(il) &
               .and. yHe > 0.0_wp
    if (with_ots) then
       gnodes = hei_ground_nodes(T)
       aHe2B  = aHe2
       a1He   = alphaA_HeII(T) - aHe2B
    end if
    !--- case-B on-the-spot He III recombination photons (par%heii_on_the_spot,
    !--- heii_recombination_photons_mod): the ground-capture shares y_H, y_He
    !--- and the excited cascade ionize H0 and He0.  Both sources join the
    !--- photoionization rates per absorber at the previous iterate, as the
    !--- He I ground captures do, except H0's share of the excited cascade,
    !--- which joins sfix (heii_recombination_photons_mod, "Absorption
    !--- limit"); He III keeps alpha_B + (y_H + y_He) alpha_1.  As a
    !--- rate per He0 the source stays finite as n(He0) -> 0 (its share falls
    !--- with it); a source count held at the previous share
    !--- could exceed the He+ -> He0 recombinations that replenish He0 and
    !--- leave the closed form without a positive solution.  In the He2+
    !--- zone, where H0 and He0 are both trace absorbers, even the lagged
    !--- rate per He0 let the two shares trade places from pass to pass, so
    !--- there the two neutral fractions are solved together
    !--- (joint_neutral_fractions, below).
    with_heii_ots = par%heii_on_the_spot .and. (.not. caseA) .and. &
                    present(il) .and. yHe > 0.0_wp
    if (with_heii_ots) then
       hnodes = heii_ots_nodes(T)
       aHe3B  = aHe3
    end if

    xHeIII = max(0.0_wp, 1.0_wp - xHeI - xHeII)
    ne = nH * ((1.0_wp - xHI) + yHe*(xHeII + 2.0_wp*xHeIII))
    ne = max(ne, 1.0e-12_wp*nH)
    nem = 0.0_wp;  r_cx_in = 0.0_wp;  r_cx_out = 0.0_wp
    sfix = 0.0_wp;  afix = 0.0_wp;  s0fix = 0.0_wp
    do it = 1, 200
       ne_old = ne
       xHI_old = xHI;  xHeII_old = xHeII
       if (with_sec) then
          call secondary_ionization(il, ionized_fraction_svs(ne, nH, yHe*nH), &
               nH*xHI, nH*yHe*xHeI, nH*yHe*xHeII, dgHI, dgHeI, dh1, dh2, dh3)
          gHt = gH + dgHI;  gHe1t = gHe1 + dgHeI
       end if
       if (with_ots) then
          call hei_on_the_spot_sources(il, T, nH, ne, xHI, xHeI, xHeII, gnodes, &
                                       y_ots, s_gnd, s_exc)
          gHt  = gH + s_gnd/(nH*max(xHI, 1.0e-12_wp))
          if (with_sec) gHt = gHt + dgHI
          aHe2 = aHe2B + y_ots*a1He
          sfix = s_exc/nH
       end if
       joint = .false.
       if (with_heii_ots) then
          call heii_on_the_spot_sources(hnodes, ne, nH*xHI, yHe*nH*xHeI, &
               yHe*nH*xHeII, yHe*nH*xHeIII, yH3, yHe3, sH3, sHe3, sH3o)
          if (.not. (with_sec .or. with_ots)) gHt = gH
          if (.not. with_sec) gHe1t = gHe1
          gHe1b = gHe1t
          gHt   = gHt   + sH3/(nH*max(xHI, 1.0e-12_wp))
          gHe1t = gHe1t + sHe3/(yHe*nH*max(xHeI, 1.0e-12_wp))
          aHe3  = aHe3B + (yH3 + yHe3)*hnodes%a1
          !--- in the He2+ zone H0 and He0 are both trace absorbers of the
          !--- same cascade, and their shares are solved together
          !--- (joint_neutral_fractions) instead of He0's being lagged
          joint = xHeIII > 1.0e-3_wp
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
       if (with_meta) then
          !--- Penning ionization of H by He(2^3S), n_3 from this iterate's
          !--- n_e, x_HI and the previous iterate's helium (lagged, like the
          !--- rhe terms, which are zero here without helium charge exchange);
          !--- an ionization of H0, so it joins Gamma_H.
          r_pen = hr%pen*hei_triplet_population(hr, ne, nH*xHI, &
                     yHe*nH*xHeI, yHe*nH*xHeII, phi3)
       end if
       if (with_ots .or. with_heii_ots) then
          !--- the on-the-spot excited cascades as an absolute source: the
          !--- same balance as the branches below, x B = A - s(x), with H0's
          !--- share of the He II cascade evaluated at every trial x
          !--- (hydrogen_fraction_with_cascade)
          if (with_meta .and. with_metal_cx) then
             aHb = aH*ne + rhe_in
             bHb = gHt + (cH + aH)*ne + rhe_in + rhe_out + r_pen
          else if (with_meta) then
             aHb = aH*ne + r_cx_in + rhe_in
             bHb = gHt + (cH + aH)*ne + r_cx_in + r_cx_out + rhe_in &
                   + rhe_out + r_pen
          else if (with_metal_cx .and. with_he_cx) then
             aHb = aH*ne + rhe_in
             bHb = gHt + (cH + aH)*ne + rhe_in + rhe_out
          else if (with_metal_cx) then
             aHb = aH*ne
             bHb = gHt + (cH + aH)*ne
          else if (with_he_cx) then
             aHb = aH*ne + r_cx_in + rhe_in
             bHb = gHt + (cH + aH)*ne + r_cx_in + r_cx_out + rhe_in + rhe_out
          else
             aHb = aH*ne + r_cx_in
             bHb = gHt + (cH + aH)*ne + r_cx_in + r_cx_out
          end if
          if (joint) then
             call joint_neutral_fractions()
          else
             xHI = hydrogen_fraction_with_cascade(nH, ne, aHb, bHb, sfix, &
                      with_heii_ots, hnodes, yHe*nH*xHeI, yHe*nH*xHeIII, &
                      with_metal_cx, afix, s0fix)
          end if
       else if (with_meta) then
          if (with_metal_cx) then
             xHI = hydrogen_neutral_fraction(nH, ne, aH*ne + rhe_in, &
                      gHt + (cH + aH)*ne + rhe_in + rhe_out + r_pen)
          else
             denH = gHt + (cH + aH)*ne + r_cx_in + r_cx_out + rhe_in &
                    + rhe_out + r_pen
             xHI = (aH*ne + r_cx_in + rhe_in) / denH    ! denH >= aH ne > 0
          end if
       else if (with_metal_cx .and. with_he_cx) then
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
       if (with_heii_ots) then
          !--- the same balances as the branches below, assembled term by term
          !--- so that the He III photons (already in gHe1t and aHe3) leave the
          !--- other paths' arithmetic untouched
          num1 = gHe1t + cHe1*ne
          den1 = aHe2*ne
          num2 = gHe2 + cHe2*ne
          den2 = aHe3*ne
          if (with_meta) then
             u3   = (ne*hr%c3 + phi3)/hei_triplet_destruction(hr, ne, nH*xHI, phi3)
             num1 = num1 + u3*(hr%q13 + hr%qfeed)*ne
             den1 = den1 - u3*hr%alpha3*ne
          end if
          if (with_he_cx) then
             num1 = num1 + k12*nH*(1.0_wp - xHI)
             den1 = den1 + k21*nH*xHI
             num2 = num2 + k23*nH*(1.0_wp - xHI)
             den2 = den2 + k32*nH*xHI
          end if
          r1 = min(num1/den1, 1.0e150_wp)
          r2 = min(num2/den2, 1.0e150_wp)
       else if (with_meta) then
          !--- He0 -> He+ through 2^3S, solved jointly (header comment)
          u3   = (ne*hr%c3 + phi3)/hei_triplet_destruction(hr, ne, nH*xHI, phi3)
          num1 = gHe1t + cHe1*ne + u3*(hr%q13 + hr%qfeed)*ne
          den1 = (aHe2 - u3*hr%alpha3)*ne
          if (with_he_cx) then
             num1 = num1 + k12*nH*(1.0_wp - xHI)
             den1 = den1 + k21*nH*xHI
             r2 = min((gHe2 + cHe2*ne + k23*nH*(1.0_wp - xHI)) &
                      / (aHe3*ne + k32*nH*xHI), 1.0e150_wp)
          else
             r2 = min((gHe2 + cHe2*ne) / (aHe3*ne), 1.0e150_wp)
          end if
          r1 = min(num1/den1, 1.0e150_wp)
       else if (with_he_cx) then
          r1 = min((gHe1t + cHe1*ne + k12*nH*(1.0_wp - xHI)) &
                   / (aHe2*ne + k21*nH*xHI), 1.0e150_wp)
          r2 = min((gHe2 + cHe2*ne + k23*nH*(1.0_wp - xHI)) &
                   / (aHe3*ne + k32*nH*xHI), 1.0e150_wp)
       else
          r1 = min((gHe1t + cHe1*ne) / (aHe2*ne), 1.0e150_wp)
          r2 = min((gHe2 + cHe2*ne) / (aHe3*ne), 1.0e150_wp)
       end if
       if (.not. joint) then
          xHeI   = 1.0_wp / (1.0_wp + r1 + r1*r2)
          xHeII  = xHeI * r1
          xHeIII = xHeII * r2
       end if
       ne = nH * ((1.0_wp - xHI) + yHe*(xHeII + 2.0_wp*xHeIII))
       if (with_metal_ne) ne = ne + nem
       ne = max(0.5_wp*(ne + ne_old), 1.0e-12_wp*nH)
       settled = .true.
       if (with_he_cx .or. with_sec .or. with_meta .or. with_ots &
           .or. with_heii_ots) settled = &
          abs(xHI - xHI_old) <= 1.0e-10_wp*min(xHI, 1.0_wp - xHI) + 1.0e-15_wp &
          .and. abs(xHeII - xHeII_old) <= 1.0e-10_wp*xHeII + 1.0e-15_wp
       if (abs(ne - ne_old) <= 1.0e-10_wp*ne .and. settled) exit
    end do
    !--- the share of the on-the-spot photons given to H0 that it can absorb
    !--- at the state reached: no more than the recombinations (and charge
    !--- exchanges) that make H0.  1 wherever the balance closed; below 1
    !--- only where the helium cascades outnumber that supply, where no local
    !--- balance exists (heii_recombination_photons_mod, "Absorption limit")
    !--- and the photons beyond it are dropped from the heat too.
    if (present(f_h_only)) then
       f_h_only = 1.0_wp
       if (with_ots .or. with_heii_ots) then
          sabs = 0.0_wp
          if (with_ots) then
             call hei_on_the_spot_sources(il, T, nH, ne, xHI, xHeI, xHeII, &
                                          gnodes, y_ots, s_gnd, s_exc)
             sabs = s_gnd + s_exc
          end if
          if (with_heii_ots) then
             xHeIII = max(0.0_wp, 1.0_wp - xHeI - xHeII)
             call heii_on_the_spot_sources(hnodes, ne, nH*xHI, yHe*nH*xHeI, &
                  yHe*nH*xHeII, yHe*nH*xHeIII, yH3, yHe3, sH3, sHe3, sH3o)
             sabs = sabs + sH3 + sH3o
          end if
          supply = (aH*ne + r_cx_in + rhe_in)*nH*(1.0_wp - xHI)
          if (sabs > supply) f_h_only = max(supply, 0.0_wp)/sabs
       end if
    end if

  contains

    !=======================================================================
    ! H0 and He0 in the He2+ zone, where both are trace absorbers of the
    ! He II on-the-spot cascade (the parts above the He0 edge are shared by
    ! opacity, n sigma) and of the He2+ ground captures (shared with He+).
    ! With He0's share entering its balance as a rate per He0 at the
    ! previous iterate, the two shares traded places from pass to pass and
    ! a cell could leave the loop with hydrogen's balance open by O(1)
    ! (He/H = 0.1: 0.89 in the He III core of the G2a run; He/H = 0.2 in
    ! tests/heii_on_the_spot (g)).  Here, at this pass's n_e and He2+/He+
    ! ratio r2, x_HeI is the root of helium's balance
    !
    !    F(x) = N1 x + s_He(x_HI(x), x)/(y n_H) - D1 (1 - x)/(1 + r2) ,
    !
    ! N1 the ionization rate of He0 per He0 without the He II photons, D1
    ! the recombination rate of He+ per He+, s_He the He0 share of the He II
    ! photons at the trial state, and x_HI(x) hydrogen's own root at that
    ! n(He0) (hydrogen_fraction_with_cascade).  F(0) < 0 because the He0
    ! share vanishes with n(He0), and F rises with x (He0's share grows,
    ! H0's falls), so the root is unique; Illinois in ln x, which can be
    ! 1e-12.  The charge-exchange and 2^3S terms of N1 and D1 are taken at
    ! the trial x_HI, the He2+/He+ ratio at the pass's (lagged) state.
    !=======================================================================
    subroutine joint_neutral_fractions()
      real(kind=wp) :: ta, tb, tx, Fa, Fb, Fx, scale
      integer :: k
      r2 = min(helium_upper_ratio(), 1.0e150_wp)
      ta = log(1.0e-30_wp);  Fa = helium_residual(ta)
      tb = 0.0_wp;           Fb = helium_residual(tb)
      if (Fb <= 0.0_wp) then
         tx = tb
      else
         scale = abs(Fa) + abs(Fb)
         tx = tb
         do k = 1, 100
            if (Fb == Fa) exit
            tx = tb - Fb*(tb - ta)/(Fb - Fa)
            tx = min(max(tx, min(ta, tb)), max(ta, tb))
            Fx = helium_residual(tx)
            if (abs(Fx) <= 1.0e-14_wp*scale) exit
            if (abs(tb - ta) <= 1.0e-12_wp) exit
            if (Fx*Fb < 0.0_wp) then
               ta = tb;  Fa = Fb
            else
               Fa = 0.5_wp*Fa         ! Illinois
            end if
            tb = tx;  Fb = Fx
         end do
      end if
      Fx = helium_residual(tx)        ! leaves x_HI and helium at the root
    end subroutine joint_neutral_fractions

    !--- He2+/He+ at the pass's state (the same terms as the r2 of the
    !--- closed forms below)
    real(kind=wp) function helium_upper_ratio() result(r)
      r = (gHe2 + cHe2*ne)/(aHe3*ne)
      if (with_he_cx) r = (gHe2 + cHe2*ne + k23*nH*(1.0_wp - xHI)) &
                          /(aHe3*ne + k32*nH*xHI)
    end function helium_upper_ratio

    !--- F at x_HeI = exp(lnx); sets xHI, xHeI, xHeII, xHeIII to that state
    real(kind=wp) function helium_residual(lnx) result(F)
      real(kind=wp), intent(in) :: lnx
      real(kind=wp) :: n1, d1, y_h, y_he, s_h, s_he, s_ho
      xHeI   = exp(lnx)
      xHeII  = (1.0_wp - xHeI)/(1.0_wp + r2)
      xHeIII = xHeII*r2
      xHI = hydrogen_fraction_with_cascade(nH, ne, aHb, bHb, sfix, &
               .true., hnodes, yHe*nH*xHeI, yHe*nH*xHeIII, &
               with_metal_cx, afix, s0fix)
      call heii_on_the_spot_sources(hnodes, ne, nH*xHI, yHe*nH*xHeI, &
           yHe*nH*xHeII, yHe*nH*xHeIII, y_h, y_he, s_h, s_he, s_ho)
      n1 = gHe1b + cHe1*ne
      d1 = aHe2*ne
      if (with_meta) then
         u3 = (ne*hr%c3 + phi3)/hei_triplet_destruction(hr, ne, nH*xHI, phi3)
         n1 = n1 + u3*(hr%q13 + hr%qfeed)*ne
         d1 = d1 - u3*hr%alpha3*ne
      end if
      if (with_he_cx) then
         n1 = n1 + k12*nH*(1.0_wp - xHI)
         d1 = d1 + k21*nH*xHI
      end if
      F = n1*xHeI + s_he/(yHe*nH) - d1*xHeII
    end function helium_residual

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
  ! Charge-exchange heating of the four He-H reactions [erg cm^-3 s^-1]
  ! (par%charge_exchange_heating with par%helium_charge_exchange): each
  ! reaction's rate times the kinetic energy one event gives the gas, the
  ! ionization-potential difference less what leaves as radiation.
  !  - He+ + H -> He + H+, non-radiative: Zygelman et al. (1989, PRA 40,
  !    2340; read) couple only the A 1Sigma+ entrance and X 1Sigma+ exit
  !    states, which separate to He(1s2 1S) + H+, so the whole
  !    IP(He) - IP(H) = 10.989 eV goes to the gas (He 2 3S at 19.82 eV is
  !    out of reach).
  !  - He+ + H -> He + H+ + photon, radiative: the photon carries the
  !    defect (Zygelman et al.: emission peaked near 113 nm, the A-X
  !    separation at large R); nothing to the gas.
  !  - He + H+ -> He+ + H, the detailed-balance reverse of the non-radiative
  !    channel: takes the 10.989 eV back, so the pair deposits nothing when
  !    the two balance.
  !  - He++ + H -> He+ + H+: the KF96 rate is their 1e-14 cm^3 s^-1 value,
  !    which KF96 (Sec. 3.1) sets for reactions that proceed by radiative
  !    charge transfer; the photon then carries the defect.  Capture into
  !    He+(n=2), resonant with H(1s), would give the gas only the 0.006 eV
  !    reduced-mass difference and radiate the rest as He II 304 A.  West,
  !    Lane & Cohen (1982, PRA 26, 3164; read), the paper KF96 cite, find
  !    the radiative channel (He+(1s) + H+ + photon) dominant at low
  !    energy.  APPROXIMATION, marked: the heat of that channel is set to
  !    zero, but it is not zero.  The photon is emitted at an internuclear
  !    distance R and carries the electronic separation Delta E(R) of the
  !    2p sigma and 1s sigma states; the products (He+ + H+, which repel)
  !    keep Delta E(infinity) - Delta E(R) as kinetic energy.  Their Fig. 2
  !    has Delta E(R) at its minimum, ~1.15 a.u., near R = 2.7 a0 and
  !    ~1.4 a.u. at 10 a0, against 1.5 a.u. (40.8 eV) at infinity, so the
  !    gas keeps ~3-10 eV per event for emission anywhere in that range.
  !    They give no emission spectrum (their method "bypasses a
  !    determination of the emission spectrum"), so no mean is available
  !    to adopt.  The same holds for the radiative X2+ + H rows of the
  !    metals (make_element_data.py).  Where He2+ meets H0 at all these
  !    reactions are a small term of the heat.
  !  - He+ + H+ -> He++ + H, cx_HeII_Hp: from ground-state reactants the
  !    reaction must take the whole IP(He+) - IP(H) = 40.819 eV from the
  !    gas.  That rate is the detailed-balance reverse of the ground channel
  !    of He++ + H, which is not the channel the forward rate describes, so
  !    this pair does not cancel; at 1e4 K the rate is ~1e-35 cm^3 s^-1.
  !=========================================================================
  !--- (He++ + H gives nothing, so n(He++) does not enter.)
  real(kind=wp) function helium_charge_exchange_heating(T, nHI, nHII, nHeI, &
                                                        nHeII) result(heat)
    real(kind=wp), intent(in) :: T, nHI, nHII, nHeI, nHeII
    heat = ( (nonradiative_HeII_H(T)*nHeII*nHI - cx_HeI_Hp(T)*nHeI*nHII)*DE_HEI_H &
           - cx_HeII_Hp(T)*nHeII*nHII*DE_HEII_H ) * ev2erg
  end function helium_charge_exchange_heating

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
  !=========================================================================
  ! Hydrogen's neutral fraction with the on-the-spot excited cascades as an
  ! absolute source (solve_ion_cell):
  !    x B + x (rin + rout) = A + rin - s(x),     s = s1 + s2(x),
  ! per n_H, where s1 is the He I cascade (lagged: it moves only through
  ! p_584) and s2(x) = H0's share of the He II cascade at n(H0) = n_H x,
  ! evaluated at every trial because that share turns on steeply as n(H0)
  ! passes n(He0) in the He2+ zone, where a lagged share oscillates.  rin
  ! and rout are the metal charge-exchange sums at the trial x when
  ! metal_cx (as hydrogen_neutral_fraction), else zero (the closed-form
  ! branches carry them in A and B).  The residual rises monotonically in
  ! x (s2 does), so the root is unique; where the source exceeds the H0
  ! supply of a fully ionized cell, F(0) >= 0, there is none in (0,1] and
  ! x = 0: the absorption limit (heii_recombination_photons_mod).  supply0
  ! is A + rin at x = 0 and s0 the source there, for the absorbed fraction.
  !=========================================================================
  real(kind=wp) function hydrogen_fraction_with_cascade(nH, ne, A, B, s1, &
       with_heii, g, nHeI, nHeIII, metal_cx, supply0, s0) result(x)
    use species_mod, only : species_metal_sums
    use heii_recombination_photons_mod, only : heii_ots_nodes_type, &
                                               heii_cascade_to_H
    real(kind=wp),             intent(in)  :: nH, ne, A, B, s1, nHeI, nHeIII
    logical,                   intent(in)  :: with_heii, metal_cx
    type(heii_ots_nodes_type), intent(in)  :: g
    real(kind=wp),             intent(out) :: supply0, s0
    real(kind=wp) :: xa, xb, Ra, Rb, Rx, scale
    integer :: k

    xa = 0.0_wp
    Ra = residual(xa, supply0, s0)
    if (Ra >= 0.0_wp) then
       x = 0.0_wp
       return
    end if
    xb = 1.0_wp
    Rb = residual(xb)
    scale = A + B + abs(Ra) + abs(Rb)
    x = xb
    do k = 1, 100
       if (Rb == Ra) exit
       x = xb - Rb*(xb - xa)/(Rb - Ra)
       x = min(max(x, min(xa, xb)), max(xa, xb))
       Rx = residual(x)
       if (abs(Rx) <= 1.0e-15_wp*scale) exit
       if (abs(xb - xa) <= 1.0e-15_wp) exit
       if (Rx*Rb < 0.0_wp) then
          xa = xb;  Ra = Rb
       else
          Ra = 0.5_wp*Ra          ! Illinois: halve the stale endpoint
       end if
       xb = x;  Rb = Rx
    end do

  contains

    real(kind=wp) function residual(xt, sup, sx) result(r)
      real(kind=wp), intent(in)            :: xt
      real(kind=wp), intent(out), optional :: sup, sx
      real(kind=wp) :: nem, rin, rout, s
      rin = 0.0_wp;  rout = 0.0_wp
      if (metal_cx) call species_metal_sums(nH, ne, nH*xt, nH*(1.0_wp - xt), &
                                            nem, rin, rout)
      s = s1
      if (with_heii) s = s + heii_cascade_to_H(g, ne, nH*xt, nHeI, nHeIII)/nH
      r = xt*(B + rin + rout) - (A + rin) + s
      if (present(sup)) sup = A + rin
      if (present(sx))  sx  = s
    end function residual
  end function hydrogen_fraction_with_cascade

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
  ! here as well, the He-H terms follow par%helium_charge_exchange, and the
  ! Penning term n_3 Q_31 follows par%hei_metastable_coupling, so the
  ! residual always measures the equation the run actually solves.
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
    real(kind=wp) :: xHeIII, dgHI, dgHeI, dh1, dh2, dh3, s_only

    r_cx_in  = 0.0_wp
    r_cx_out = 0.0_wp
    if (par%use_metals .and. n_elements > 0 .and. &
        trim(par%charge_exchange) == 'two_way') then
       call species_ne_prepare(il, T)
       call species_metal_sums(nH, ne, nH*xHI, nH*(1.0_wp - xHI), &
                               nem, r_cx_in, r_cx_out)
    end if
    if (trim(par%recombination_continuum) /= 'on_the_spot') then
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
    !--- H0's share of the on-the-spot excited cascades
    s_only = 0.0_wp
    !--- Penning ionization by He(2^3S) on the written state
    !--- (par%hei_metastable_coupling)
    if (par%hei_metastable_coupling .and. par%abund_He > 0.0_wp) then
       block
         use hei_metastable_mod, only : hei_triplet_rates_type, &
              hei_triplet_rates, hei_triplet_population, hei_phi3_rate
         type(hei_triplet_rates_type) :: hr
         hr = hei_triplet_rates(T)
         loss = loss + nHI*hr%pen*hei_triplet_population(hr, ne, nHI, &
                par%abund_He*nH*gas_xHeI(il), par%abund_He*nH*gas_xHeII(il), &
                hei_phi3_rate(il))
       end block
    end if
    !--- on-the-spot He I recombination photons (par%hei_on_the_spot)
    if (par%hei_on_the_spot .and. &
        trim(par%recombination_continuum) == 'on_the_spot' .and. &
        par%abund_He > 0.0_wp) then
       block
         use hei_recombination_photons_mod, only : hei_ground_nodes, &
              hei_on_the_spot_sources
         real(kind=wp) :: y_ots, s_gnd, s_exc
         call hei_on_the_spot_sources(il, T, nH, ne, xHI, gas_xHeI(il), &
              gas_xHeII(il), hei_ground_nodes(T), y_ots, s_gnd, s_exc)
         loss = loss + s_gnd
         s_only = s_only + s_exc
       end block
    end if
    !--- on-the-spot He III recombination photons (par%heii_on_the_spot)
    if (par%heii_on_the_spot .and. &
        trim(par%recombination_continuum) == 'on_the_spot' .and. &
        par%abund_He > 0.0_wp) then
       block
         use heii_recombination_photons_mod, only : heii_ots_nodes, &
              heii_on_the_spot_sources
         real(kind=wp) :: yH3, yHe3, sH3, sHe3, sH3o, nHe
         nHe = par%abund_He*nH
         call heii_on_the_spot_sources(heii_ots_nodes(T), ne, nHI, &
              nHe*gas_xHeI(il), nHe*gas_xHeII(il), &
              nHe*max(0.0_wp, 1.0_wp - gas_xHeI(il) - gas_xHeII(il)), &
              yH3, yHe3, sH3, sHe3, sH3o)
         loss = loss + sH3
         s_only = s_only + sH3o
       end block
    end if
    !--- counted in full: where the helium cascades outnumber hydrogen's
    !--- recombinations no local balance exists, and this residual is what
    !--- shows it (heii_recombination_photons_mod, "Absorption limit")
    gain = (aH*ne + r_cx_in)*nHII
    loss = loss + s_only
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
    real(kind=wp) :: eps, vol, sums(2)
    integer :: il, ierr

    !--- every rank evaluates its own leaves (gas_state_mod, "Where the cell
    !--- solves run"); the sums and the maximum are reduced over
    !--- MPI_COMM_WORLD and global rank 0's values broadcast
    eps_max = 0.0_wp
    sums = 0.0_wp
    do il = gas_leaf_owner_first(), gas_nleaf, gas_leaf_owner_stride()
       if (gas_nH(il) <= 0.0_wp) cycle
       if (1.0_wp - gas_xHI(il) < XHII_FLOOR) cycle
       eps = hydrogen_balance_residual(il, gas_Te(il), gas_nH(il), &
                                       gas_xHI(il), gas_ne(il))
       eps_max = max(eps_max, eps)
       vol = (2.0_wp*leaf_half(il))**3
       sums(1) = sums(1) + vol*eps
       sums(2) = sums(2) + vol
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE, eps_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, sums, 2, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    eps_vol = sums(1) / max(sums(2), tinest)
    call MPI_BCAST(eps_max, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(eps_vol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  end subroutine hydrogen_balance_residual_grid

  !=========================================================================
  ! Returns both convergence measures (par%conv_crit picks which gates the
  ! iteration): max_dx = max leaf |delta x_HII| (stalls at the front-cell
  ! Monte Carlo noise floor); dx_vol = sum V |delta x_HII| / sum V x_HII,
  ! the volume-weighted L1 change of the ionized volume — front jitter
  ! occupies little volume, so this measure converges when the global
  ! solution has.  dhe_vol(2): the same measure for helium
  ! (helium_change_measures); with conv_crit = 'vol' it enters the stopping
  ! test beside dx_vol (main.f90::composition_converged).  drift_vol(3):
  ! the signed counterparts (signed_change_measures).
  subroutine gas_equilibrium_update(max_dx, dx_vol, dhe_vol, drift_vol)
    use mpi
    use octree_mod, only : amr_grid, leaf_half
    implicit none
    real(kind=wp), intent(out) :: max_dx, dx_vol, dhe_vol(2), drift_vol(3)

    real(kind=wp), allocatable :: xHI_new(:), xHeI_new(:), xHeII_new(:), ne_new(:)
    real(kind=wp) :: T, yHe, w, vol, sums(9), xHeIII_old
    real(kind=wp) :: nH, ne, xHI, xHeI, xHeII, xHeIII
    logical :: caseA
    integer :: il, ierr

    T   = par%te_fixed
    yHe = par%abund_He
    w   = par%ion_relax
    caseA = trim(par%recombination_continuum) /= 'on_the_spot'

    !--- every rank solves its own leaves (gas_state_mod, "Where the cell
    !--- solves run"; until 2026-09-26 each node's rank 0 solved all of them
    !--- alone), into arrays that are zero elsewhere, and gas_state_publish
    !--- sums them over MPI_COMM_WORLD into every node's shared state
    allocate(xHI_new(gas_nleaf), xHeI_new(gas_nleaf), xHeII_new(gas_nleaf), &
             ne_new(gas_nleaf))
    xHI_new = 0.0_wp;  xHeI_new = 0.0_wp;  xHeII_new = 0.0_wp;  ne_new = 0.0_wp

    max_dx  = 0.0_wp
    sums    = 0.0_wp
    do il = gas_leaf_owner_first(), gas_nleaf, gas_leaf_owner_stride()
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
       sums(1) = sums(1) + vol*abs(xHI_new(il) - gas_xHI(il))
       sums(2) = sums(2) + vol*(1.0_wp - xHI_new(il))
       sums(3) = sums(3) + vol*abs(xHeI_new(il) - gas_xHeI(il))
       sums(4) = sums(4) + vol*(1.0_wp - xHeI_new(il))
       xHeIII_old = max(0.0_wp, 1.0_wp - gas_xHeI(il) - gas_xHeII(il))
       sums(5) = sums(5) + vol*abs(xHeIII - xHeIII_old)
       sums(6) = sums(6) + vol*xHeIII
       sums(7) = sums(7) + vol*(gas_xHI(il) - xHI_new(il))
       sums(8) = sums(8) + vol*(gas_xHeI(il) - xHeI_new(il))
       sums(9) = sums(9) + vol*(xHeIII - xHeIII_old)
    end do
    !--- the metrics compare with the state before this update, so they are
    !--- reduced before the new state goes into the shared arrays
    call MPI_ALLREDUCE(MPI_IN_PLACE, max_dx, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, sums, 9, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    call gas_state_publish(xHI_new, xHeI_new, xHeII_new, ne_new)
    deallocate(xHI_new, xHeI_new, xHeII_new, ne_new)
    dx_vol  = sums(1) / max(sums(2), tinest)
    dhe_vol = helium_change_measures(sums(3:6))
    drift_vol = signed_change_measures(sums([7, 2, 8, 4, 9, 6]), dhe_vol)
    !--- global rank 0's metrics everywhere, so every rank makes the same
    !--- loop-exit decision
    call MPI_BCAST(max_dx, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(dx_vol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(dhe_vol, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(drift_vol, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
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

  !=========================================================================
  ! Signed counterparts of dx_vol and helium_change_measures, from
  !   s = [sum V (x_HII,new - x_HII,old), sum V x_HII,
  !        sum V (x_HeI,old - x_HeI,new), sum V (1 - x_HeI),
  !        sum V (x_HeIII,new - x_HeIII,old), sum V x_HeIII]
  ! (denominators at the new state, as there): the net change of the
  ! ionized volume of each front over one iteration, positive when it
  ! grows.  The third is 0 where the He+ <-> He2+ step is not measured
  ! (dhe(2) < 0).  fixed_point_distance_mod uses them to tell a front that
  ! still moves from one fluctuating about its fixed point.
  !=========================================================================
  pure function signed_change_measures(s, dhe) result(d)
    real(kind=wp), intent(in) :: s(6), dhe(2)
    real(kind=wp) :: d(3)
    d(1) = s(1)/max(s(2), tinest)
    d(2) = s(3)/max(s(4), tinest)
    if (dhe(2) >= 0.0_wp) then
       d(3) = s(5)/max(s(6), tinest)
    else
       d(3) = 0.0_wp
    end if
  end function signed_change_measures

end module ion_balance_mod
