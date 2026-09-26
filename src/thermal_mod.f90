module thermal_mod
!---------------------------------------------------------------------------
! MoCHII: thermal balance — solve T_e with the ionization state.
!
! Per leaf, find T_e such that photoheating balances cooling:
!     H(T) = n_HI H_HI + n_HeI H_HeI + n_HeII H_HeII     (rate integrals,
!            per particle, from gas_rates_mod — T-independent)
!     C(T) = cooling_total(T, ...)                        (cooling_mod)
! where the ionization fractions at each trial T come from solve_ion_cell
! (photo + collisional vs recombination at that T).  Since C rises and the
! ionization response is mild, the net function H - C is monotone through
! the root: bisection on log T in [par%te_min, par%te_max] (~40 iterations
! to 1e-4 relative).
!
! gas_thermal_update replaces gas_equilibrium_update when par%solve_te:
! it writes x, n_e AND T_e (under-relaxation on x) and returns
! max |delta x_HII| and max |delta T_e|/T_e.
!---------------------------------------------------------------------------
  use define
  use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII, gas_ne, &
                            gas_Te, gas_nleaf, gas_state_sync, &
                            gas_state_publish, gas_leaf_owner_first, &
                            gas_leaf_owner_stride
  use gas_rates_mod, only : gamma_HI, gamma_HeI, gamma_HeII, &
                            heat_HI, heat_HeI, heat_HeII, secondary_ionization
  use photoelectron_degradation_mod, only : ionized_fraction_svs
  use ion_balance_mod, only : solve_ion_cell, charge_neutrality_ne, &
                              helium_change_measures, signed_change_measures, &
                              helium_charge_exchange_heating
  use cooling_mod,     only : cooling_total
  implicit none
  private

  public :: gas_thermal_update

contains

  !=========================================================================
  ! Net heating - cooling [erg cm^-3 s^-1] at trial temperature T for one
  ! leaf; the ionization state at T is solved self-consistently.
  !=========================================================================
  subroutine net_rate(il, T, caseA, xHI, xHeI, xHeII, ne, net, &
                      heat_star, heat_heii_ots)
    integer,       intent(in)    :: il
    real(kind=wp), intent(in)    :: T
    logical,       intent(in)    :: caseA
    real(kind=wp), intent(inout) :: xHI, xHeI, xHeII
    real(kind=wp), intent(out)   :: ne, net
    !--- two terms of the heating, for gas_thermal_update's on-the-spot
    !--- validity count: the photoheating by the transported field and the
    !--- heat of the on-the-spot He II cascade [erg cm^-3 s^-1]
    real(kind=wp), intent(out), optional :: heat_star, heat_heii_ots
    real(kind=wp) :: nH, heat, yHe, dgHI, dgHeI, dhHI, dhHeI, dhHeII, f_h_only

    nH = gas_nH(il)
    yHe = par%abund_He
    call solve_ion_cell(gamma_HI(il), gamma_HeI(il), gamma_HeII(il), &
                        nH, T, caseA, xHI, xHeI, xHeII, ne, il, f_h_only)
    !--- photoheating; with par%use_sec_ion the part of each fast
    !--- photoelectron's energy that ionizes or excites instead of heating is
    !--- taken off, at the state the solve just reached.
    dhHI = 0.0_wp;  dhHeI = 0.0_wp;  dhHeII = 0.0_wp
    if (par%use_sec_ion) call secondary_ionization(il, &
         ionized_fraction_svs(ne, nH, yHe*nH), nH*xHI, nH*yHe*xHeI, &
         nH*yHe*xHeII, dgHI, dgHeI, dhHI, dhHeI, dhHeII)
    heat = nH*( xHI*(heat_HI(il) - dhHI) &
           + yHe*(xHeI*(heat_HeI(il) - dhHeI) + xHeII*(heat_HeII(il) - dhHeII)) )
    net = heat - cooling_total(T, nH, ne, xHI, xHeI, xHeII)
    if (present(heat_star)) heat_star = heat
    if (present(heat_heii_ots)) heat_heii_ots = 0.0_wp
    !--- trace-metal line cooling (registry)
    if (par%use_metals) then
       block
         use species_mod, only : metal_cooling, metal_heating, metal_freefree, &
                                 metal_cooling_H, metal_ionization_cooling, &
                                 metal_charge_exchange_heating
         net = net - metal_cooling(il, T, nH, ne, nH*xHI, nH*(1.0_wp-xHI))
         !--- atomic-hydrogen impact on the C I / C II / O I ground-term
         !--- fine structure: the increment over the electron-only
         !--- populations of the table above, from one statistical
         !--- equilibrium with both colliders; the dominant coolant wherever
         !--- the gas is neutral, whatever heats it.
         net = net - metal_cooling_H(il, T, nH, ne, nH*xHI, nH*(1.0_wp-xHI))
         !--- metal free-free cooling, net-charge weighted (par%metal_freefree)
         if (par%metal_freefree) &
            net = net - metal_freefree(il, T, nH, ne, nH*xHI, nH*(1.0_wp-xHI))
         !--- metal ionization energy budget (par%metal_heat): photoheating
         !--- in; recombination kinetic energy and collisional-ionization
         !--- energy out.
         if (par%metal_heat) &
            net = net + metal_heating(il, T, nH, ne, nH*xHI, nH*(1.0_wp-xHI)) &
                  - metal_ionization_cooling(il, T, nH, ne, nH*xHI, nH*(1.0_wp-xHI))
         !--- kinetic energy exchanged by the metal-H charge-exchange
         !--- reactions (par%charge_exchange_heating)
         if (par%charge_exchange_heating) &
            net = net + metal_charge_exchange_heating(il, T, nH, ne, &
                                                      nH*xHI, nH*(1.0_wp-xHI))
       end block
    end if
    !--- the same for the He-H reactions of par%helium_charge_exchange
    if (par%charge_exchange_heating .and. par%helium_charge_exchange &
        .and. yHe > 0.0_wp) &
       net = net + helium_charge_exchange_heating(T, nH*xHI, nH*(1.0_wp-xHI), &
                                                  yHe*nH*xHeI, yHe*nH*xHeII)
    !--- He I 2^3S processes (par%hei_metastable_coupling): excitation into
    !--- the triplets, the collisional transfers out of 2^3S, its electron-
    !--- impact ionization, its photoionization heating and Penning
    !--- ionization of H, on the state just solved (the energy ledger is
    !--- derived in hei_metastable_mod's header).
    if (par%hei_metastable_coupling .and. yHe > 0.0_wp) then
       block
         use hei_metastable_mod, only : hei_triplet_rates_type, &
              hei_triplet_rates, hei_triplet_population, &
              hei_triplet_net_heating, hei_phi3_rate, hei_phi3_heating
         type(hei_triplet_rates_type) :: hr
         real(kind=wp) :: n3
         hr = hei_triplet_rates(T)
         n3 = hei_triplet_population(hr, ne, nH*xHI, nH*yHe*xHeI, &
                                     nH*yHe*xHeII, hei_phi3_rate(il))
         net = net + hei_triplet_net_heating(hr, ne, nH*xHI, nH*yHe*xHeI, &
                                             n3, hei_phi3_heating(il))
       end block
    end if
    !--- case-B on-the-spot He I recombination photons absorbed by H0
    !--- (par%hei_on_the_spot; the heat of each absorption is derived in
    !--- hei_recombination_photons_mod's header).  f_h_only is the share of
    !--- the photons only H0 can take that it did take (solve_ion_cell):
    !--- photons beyond its absorption limit heat nothing.
    if (par%hei_on_the_spot .and. (.not. caseA) .and. yHe > 0.0_wp) then
       block
         use hei_recombination_photons_mod, only : hei_ground_nodes, &
              hei_on_the_spot_heating
         net = net + hei_on_the_spot_heating(il, T, nH, ne, xHI, xHeI, xHeII, &
                                             hei_ground_nodes(T), f_h_only)
       end block
    end if
    !--- case-B on-the-spot He III recombination photons absorbed by H0 and
    !--- He0 (par%heii_on_the_spot; the heat of each absorption is derived in
    !--- heii_recombination_photons_mod's header)
    if (par%heii_on_the_spot .and. (.not. caseA) .and. yHe > 0.0_wp) then
       block
         use heii_recombination_photons_mod, only : heii_ots_nodes, &
              heii_on_the_spot_heating
         real(kind=wp) :: h2
         h2 = heii_on_the_spot_heating(heii_ots_nodes(T), ne, nH*xHI, &
                     yHe*nH*xHeI, yHe*nH*xHeII, &
                     yHe*nH*max(0.0_wp, 1.0_wp - xHeI - xHeII), f_h_only)
         net = net + h2
         if (present(heat_heii_ots)) heat_heii_ots = h2
       end block
    end if
    !--- PDR: grain photoelectric heating + grain recombination cooling
    !--- (Bakes & Tielens 1994 fits) driven by the local FUV field.
    !--- Two guards, both found by the PDR gate (runaway to 12-35 kK
    !--- without them): (i) the BT94 fits are evaluated at
    !--- min(T, 1e4 K) — outside their fit range the T^0.7 term of
    !--- epsilon diverges and heating runs away through the H
    !--- collisional-ionization feedback; (ii) both terms carry an
    !--- x_HI weight — the BT94 charging parameter knows only the FUV
    !--- field, so in ionized gas (EUV-charged grains, and the grain
    !--- energy already routed to the IR via Heat_dust) it would badly
    !--- overestimate; grain photoelectric GAS heating belongs to the
    !--- predominantly neutral zone.
    if (par%grain_pe) then
       block
         use gas_rates_mod, only : g0_fuv
         use octree_mod,    only : amr_grid
         !--- MW reference dust extinction per H (D03 V band) — the
         !--- normalization of the BT94 fits.
         real(kind=wp), parameter :: CEXT_MW_V = 4.868e-22_wp
         real(kind=wp) :: g0, xpe, eps, beta, dscale, Tpe
         g0 = g0_fuv(il)
         if (g0 > 0.0_wp) then
            !--- dust amount relative to MW from the leaf opacity itself
            !--- (robust against the DGR-vs-cext_dust input conventions;
            !--- the first gate run had DGR=1 with the MW extinction
            !--- already folded into cext_dust — the naive DGR/0.01
            !--- scale made the PE heating 100x too strong and drove
            !--- the PDR zone to a false 12 kK Ly-alpha-capped
            !--- equilibrium).
            dscale = par%pe_scale * xHI &
                     * (amr_grid%rhokap(il)/par%distance2cm) &
                     / (max(nH, tinest)*CEXT_MW_V)
            Tpe = min(T, 1.0e4_wp)
            xpe = g0*sqrt(Tpe)/max(ne, 1.0e-12_wp*nH)
            eps = 4.87e-2_wp/(1.0_wp + 4.0e-3_wp*xpe**0.73_wp) &
                + 3.65e-2_wp*(Tpe/1.0e4_wp)**0.7_wp/(1.0_wp + 2.0e-4_wp*xpe)
            net = net + 1.3e-24_wp*eps*g0*nH*dscale
            beta = 0.74_wp/Tpe**0.068_wp
            net = net - 4.65e-30_wp*Tpe**0.94_wp*xpe**beta*ne*nH*dscale
         end if
       end block
    end if
  end subroutine net_rate

  !=========================================================================
  ! Returns cell-max AND volume-integrated convergence measures
  ! (par%conv_crit picks which pair gates the iteration):
  !   dx_vol  = sum V |delta x_HII| / sum V x_HII
  !   dte_vol = sum V n_e |delta T_e| / sum V n_e T_e
  ! The n_e weight removes the cells with vanishing heating (deep PDR /
  ! vacuum), whose te_min-pinned oscillation dominates max_dte.
  ! dhe_vol(2), the same measure for helium, which enters the 'vol'
  ! stopping test beside dx_vol (main.f90::composition_converged):
  !   dhe_vol(1) = sum V |delta x_HeI| / sum V (1 - x_HeI)    He0 <-> He+
  !   dhe_vol(2) = sum V |delta x_HeIII| / sum V x_HeIII      He+ <-> He2+
  ! (ion_balance_mod::helium_change_measures; dhe_vol(2) = -1 where He2+
  ! fills a negligible share of the volume).  drift_vol(3) are the signed
  ! counterparts of dx_vol, dhe_vol(1), dhe_vol(2) (the net change of
  ! x_HII, 1 - x_HeI, x_HeIII over the same denominators), which tell a
  ! moving front from the Monte Carlo floor (fixed_point_distance_mod).
  !
  ! Two validity counts are printed when not empty: the leaves whose thermal
  ! balance has no root below par%te_max (heating wins there, T is held at
  ! te_max), and, in the on-the-spot mode, the leaves where the He II
  ! on-the-spot photons (par%heii_on_the_spot: the He2+ ground captures and,
  ! unless par%heii_diffuse transports it, the cascade) out-heat the
  ! transported field in gas optically thin at 40.8 eV (tau across the leaf
  ! < 0.01).
  ! There the cascade photons in fact travel on, and depositing their
  ! E - IP (up to ~27 eV each) in the leaf is not what on-the-spot
  ! absorption assumes: in an optically thin He III core it pins T at
  ! te_max (the G2a run, 2026-09-26).  The mode is left as defined;
  ! par%heii_diffuse = .true. transports that cascade alone (channel 5,
  ! heii_cascade_mod::heii_cascade_transported), and 'transported' carries
  ! every recombination photon.
  subroutine gas_thermal_update(max_dx, max_dte, dx_vol, dte_vol, dhe_vol, &
                                drift_vol)
    use mpi
    use octree_mod, only : leaf_half
    use photo_xsec, only : sigma_HI, sigma_HeI
    implicit none
    real(kind=wp), intent(out) :: max_dx, max_dte, dx_vol, dte_vol, dhe_vol(2)
    real(kind=wp), intent(out) :: drift_vol(3)

    real(kind=wp), allocatable :: xHI_new(:), xHeI_new(:), xHeII_new(:), &
                                  ne_new(:), te_new(:)
    real(kind=wp) :: nH, ne, xHI, xHeI, xHeII, xHeIII, w, vol
    real(kind=wp) :: max_local(2), max_host(2)
    real(kind=wp) :: sum_local(11), sum_host(11), metrics(9), xHeIII_old
    real(kind=wp) :: tlo, thi, tmid, net_lo, net_hi, net_mid, te
    real(kind=wp) :: h_star, h_ots, tau_leaf, flag_local(5), flag_host(5)
    real(kind=wp), parameter :: E_HEII_LYA = 40.813_wp, TAU_THIN = 1.0e-2_wp
    logical :: caseA, at_te_max
    integer :: il, it, ierr

    caseA = trim(par%recombination_continuum) /= 'on_the_spot'
    w     = par%ion_relax
    allocate(xHI_new(gas_nleaf), xHeI_new(gas_nleaf), xHeII_new(gas_nleaf), &
             ne_new(gas_nleaf), te_new(gas_nleaf))
    !--- zero off this rank's leaves: gas_state_publish sums over the ranks
    xHI_new = 0.0_wp;  xHeI_new = 0.0_wp;  xHeII_new = 0.0_wp
    ne_new  = 0.0_wp;  te_new   = 0.0_wp

    max_local = 0.0_wp
    sum_local = 0.0_wp
    !--- validity counts: gas volume; leaves and volume at te_max; leaves and
    !--- volume where the on-the-spot He II cascade out-heats the field in
    !--- optically thin gas
    flag_local = 0.0_wp
    !--- Cyclic ownership over every rank of MPI_COMM_WORLD (gas_state_mod,
    !--- "Where the cell solves run") mixes cheap boundary-clamped cells with
    !--- expensive bracketed cells, and adding nodes divides the work.
    do il = gas_leaf_owner_first(), gas_nleaf, gas_leaf_owner_stride()
       nH = gas_nH(il)
       if (nH <= 0.0_wp) then
          xHI_new(il) = gas_xHI(il);  xHeI_new(il) = gas_xHeI(il)
          xHeII_new(il) = gas_xHeII(il);  ne_new(il) = 0.0_wp
          te_new(il) = gas_Te(il)
          cycle
       end if

       !--- bisection on log T; ionization re-solved at each trial T.
       xHI = gas_xHI(il);  xHeI = gas_xHeI(il);  xHeII = gas_xHeII(il)
       tlo = par%te_min;  thi = par%te_max
       call net_rate(il, tlo, caseA, xHI, xHeI, xHeII, ne, net_lo)
       call net_rate(il, thi, caseA, xHI, xHeI, xHeII, ne, net_hi)
       at_te_max = .false.
       if (net_lo <= 0.0_wp) then
          te = tlo               ! cooling wins even at te_min (no heating)
       else if (net_hi >= 0.0_wp) then
          te = thi               ! heating wins even at te_max
          at_te_max = .true.
       else
          do it = 1, 60
             tmid = sqrt(tlo*thi)
             call net_rate(il, tmid, caseA, xHI, xHeI, xHeII, ne, net_mid)
             if (net_mid > 0.0_wp) then
                tlo = tmid
             else
                thi = tmid
             end if
             if (thi/tlo - 1.0_wp < 1.0e-5_wp) exit
          end do
          te = sqrt(tlo*thi)
       end if
       !--- final consistent state at te.
       call net_rate(il, te, caseA, xHI, xHeI, xHeII, ne, net_mid, &
                     h_star, h_ots)
       vol = (2.0_wp*leaf_half(il))**3
       flag_local(1) = flag_local(1) + vol
       if (at_te_max) then
          flag_local(2) = flag_local(2) + 1.0_wp
          flag_local(3) = flag_local(3) + vol
       end if
       if (h_ots > h_star) then
          tau_leaf = 2.0_wp*leaf_half(il)*par%distance2cm*nH &
                     *(xHI*sigma_HI(E_HEII_LYA) &
                       + par%abund_He*xHeI*sigma_HeI(E_HEII_LYA))
          if (tau_leaf < TAU_THIN) then
             flag_local(4) = flag_local(4) + 1.0_wp
             flag_local(5) = flag_local(5) + vol
          end if
       end if

       !--- under-relaxation on x; te taken directly.
       xHI_new(il)   = (1.0_wp - w)*gas_xHI(il)   + w*xHI
       xHeI_new(il)  = (1.0_wp - w)*gas_xHeI(il)  + w*xHeI
       xHeII_new(il) = (1.0_wp - w)*gas_xHeII(il) + w*xHeII
       xHeIII = max(0.0_wp, 1.0_wp - xHeI_new(il) - xHeII_new(il))
       ne_new(il) = nH * ((1.0_wp - xHI_new(il)) &
                    + par%abund_He*(xHeII_new(il) + 2.0_wp*xHeIII))
       !--- metal electrons in the WRITTEN state too (the solve already
       !--- had them in its closure; the state n_e drives the n-level
       !--- excitation and the outputs — in the PDR zone it IS the
       !--- metal contribution, n_e ~ A_C n_H).  Charge neutrality is
       !--- re-solved at the under-relaxed H/He fractions and at the newly
       !--- solved te (the temperature the state is written at), seeded
       !--- from the n_e the cell solve converged to.
       if (par%metal_ne .and. par%use_metals) then
          block
            use species_mod, only : n_elements
            if (n_elements > 0) ne_new(il) = &
               charge_neutrality_ne(il, te, nH, ne_new(il), xHI_new(il), ne)
          end block
       end if
       te_new(il) = te
       max_local(1) = max(max_local(1), abs(xHI_new(il) - gas_xHI(il)))
       max_local(2) = max(max_local(2), &
                          abs(te_new(il) - gas_Te(il))/gas_Te(il))
       sum_local(1) = sum_local(1) + vol*abs(xHI_new(il) - gas_xHI(il))
       sum_local(2) = sum_local(2) + vol*(1.0_wp - xHI_new(il))
       sum_local(3) = sum_local(3) &
                    + vol*ne_new(il)*abs(te_new(il) - gas_Te(il))
       sum_local(4) = sum_local(4) + vol*ne_new(il)*te_new(il)
       sum_local(5) = sum_local(5) + vol*abs(xHeI_new(il) - gas_xHeI(il))
       sum_local(6) = sum_local(6) + vol*(1.0_wp - xHeI_new(il))
       xHeIII_old = max(0.0_wp, 1.0_wp - gas_xHeI(il) - gas_xHeII(il))
       sum_local(7) = sum_local(7) + vol*abs(xHeIII - xHeIII_old)
       sum_local(8) = sum_local(8) + vol*xHeIII
       sum_local(9)  = sum_local(9)  + vol*(gas_xHI(il) - xHI_new(il))
       sum_local(10) = sum_local(10) + vol*(gas_xHeI(il) - xHeI_new(il))
       sum_local(11) = sum_local(11) + vol*(xHeIII - xHeIII_old)
    end do

    call MPI_ALLREDUCE(max_local, max_host, 2, MPI_DOUBLE_PRECISION, MPI_MAX, &
                       MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(sum_local, sum_host, 11, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(flag_local, flag_host, 5, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    if (mpar%p_rank == 0 .and. flag_host(1) > 0.0_wp) then
       if (flag_host(2) > 0.0_wp) write(*,'(a,i0,a,es9.2,a,f9.1,a)') &
          '                  NOTE: ', nint(flag_host(2)), ' leaves (', &
          flag_host(3)/flag_host(1), ' of the gas volume) at te_max = ', &
          par%te_max, ' K: heating exceeds cooling there'
       if (flag_host(4) > 0.0_wp) then
          write(*,'(a,i0,a,es9.2,a)') &
             '                  NOTE: ', nint(flag_host(4)), ' leaves (', &
             flag_host(5)/flag_host(1), ' of the gas volume) where the on-the-spot'// &
             ' He II recombination photons out-heat the field in gas thin at'// &
             ' 40.8 eV (tau_leaf < 0.01): they leave those leaves, which'// &
             ' on-the-spot absorption does not describe'
          if (par%heii_diffuse) then
             write(*,'(a)') '                        (the cascade is'// &
                ' transported; these are the He2+ ground captures, which'// &
                ' recombination_continuum = ''transported'' also carries)'
          else
             write(*,'(a)') '                        (heii_diffuse = .true.'// &
                ' transports the cascade, recombination_continuum ='// &
                ' ''transported'' every recombination photon)'
          end if
       end if
    end if
    max_dx  = max_host(1)
    max_dte = max_host(2)
    dx_vol  = sum_host(1) / max(sum_host(2), tinest)
    dte_vol = sum_host(3) / max(sum_host(4), tinest)
    dhe_vol = helium_change_measures(sum_host(5:8))
    drift_vol = signed_change_measures(sum_host([9, 2, 10, 6, 11, 8]), dhe_vol)

    !--- every rank's leaves into every node's shared state
    call gas_state_publish(xHI_new, xHeI_new, xHeII_new, ne_new, te_new)
    deallocate(xHI_new, xHeI_new, xHeII_new, ne_new, te_new)

    !--- The metrics are reduced over MPI_COMM_WORLD; global rank 0's copy is
    !--- broadcast so that every rank makes the same loop-exit decision
    !--- whatever the reduction algorithm returns on each.
    metrics = [max_dx, max_dte, dx_vol, dte_vol, dhe_vol(1), dhe_vol(2), &
               drift_vol(1), drift_vol(2), drift_vol(3)]
    call MPI_BCAST(metrics, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    max_dx  = metrics(1)
    max_dte = metrics(2)
    dx_vol  = metrics(3)
    dte_vol = metrics(4)
    dhe_vol = metrics(5:6)
    drift_vol = metrics(7:9)
  end subroutine gas_thermal_update

end module thermal_mod
