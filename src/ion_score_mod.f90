module ion_score_mod
!---------------------------------------------------------------------------
! Common ionizing-packet path scorer.
!
! Every direct analytic segment and every sampled scattered segment supplies
! the same path quantity
!
!   path_lum = Lpacket * weight * effective_path_length
!
! [erg/s * code length].  Every rate coefficient is scored at the packet's own
! photon%energy_eV; those reduced arrays are authoritative for gas_rates_mod,
! while the jt_ion(inu,leaf) tally is only a diagnostic spectrum.
!
! Under par%diffuse_launch = 'absorption_reemission' (plan
! md/ABSORPTION_REEMISSION_PLAN_2026-11-13.md, Section 12, r7):
! - the photoionizations of each re-emitting absorber s = H0, He0, He+ of
!   leaf i in the field of the last reduced pass with the densities of the
!   CURRENT state,
!     N_abs,s = V n_s Gamma_s = D n_s Sum(path_lum sigma_s/(E erg))   [1/s]
!   (previous_photoionizations), from the raw H/He rate path sums, which
!   diffuse_build reads before ion_score_reset clears them;
! - the gas absorption measured in the pass on its own state,
!     A_i = Sum path_lum * kap_gas(E, i)   [erg/s],
!   kap_gas = ion_packet_mod::ion_packet_gas_opacity, rank local, reduced
!   over MPI_COMM_WORLD by ion_score_reduce and published node-shared as
!   gas_absorbed_power (GasAbsorbedPower): a diagnostic.
! A new grid (ion_score_setup/resize) clears the path sums and the flag
! that marks them valid, so its first pass has no previous field.
!---------------------------------------------------------------------------
  use utility, only : abort_run
  use define
  use memory_mod, only : shared_host_rank, shared_host_size, shared_host_comm, &
                         create_shared_mem, destroy_shared_mem, get_window, &
                         sum_over_ranks
  use photoelectron_degradation_mod, only : N_EDEP, E_SECONDARY, &
                                            photoelectron_energy_weights
  implicit none
  private

  public :: ion_score_setup, ion_score_resize, ion_score_reset
  public :: score_ion_path, ion_score_reduce
  public :: ion_score_apply_hhe
  public :: ion_score_apply_dust_heat
  public :: band_path_luminosity
  public :: ion_score_apply_sec
  public :: ion_score_energy_packet, ion_score_energy_totals
  public :: ion_score_scattered_escape, ion_score_lost_weight
  public :: ion_score_report_energy_closure, ion_score_report_metal_rates
  public :: ion_score_report_ledger, ion_score_ledger_leaf
  public :: LEDGER_NSRC
  !--- interactions of the sampled histories (raytrace_amr_mod) and the
  !--- quantities of par%diffuse_launch = 'absorption_reemission'
  public :: ion_score_forced_interaction, ion_score_interaction
  public :: ion_score_interaction_absorption, ion_score_reemission
  public :: ion_score_roulette
  public :: gas_absorbed_power, previous_photoionizations
  public :: ion_score_continuation_factor, ion_score_gain_history_weight
  public :: ion_score_reemission_totals, ion_score_report_reemission
  public :: ion_score_absorbed_power_leaf, ion_score_reemission_leaf, &
            ion_score_interaction_energy_totals, ion_score_ledger_channels
  public :: NDIFFUSE_CHANNEL, LEDGER_NABS

  !--- how a walk of raytrace_amr_mod ended, for the energy totals below.
  !--- FLIGHT_STOPPED_IN_GRID: a sampled flight reached its optical depth
  !--- inside the grid (the packet interacts there).  FLIGHT_ESCAPED: the
  !--- flight left through an open face (any face of a closed box, a z face
  !--- of the xy-periodic slab); for the edge walk, the direct transmission
  !--- e^-tau_edge leaves the grid (it is below double precision once the
  !--- walk stops at its opaque cutoff).  FLIGHT_LOST: the weight left the
  !--- calculation without leaving the nebula -- MAX_WALK_STEPS reached, or
  !--- no leaf found on re-entering the xy-periodic slab, or a start outside
  !--- the grid.  Russian roulette is not a loss: the survivors' 1/PSURV
  !--- keeps the expectation.
  integer, parameter, public :: FLIGHT_STOPPED_IN_GRID = 0
  integer, parameter, public :: FLIGHT_ESCAPED = 1
  integer, parameter, public :: FLIGHT_LOST = 2

  integer, parameter :: I_GHI = 1, I_GHEI = 2, I_GHEII = 3
  integer, parameter :: I_HHI = 4, I_HHEI = 5, I_HHEII = 6
  integer, parameter :: NHHE = 6
  !--- secondary-ionization (par%use_sec_ion) path-estimator layout: for each
  !--- absorber (H I, He I, He II) the photoelectron energy deposited above
  !--- E_SECONDARY, sigma (1 - E_th/E), split onto the N_EDEP energy nodes of
  !--- photoelectron_degradation_mod by the x-independent weights of E0, so
  !--- the state-dependent partition is applied afterwards exactly.  Row
  !--- (a-1)*N_EDEP + k is absorber a, node k.
  integer, parameter :: NSEC = 3*N_EDEP
  integer, parameter :: ENERGY_ACC_WP = selected_real_kind(18,300)

  real(kind=wp), allocatable :: hhe_rate_path(:,:) ! (NHHE,nleaf), rank local then
                                                   ! reduced; Sum(path_lum * coeff(E))
                                                   ! for the three photoionization
                                                   ! rates and the three photoheating
                                                   ! integrals.
  real(kind=wp), allocatable :: band_path(:)       ! (nleaf), rank local then reduced;
                                                   ! Sum(path_lum) over every scored
                                                   ! segment of every band the gas
                                                   ! iteration transports (ionizing,
                                                   ! FUV), all energies: positive
                                                   ! exactly when a packet crossed
                                                   ! the leaf in the pass.
  real(kind=wp), allocatable :: dust_abs_path(:)   ! (nleaf), rank local then reduced;
                                                   ! Sum(path_lum * sigma_abs,dust(E))
                                                   ! per leaf, the EUV grain-heating
                                                   ! integral.
  real(kind=wp), allocatable :: sec_ion_path(:,:)  ! (NSEC,nleaf), rank local then
                                                   ! reduced; the photoelectron
                                                   ! energy deposited above 30 eV,
                                                   ! by absorber and energy node.
  integer :: score_nleaf = 0
  integer(kind=int64) :: direct_segments = 0_int64
  !--- segments of the sampled flights (photon%nevent > 0)
  integer(kind=int64) :: sampled_segments = 0_int64
  !--- energy totals of the transported packets [erg/s], rank local in the
  !--- extended kind (_acc), then reduced.
  !--- (a) The direct split: each packet's luminosity L at its edge
  !---     transmission T = e^-tau_edge, extinguished L(1-T) and transmitted
  !---     L T.  emitted = extinguished + transmitted holds packet by packet,
  !---     so it checks the accumulation and reduction alone
  !---     (ion_score_report_energy_closure aborts past 5e-13).
  !--- (b) What the direct split cannot see: the weight of the sampled
  !---     flights after an interaction (dust-scattered or, under
  !---     absorption_reemission, re-emitted) that leave the grid
  !---     (scattered_escaped), the weight lost
  !---     by such a flight (lost_scattered, which comes out of the
  !---     extinguished part) and the transmission of an edge walk that was
  !---     lost (lost_direct, which comes out of the transmitted part).
  !--- The physical totals are then
  !---     escaped  = transmitted - lost_direct + scattered_escaped
  !---     lost     = lost_direct + lost_scattered
  !---     absorbed = extinguished - scattered_escaped - lost_scattered
  !---              = emitted - escaped - lost,
  !--- formed from the extinguished part rather than from the emitted one so
  !--- that a run without scattering and without losses returns the direct
  !--- split itself, and an optically thin run keeps the precision of its
  !--- small absorbed power.  With scattering the absorbed total is exact
  !--- in expectation: the weight a_eff takes off at every interaction is
  !--- the absorbed share, and Russian roulette conserves the expectation.
  real(kind=wp) :: energy_emitted = 0.0_wp
  real(kind=wp) :: energy_extinguished = 0.0_wp
  real(kind=wp) :: energy_transmitted = 0.0_wp
  real(kind=wp) :: energy_scattered_escaped = 0.0_wp
  real(kind=wp) :: energy_lost_direct = 0.0_wp
  real(kind=wp) :: energy_lost_scattered = 0.0_wp
  !--- absorbed power by the path estimators the rates use [erg/s],
  !--- Sum path_lum kappa_abs over the scored segments (absorption_ledger
  !--- only); in expectation equal to the absorbed total above.
  real(kind=wp) :: energy_path_absorbed = 0.0_wp
  integer(kind=int64) :: lost_flights = 0_int64
  real(kind=ENERGY_ACC_WP) :: energy_emitted_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_extinguished_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_transmitted_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_scattered_escaped_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_lost_direct_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_lost_scattered_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_path_absorbed_acc = 0.0_ENERGY_ACC_WP
  !--- absorption ledger (par%absorption_ledger) [photons/s]: by source
  !--- (photon%source_channel 0..LEDGER_NSRC-1) the photons emitted, those
  !--- escaping (the direct transmission e^-tau_edge of every packet whose
  !--- edge walk was not lost, plus the weight of every scattered flight
  !--- that leaves the grid), and the
  !--- photons absorbed by each absorber, path_lum kappa_a / (E erg) on the
  !--- scored segments
  integer, parameter :: LEDGER_NSRC = 10, LEDGER_NABS = 5
  character(len=6), parameter :: LEDGER_ABS_NAME(LEDGER_NABS) = &
     ['H0    ', 'He0   ', 'He+   ', 'metals', 'dust  ']
  real(kind=ENERGY_ACC_WP) :: ledger_emitted(0:LEDGER_NSRC-1) = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: ledger_escaped(0:LEDGER_NSRC-1) = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: ledger_absorbed(LEDGER_NABS,0:LEDGER_NSRC-1) = 0.0_ENERGY_ACC_WP
  !--- the same for H0, He0, He+ in every leaf [photons/s]:
  !--- (absorber, source, leaf), rank local then reduced, written to the
  !--- rates file as Ledger_H0 / Ledger_He0 / Ledger_He+ (nleaf, source)
  integer, parameter :: LEDGER_NLEAF_ABS = 3
  real(kind=wp), allocatable :: ledger_leaf(:,:,:)
  !--- the ledger's photons lost (FLIGHT_LOST, as lost_direct and
  !--- lost_scattered above) and moved by Russian roulette (killed - gained),
  !--- by the channel the packet carried then [photons/s].  With them each
  !--- channel closes in expectation, emitted = absorbed + escaped + lost,
  !--- the roulette term being zero in expectation.
  real(kind=ENERGY_ACC_WP) :: ledger_lost(0:LEDGER_NSRC-1) = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: ledger_roulette(0:LEDGER_NSRC-1) = 0.0_ENERGY_ACC_WP

  !--- (c) The interactions of the sampled histories (transport_ion_packet;
  !---     ionizing band only; a pass samples them with dust scattering or
  !---     under par%diffuse_launch = 'absorption_reemission').  The forced
  !---     first interaction takes the weight L (1 - e^-tau_edge) of the
  !---     launched packet (energy_forced).  At each interaction the weight
  !---     that does not continue, incoming - continued, is deposited
  !---     (energy_deposited, D: the dust absorption and the gas absorption
  !---     that is not re-emitted); Russian roulette moves
  !---     Q = killed - gained (energy_roulette, zero in expectation); the
  !---     continued flights end as scattered_escaped (Cesc), lost_scattered
  !---     (Clost) or at a next interaction.  Packet by packet, to rounding,
  !---         energy_forced = D + Cesc + Clost + Q,
  !---     and energy_forced = energy_extinguished (every tau_edge > 0
  !---     forces an interaction, raytrace_amr_mod::transport_ion_packet).
  !--- (d) absorption_reemission only.  At an interaction with incoming
  !---     luminosity W in leaf i, with kap = kap_gas + kap_abs,dust +
  !---     kap_sca,dust at the packet energy:
  !---         G_event  = W kap_gas/kap          gross gas absorption
  !---         Dd_event = W kap_abs,dust/kap     dust absorption
  !---         R_exp    = W Sum_s y_s(E) k_s/kap   expected re-emission
  !---                    (y_s, k_s of raytrace_amr_mod, plan Section 12)
  !---     (energy_event_gas, energy_event_dust, energy_reemission_expected),
  !---     and Re, the continued luminosity of every selected re-emission
  !---     branch, by the channel drawn (energy_reemitted).  The physical
  !---     totals (ion_score_energy_totals, decision D6 of the plan) are then
  !---         L_EMIT = sources + explicit remainder + Re
  !---         L_ABS  = L_EMIT - L_ESC - L_LOST = D + Q + Re   (gross)
  !---     and L_NETABS = L_ABS - Re = D + Q, the power removed from the
  !---     transported field.  Sampled quantities (G_event, Re) and their
  !---     conditional means (R_exp) carry separate names.
  integer, parameter :: NDIFFUSE_CHANNEL = 8
  real(kind=wp) :: energy_forced = 0.0_wp
  real(kind=wp) :: energy_deposited = 0.0_wp
  real(kind=wp) :: energy_roulette = 0.0_wp
  real(kind=wp) :: energy_event_gas = 0.0_wp
  real(kind=wp) :: energy_event_dust = 0.0_wp
  real(kind=wp) :: energy_reemission_expected = 0.0_wp
  real(kind=wp) :: energy_reemitted(NDIFFUSE_CHANNEL) = 0.0_wp
  real(kind=ENERGY_ACC_WP) :: energy_forced_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_deposited_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_roulette_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_event_gas_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_event_dust_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_reemission_expected_acc = 0.0_ENERGY_ACC_WP
  real(kind=ENERGY_ACC_WP) :: energy_reemitted_acc(NDIFFUSE_CHANNEL) = 0.0_ENERGY_ACC_WP
  integer(kind=int64) :: interaction_count = 0_int64
  !--- the largest continuation factor t of an interaction in the pass
  !--- (absorption_reemission; t may exceed 1 where the in-band energy of a
  !--- recombination exceeds the absorbed photon's), reduced by MAX
  real(kind=wp) :: max_continuation = 0.0_wp
  !--- the largest weight a history reached after an interaction with t > 1
  !--- (after the interaction or a roulette survival), reduced by MAX
  real(kind=wp) :: max_gain_weight = 0.0_wp
  integer(kind=int64) :: reemission_count = 0_int64
  !--- the same luminosities by the root source of the history
  !--- (photon%root_channel, never relabeled) [erg/s]: launched, escaped
  !--- (direct transmission and continued flights), lost, deposited (D),
  !--- roulette (Q) and re-emitted (Re).  When every launched packet's
  !--- extinguished weight enters the sampled history (absorption_reemission)
  !--- each root closes to rounding: launched = escaped + lost + D + Q.
  integer, parameter :: ROOT_LAUNCHED = 1, ROOT_ESCAPED = 2, ROOT_LOST = 3, &
                        ROOT_DEPOSITED = 4, ROOT_ROULETTE = 5, ROOT_REEMITTED = 6
  integer, parameter :: NROOT_TERM = 6
  real(kind=wp) :: root_energy(NROOT_TERM, 0:LEDGER_NSRC-1) = 0.0_wp
  real(kind=ENERGY_ACC_WP) :: root_energy_acc(NROOT_TERM, 0:LEDGER_NSRC-1) = 0.0_ENERGY_ACC_WP

  !--- measured gas absorption of every leaf [erg/s] (absorption_reemission
  !--- only, see the header): gas_absorbed_now is this pass's rank-local
  !--- path estimate, reduced in place by ion_score_reduce; gas_absorbed_power
  !--- (node-shared, the node's ranks writing disjoint leaves) holds the
  !--- reduced values until the next pass's reduction.  Neither exists in
  !--- the default mode.
  logical :: reemission_on = .false.
  !--- the H/He and metal path sums hold the reduced field of a pass on the
  !--- present grid: set by ion_score_reduce, cleared by ion_score_setup/
  !--- resize (a new grid, a restart) and ion_score_reset (the sums are being
  !--- accumulated again).  previous_photoionizations returns valid = .false. without it.
  logical :: absorbed_field_valid = .false.
  real(kind=wp), allocatable :: gas_absorbed_now(:)
  !--- absorption_reemission leaf arrays for the output (rank local, then
  !--- reduced by ion_score_reduce) [erg/s]: the sampled re-emission by
  !--- absorbing leaf and new channel (the continued luminosity before
  !--- roulette; its sum is Re); the gross gas absorption of the
  !--- interactions, G_event, by leaf.
  real(kind=wp), allocatable :: reemitted_leaf(:,:)     ! (NDIFFUSE_CHANNEL, nleaf)
  real(kind=wp), allocatable :: event_gas_leaf(:)
  real(kind=wp), pointer, contiguous, protected :: gas_absorbed_power(:) => null()

contains

  subroutine ion_score_setup(nleaf)
    implicit none
    integer, intent(in) :: nleaf

    if (allocated(hhe_rate_path)) deallocate(hhe_rate_path)
    if (allocated(band_path)) deallocate(band_path)
    if (allocated(dust_abs_path)) deallocate(dust_abs_path)
    if (allocated(sec_ion_path)) deallocate(sec_ion_path)
    score_nleaf = nleaf
    allocate(hhe_rate_path(NHHE,nleaf))
    hhe_rate_path = 0.0_wp
    allocate(band_path(nleaf))
    band_path = 0.0_wp

    !--- EUV grain-heating integral: allocated only when dust is on.  Left
    !--- unallocated otherwise, so score_ion_path, the reduction, and the reset
    !--- are no-ops on the dust-off path.
    if (par%ion_add_dust) then
       allocate(dust_abs_path(nleaf))
       dust_abs_path = 0.0_wp
    end if

    !--- secondary-ionization integrals: allocated only when use_sec_ion is on.
    if (par%use_sec_ion) then
       allocate(sec_ion_path(NSEC,nleaf))
       sec_ion_path = 0.0_wp
    end if

    if (par%use_metals) then
       block
         use species_mod, only : species_rate_path_setup
         call species_rate_path_setup(nleaf)
       end block
    end if
    if (allocated(ledger_leaf)) deallocate(ledger_leaf)
    if (par%absorption_ledger) then
       allocate(ledger_leaf(LEDGER_NLEAF_ABS, 0:LEDGER_NSRC-1, nleaf))
       ledger_leaf = 0.0_wp
    end if
    !--- He I 2^3S photoionization rate Phi_3 and its heating (the output
    !--- diagnostic par%hei_metastable and the coupling
    !--- par%hei_metastable_coupling)
    if (par%hei_metastable .or. par%hei_metastable_coupling) then
       block
         use hei_metastable_mod, only : hei_metastable_path_setup
         call hei_metastable_path_setup(nleaf)
       end block
    end if
    !--- measured gas absorption (absorption_reemission): created anew on a
    !--- new grid (create_shared_mem zeroes the window).
    !--- Collective on the node, as every caller (main, the re-refinement)
    !--- runs on all ranks.
    reemission_on = trim(par%diffuse_launch) == 'absorption_reemission'
    absorbed_field_valid = .false.
    if (allocated(gas_absorbed_now)) deallocate(gas_absorbed_now)
    if (allocated(reemitted_leaf)) deallocate(reemitted_leaf)
    if (allocated(event_gas_leaf)) deallocate(event_gas_leaf)
    if (associated(gas_absorbed_power)) call destroy_shared_mem(gas_absorbed_power)
    if (reemission_on) then
       allocate(gas_absorbed_now(nleaf), &
                reemitted_leaf(NDIFFUSE_CHANNEL,nleaf), event_gas_leaf(nleaf))
       gas_absorbed_now = 0.0_wp
       reemitted_leaf = 0.0_wp
       event_gas_leaf = 0.0_wp
       call create_shared_mem(gas_absorbed_power, [nleaf])
    end if
  end subroutine ion_score_setup

  subroutine ion_score_resize(nleaf)
    integer, intent(in) :: nleaf
    call ion_score_setup(nleaf)
  end subroutine ion_score_resize

  subroutine ion_score_reset()
    if (allocated(hhe_rate_path)) hhe_rate_path = 0.0_wp
    if (allocated(band_path)) band_path = 0.0_wp
    if (allocated(dust_abs_path)) dust_abs_path = 0.0_wp
    if (allocated(sec_ion_path)) sec_ion_path = 0.0_wp
    direct_segments = 0_int64
    sampled_segments = 0_int64
    energy_emitted = 0.0_wp
    energy_extinguished = 0.0_wp
    energy_transmitted = 0.0_wp
    energy_scattered_escaped = 0.0_wp
    energy_lost_direct = 0.0_wp
    energy_lost_scattered = 0.0_wp
    energy_path_absorbed = 0.0_wp
    lost_flights = 0_int64
    energy_emitted_acc = 0.0_ENERGY_ACC_WP
    energy_extinguished_acc = 0.0_ENERGY_ACC_WP
    energy_transmitted_acc = 0.0_ENERGY_ACC_WP
    energy_scattered_escaped_acc = 0.0_ENERGY_ACC_WP
    energy_lost_direct_acc = 0.0_ENERGY_ACC_WP
    energy_lost_scattered_acc = 0.0_ENERGY_ACC_WP
    energy_path_absorbed_acc = 0.0_ENERGY_ACC_WP
    ledger_emitted = 0.0_ENERGY_ACC_WP
    ledger_escaped = 0.0_ENERGY_ACC_WP
    ledger_absorbed = 0.0_ENERGY_ACC_WP
    ledger_lost = 0.0_ENERGY_ACC_WP
    ledger_roulette = 0.0_ENERGY_ACC_WP
    if (allocated(ledger_leaf)) ledger_leaf = 0.0_wp
    energy_forced = 0.0_wp
    energy_deposited = 0.0_wp
    energy_roulette = 0.0_wp
    energy_event_gas = 0.0_wp
    energy_event_dust = 0.0_wp
    energy_reemission_expected = 0.0_wp
    energy_reemitted = 0.0_wp
    energy_forced_acc = 0.0_ENERGY_ACC_WP
    energy_deposited_acc = 0.0_ENERGY_ACC_WP
    energy_roulette_acc = 0.0_ENERGY_ACC_WP
    energy_event_gas_acc = 0.0_ENERGY_ACC_WP
    energy_event_dust_acc = 0.0_ENERGY_ACC_WP
    energy_reemission_expected_acc = 0.0_ENERGY_ACC_WP
    energy_reemitted_acc = 0.0_ENERGY_ACC_WP
    interaction_count = 0_int64
    reemission_count = 0_int64
    max_continuation = 0.0_wp
    max_gain_weight = 0.0_wp
    root_energy = 0.0_wp
    root_energy_acc = 0.0_ENERGY_ACC_WP
    !--- this pass's estimate only; the published value is kept until
    !--- ion_score_reduce replaces it
    if (allocated(gas_absorbed_now)) gas_absorbed_now = 0.0_wp
    absorbed_field_valid = .false.
    if (allocated(reemitted_leaf)) reemitted_leaf = 0.0_wp
    if (allocated(event_gas_leaf)) event_gas_leaf = 0.0_wp
    if (par%use_metals) then
       block
         use species_mod, only : species_rate_path_reset
         call species_rate_path_reset()
       end block
    end if
    if (par%hei_metastable .or. par%hei_metastable_coupling) then
       block
         use hei_metastable_mod, only : hei_metastable_path_reset
         call hei_metastable_path_reset()
       end block
    end if
  end subroutine ion_score_reset

  subroutine score_ion_path(photon, il, path_lum)
    use jtally_mod,       only : jt_ion, jt_ion_on
    use dust_ir_band_mod, only : transported_band, BAND_DUST_IR, &
                                 dust_ir_absorb_path
    implicit none
    type(photon_type), intent(in) :: photon
    integer,           intent(in) :: il
    real(kind=wp),     intent(in) :: path_lum
    integer :: inu
    real(kind=wp) :: energy, coeff(NHHE)

    if (path_lum == 0.0_wp) return
    !--- Dust thermal re-emission packets sit below every gas photoionization
    !--- threshold, so they carry no photoionization, no photoheating, no
    !--- secondary ionization, and their energies lie outside
    !--- the ionizing-band bin grid jt_ion is defined on.  Their only
    !--- interaction is grain absorption, which is scored into the infrared
    !--- band's own estimator; every ionizing-band estimator is left untouched.
    if (transported_band == BAND_DUST_IR) then
       call dust_ir_absorb_path(il, path_lum*photon%ionphys%dust_abs)
       return
    end if
    inu = photon%inu
    if (jt_ion_on) jt_ion(inu,il) = jt_ion(inu,il) + path_lum
    band_path(il) = band_path(il) + path_lum
    !--- EUV grain-heating integral: accumulate the packet's dust absorption
    !--- cross section, evaluated at its exact energy.
    if (allocated(dust_abs_path)) &
       dust_abs_path(il) = dust_abs_path(il) + path_lum*photon%ionphys%dust_abs
    !--- photoelectron energy above E_SECONDARY, by absorber and energy node,
    !--- at the packet's own energy.
    if (allocated(sec_ion_path)) then
       block
         real(kind=wp) :: en, sig(3), eth(3), w(N_EDEP)
         integer :: ia
         en  = photon%ionphys%energy_eV
         sig = [photon%ionphys%sigma_HI, photon%ionphys%sigma_HeI, &
                photon%ionphys%sigma_HeII]
         eth = [eth_HI, eth_HeI, eth_HeII]
         do ia = 1, 3
            if (en - eth(ia) <= E_SECONDARY) cycle
            call photoelectron_energy_weights(en - eth(ia), w)
            sec_ion_path((ia-1)*N_EDEP+1:ia*N_EDEP, il) = &
               sec_ion_path((ia-1)*N_EDEP+1:ia*N_EDEP, il) &
               + path_lum*sig(ia)*(1.0_wp - eth(ia)/en)*w
         end do
       end block
    end if
    if (photon%nevent > 0) then
       sampled_segments = sampled_segments + 1_int64
    else
       direct_segments = direct_segments + 1_int64
    end if
    energy = photon%ionphys%energy_eV
    coeff(I_GHI)   = photon%ionphys%sigma_HI/(energy*ev2erg)
    coeff(I_GHEI)  = photon%ionphys%sigma_HeI/(energy*ev2erg)
    coeff(I_GHEII) = photon%ionphys%sigma_HeII/(energy*ev2erg)
    coeff(I_HHI)   = photon%ionphys%sigma_HI * &
                     max(1.0_wp-eth_HI/energy,0.0_wp)
    coeff(I_HHEI)  = photon%ionphys%sigma_HeI * &
                     max(1.0_wp-eth_HeI/energy,0.0_wp)
    coeff(I_HHEII) = photon%ionphys%sigma_HeII * &
                     max(1.0_wp-eth_HeII/energy,0.0_wp)
    hhe_rate_path(:,il) = hhe_rate_path(:,il) + path_lum*coeff
    if (par%use_metals) then
       block
         use species_mod, only : species_rate_path_score
         call species_rate_path_score(photon%ionphys%energy_eV, &
              photon%ionphys%sigma_metal, photon%ionphys%heat_metal, &
              photon%ionphys%n_metal, il, path_lum)
       end block
    end if
    if (par%hei_metastable .or. par%hei_metastable_coupling) then
       block
         use hei_metastable_mod, only : hei_metastable_path_score
         call hei_metastable_path_score(il, photon%ionphys%energy_eV, path_lum)
       end block
    end if
    if (par%absorption_ledger) call ledger_score(photon, il, path_lum)
    !--- gas-absorbed power A of the leaf (absorption_reemission)
    if (reemission_on) then
       block
         use ion_packet_mod, only : ion_packet_gas_opacity
         gas_absorbed_now(il) = gas_absorbed_now(il) &
                                + path_lum*ion_packet_gas_opacity(photon, il)
       end block
    end if
  end subroutine score_ion_path

  !--- photons absorbed on one scored segment, by absorber, with the
  !--- opacity terms of ion_packet_opacity (dust: absorption only)
  subroutine ledger_score(photon, il, path_lum)
    use octree_mod,    only : amr_grid
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII
    type(photon_type), intent(in) :: photon
    integer,           intent(in) :: il
    real(kind=wp),     intent(in) :: path_lum
    real(kind=wp) :: kap(LEDGER_NABS), nH
    integer :: is
    nH = gas_nH(il)
    kap = 0.0_wp
    kap(1) = nH*gas_xHI(il)*photon%ionphys%sigma_HI*par%distance2cm
    kap(2) = nH*par%abund_He*gas_xHeI(il)*photon%ionphys%sigma_HeI*par%distance2cm
    kap(3) = nH*par%abund_He*gas_xHeII(il)*photon%ionphys%sigma_HeII*par%distance2cm
    if (par%use_metals .and. par%ion_metal_abs) then
       block
         use species_mod, only : species_cached_opacity
         kap(4) = nH*species_cached_opacity(photon%ionphys%sigma_metal, &
                  photon%ionphys%n_metal, il)*par%distance2cm
       end block
    end if
    if (par%ion_add_dust) kap(5) = amr_grid%rhokap(il)*photon%ionphys%dust_abs
    is = min(max(photon%source_channel, 0), LEDGER_NSRC-1)
    kap = path_lum*kap/(photon%ionphys%energy_eV*ev2erg)
    ledger_absorbed(:,is) = ledger_absorbed(:,is) + real(kap, ENERGY_ACC_WP)
    ledger_leaf(:,is,il) = ledger_leaf(:,is,il) + kap(1:LEDGER_NLEAF_ABS)
    !--- the same absorptions in energy [erg/s], for the comparison with the
    !--- physical absorbed total (ion_score_report_energy_closure)
    energy_path_absorbed_acc = energy_path_absorbed_acc &
       + real(sum(kap), ENERGY_ACC_WP)*real(photon%ionphys%energy_eV*ev2erg, ENERGY_ACC_WP)
  end subroutine ledger_score

  !--- photons/s absorbed by absorber ia (1 H0, 2 He0, 3 He+) in every leaf
  !--- from source is, after ion_score_reduce
  subroutine ion_score_ledger_leaf(ia, is, out)
    integer,       intent(in)  :: ia, is
    real(kind=wp), intent(out) :: out(:)
    out = ledger_leaf(ia, is, 1:size(out))
  end subroutine ion_score_ledger_leaf

  !--- one launched packet at its edge walk: the direct split of its
  !--- luminosity by the edge transmission e^-tau_edge.  edge_flight is how
  !--- the edge walk ended; FLIGHT_LOST (MAX_WALK_STEPS, or no leaf on
  !--- re-entering the xy-periodic slab) books the transmitted part as lost
  !--- instead of escaped.  That cannot undo the path estimates the walk has
  !--- already scored up to where it stopped, nor the forced interaction
  !--- drawn on its truncated tau_edge; it only keeps the truncated
  !--- transmission out of the escape and makes the event visible (lost
  !--- weight and count in the energy report, slab_walk_kills).
  subroutine ion_score_energy_packet(photon, tau_edge, edge_flight)
    use dust_ir_band_mod, only : transported_band, BAND_DUST_IR, &
                                 dust_ir_escape_packet
    type(photon_type), intent(in) :: photon
    real(kind=wp), intent(in) :: tau_edge
    integer, intent(in) :: edge_flight
    real(kind=wp) :: luminosity, transmission
    logical :: edge_lost
    edge_lost = edge_flight == FLIGHT_LOST
    !--- an infrared packet belongs to the dust thermal band's own ledger, not
    !--- to the ionizing band's.
    if (transported_band == BAND_DUST_IR) then
       call dust_ir_escape_packet(photon, tau_edge, edge_lost)
       return
    end if
    luminosity = photon%Lpacket*photon%wgt
    transmission = exp(-max(tau_edge,0.0_wp))
    energy_emitted_acc = energy_emitted_acc + real(luminosity, ENERGY_ACC_WP)
    !--- L (1 - T) at full relative accuracy for every tau_edge
    !--- (define::extinguished_fraction; L - L T vanished in thin media)
    energy_extinguished_acc = energy_extinguished_acc + real( &
       luminosity*extinguished_fraction(max(tau_edge,0.0_wp)), ENERGY_ACC_WP)
    energy_transmitted_acc = energy_transmitted_acc + real( &
       luminosity*transmission, ENERGY_ACC_WP)
    if (edge_lost) then
       energy_lost_direct_acc = energy_lost_direct_acc + real( &
          luminosity*transmission, ENERGY_ACC_WP)
       lost_flights = lost_flights + 1_int64
    end if
    if (par%absorption_ledger) then
       block
         integer :: is
         real(kind=wp) :: nphot
         is = min(max(photon%source_channel, 0), LEDGER_NSRC-1)
         nphot = luminosity/(photon%energy_eV*ev2erg)
         ledger_emitted(is) = ledger_emitted(is) + real(nphot, ENERGY_ACC_WP)
         if (.not. edge_lost) ledger_escaped(is) = ledger_escaped(is) &
            + real(nphot*transmission, ENERGY_ACC_WP)
         if (edge_lost) ledger_lost(is) = ledger_lost(is) &
            + real(nphot*transmission, ENERGY_ACC_WP)
       end block
    end if
    block
      integer :: ir
      ir = root_index(photon)
      root_energy_acc(ROOT_LAUNCHED,ir) = root_energy_acc(ROOT_LAUNCHED,ir) &
         + real(luminosity, ENERGY_ACC_WP)
      if (edge_lost) then
         root_energy_acc(ROOT_LOST,ir) = root_energy_acc(ROOT_LOST,ir) &
            + real(luminosity*transmission, ENERGY_ACC_WP)
      else
         root_energy_acc(ROOT_ESCAPED,ir) = root_energy_acc(ROOT_ESCAPED,ir) &
            + real(luminosity*transmission, ENERGY_ACC_WP)
      end if
    end block
  end subroutine ion_score_energy_packet

  !--- ledger row of the channel a packet carries now, and of its root
  integer function channel_index(photon) result(is)
    type(photon_type), intent(in) :: photon
    is = min(max(photon%source_channel, 0), LEDGER_NSRC-1)
  end function channel_index

  integer function root_index(photon) result(ir)
    type(photon_type), intent(in) :: photon
    ir = min(max(photon%root_channel, 0), LEDGER_NSRC-1)
  end function root_index

  !=========================================================================
  ! Interactions of a sampled history (raytrace_amr_mod::transport_ion_packet),
  ! ionizing band only (an infrared packet returns at once): the definitions
  ! are at the module variables, (c) and (d).
  !=========================================================================
  !--- the forced first interaction: photon%wgt already carries the factor
  !--- 1 - e^-tau_edge, so L wgt is the weight entering the history
  subroutine ion_score_forced_interaction(photon)
    use dust_ir_band_mod, only : transported_band, BAND_DUST_IR
    type(photon_type), intent(in) :: photon
    if (transported_band == BAND_DUST_IR) return
    energy_forced_acc = energy_forced_acc &
       + real(photon%Lpacket*photon%wgt, ENERGY_ACC_WP)
  end subroutine ion_score_forced_interaction

  !--- one interaction: wgt_in the incoming weight, photon%wgt the continued
  !--- one (zero when the history ends there); the difference is deposited
  subroutine ion_score_interaction(photon, wgt_in)
    use dust_ir_band_mod, only : transported_band, BAND_DUST_IR
    type(photon_type), intent(in) :: photon
    real(kind=wp),     intent(in) :: wgt_in
    real(kind=wp) :: deposit
    integer :: ir
    if (transported_band == BAND_DUST_IR) return
    deposit = photon%Lpacket*wgt_in - photon%Lpacket*photon%wgt
    energy_deposited_acc = energy_deposited_acc + real(deposit, ENERGY_ACC_WP)
    ir = root_index(photon)
    root_energy_acc(ROOT_DEPOSITED,ir) = root_energy_acc(ROOT_DEPOSITED,ir) &
       + real(deposit, ENERGY_ACC_WP)
    interaction_count = interaction_count + 1_int64
  end subroutine ion_score_interaction

  !--- absorption_reemission: the gross gas and dust absorption and the
  !--- expected re-emission of one interaction with incoming weight wgt_in,
  !--- from the shares kap_gas/kap, kap_abs,dust/kap and Sum_s y_s k_s/kap
  subroutine ion_score_interaction_absorption(photon, wgt_in, gas_share, &
                                              dust_share, reemit_share)
    use dust_ir_band_mod, only : transported_band, BAND_DUST_IR
    type(photon_type), intent(in) :: photon
    real(kind=wp),     intent(in) :: wgt_in, gas_share, dust_share, reemit_share
    real(kind=wp) :: w
    if (transported_band == BAND_DUST_IR) return
    w = photon%Lpacket*wgt_in
    energy_event_gas_acc  = energy_event_gas_acc  + real(w*gas_share, ENERGY_ACC_WP)
    energy_event_dust_acc = energy_event_dust_acc + real(w*dust_share, ENERGY_ACC_WP)
    energy_reemission_expected_acc = energy_reemission_expected_acc &
       + real(w*reemit_share, ENERGY_ACC_WP)
    if (allocated(event_gas_leaf) .and. photon%icell_amr > 0) &
       event_gas_leaf(photon%icell_amr) = event_gas_leaf(photon%icell_amr) + w*gas_share
  end subroutine ion_score_interaction_absorption

  !--- a selected re-emission branch: photon carries its new channel and
  !--- energy and the continued weight (before Russian roulette).  Booked as
  !--- re-emitted luminosity Re of that channel and, with
  !--- par%absorption_ledger, as photons emitted in it, L wgt / (E_new erg).
  subroutine ion_score_reemission(photon)
    use dust_ir_band_mod, only : transported_band, BAND_DUST_IR
    type(photon_type), intent(in) :: photon
    real(kind=wp) :: luminosity
    integer :: is, ir
    if (transported_band == BAND_DUST_IR) return
    luminosity = photon%Lpacket*photon%wgt
    is = channel_index(photon)
    if (is >= 1 .and. is <= NDIFFUSE_CHANNEL) then
       energy_reemitted_acc(is) = energy_reemitted_acc(is) &
                                  + real(luminosity, ENERGY_ACC_WP)
       !--- by absorbing leaf (the packet has not moved since the absorption)
       if (allocated(reemitted_leaf) .and. photon%icell_amr > 0) &
          reemitted_leaf(is,photon%icell_amr) = reemitted_leaf(is,photon%icell_amr) &
                                                + luminosity
    end if
    ir = root_index(photon)
    root_energy_acc(ROOT_REEMITTED,ir) = root_energy_acc(ROOT_REEMITTED,ir) &
       + real(luminosity, ENERGY_ACC_WP)
    reemission_count = reemission_count + 1_int64
    if (par%absorption_ledger) ledger_emitted(is) = ledger_emitted(is) &
       + real(luminosity/(photon%energy_eV*ev2erg), ENERGY_ACC_WP)
  end subroutine ion_score_reemission

  !--- the continuation factor t of one interaction (its maximum is reported)
  subroutine ion_score_continuation_factor(t)
    real(kind=wp), intent(in) :: t
    max_continuation = max(max_continuation, t)
  end subroutine ion_score_continuation_factor

  !--- a weight reached by a history that has had an interaction with t > 1
  subroutine ion_score_gain_history_weight(w)
    real(kind=wp), intent(in) :: w
    max_gain_weight = max(max_gain_weight, w)
  end subroutine ion_score_gain_history_weight

  !--- Russian roulette: wgt_in the weight before it; a survivor carries
  !--- the raised weight in photon%wgt, a killed packet leaves photon%wgt
  !--- as it was.  Q = killed - gained, in energy and, with
  !--- par%absorption_ledger, in photons of the current channel.
  subroutine ion_score_roulette(photon, wgt_in, survived)
    use dust_ir_band_mod, only : transported_band, BAND_DUST_IR
    type(photon_type), intent(in) :: photon
    real(kind=wp),     intent(in) :: wgt_in
    logical,           intent(in) :: survived
    real(kind=wp) :: q
    integer :: ir, is
    if (transported_band == BAND_DUST_IR) return
    if (survived) then
       q = -(photon%Lpacket*photon%wgt - photon%Lpacket*wgt_in)
    else
       q = photon%Lpacket*wgt_in
    end if
    energy_roulette_acc = energy_roulette_acc + real(q, ENERGY_ACC_WP)
    ir = root_index(photon)
    root_energy_acc(ROOT_ROULETTE,ir) = root_energy_acc(ROOT_ROULETTE,ir) &
       + real(q, ENERGY_ACC_WP)
    if (par%absorption_ledger) then
       is = channel_index(photon)
       ledger_roulette(is) = ledger_roulette(is) &
          + real(q/(photon%energy_eV*ev2erg), ENERGY_ACC_WP)
    end if
  end subroutine ion_score_roulette

  !--- a sampled flight after an interaction (dust scattering or
  !--- re-emission) left the grid through an open face
  !--- (raytrace_amr_mod::transport_ion_packet): its remaining weight
  !--- escapes.  The one entry point for both bands: an infrared packet goes
  !--- to the dust thermal band's escape and escaping spectrum; an ionizing
  !--- one to the escaped energy total always and, with
  !--- par%absorption_ledger, to its source's ledger escape [photons/s].
  subroutine ion_score_scattered_escape(photon)
    use dust_ir_band_mod, only : transported_band, BAND_DUST_IR, &
                                 dust_ir_scattered_escape
    type(photon_type), intent(in) :: photon
    real(kind=wp) :: luminosity
    integer :: is
    if (transported_band == BAND_DUST_IR) then
       call dust_ir_scattered_escape(photon)
       return
    end if
    luminosity = photon%wgt*photon%Lpacket
    energy_scattered_escaped_acc = energy_scattered_escaped_acc &
       + real(luminosity, ENERGY_ACC_WP)
    root_energy_acc(ROOT_ESCAPED,root_index(photon)) = &
       root_energy_acc(ROOT_ESCAPED,root_index(photon)) + real(luminosity, ENERGY_ACC_WP)
    if (.not. par%absorption_ledger) return
    is = min(max(photon%source_channel, 0), LEDGER_NSRC-1)
    ledger_escaped(is) = ledger_escaped(is) + real( &
       luminosity/(photon%energy_eV*ev2erg), ENERGY_ACC_WP)
  end subroutine ion_score_scattered_escape

  !--- the remaining weight of a scattered flight left the calculation
  !--- without leaving the nebula (FLIGHT_LOST, or a flight the scattering
  !--- loop cannot continue): booked as lost in the band it belongs to.
  !--- Normally zero; nonzero values are reported, not hidden.
  subroutine ion_score_lost_weight(photon)
    use dust_ir_band_mod, only : transported_band, BAND_DUST_IR, &
                                 dust_ir_lost_weight
    type(photon_type), intent(in) :: photon
    if (transported_band == BAND_DUST_IR) then
       call dust_ir_lost_weight(photon%wgt*photon%Lpacket)
       return
    end if
    energy_lost_scattered_acc = energy_lost_scattered_acc &
       + real(photon%wgt*photon%Lpacket, ENERGY_ACC_WP)
    lost_flights = lost_flights + 1_int64
    root_energy_acc(ROOT_LOST,root_index(photon)) = &
       root_energy_acc(ROOT_LOST,root_index(photon)) &
       + real(photon%wgt*photon%Lpacket, ENERGY_ACC_WP)
    if (par%absorption_ledger) ledger_lost(channel_index(photon)) = &
       ledger_lost(channel_index(photon)) &
       + real(photon%wgt*photon%Lpacket/(photon%energy_eV*ev2erg), ENERGY_ACC_WP)
  end subroutine ion_score_lost_weight

  subroutine ion_score_reduce()
    use mpi
    implicit none
    integer :: ierr
    integer(kind=int64) :: counts(3)
    real(kind=wp) :: energy_totals(7)

    call sum_over_ranks(hhe_rate_path)
    call sum_over_ranks(band_path)
    if (allocated(dust_abs_path)) call sum_over_ranks(dust_abs_path)
    if (allocated(sec_ion_path))  call sum_over_ranks(sec_ion_path)
    counts = [direct_segments, sampled_segments, lost_flights]
    call MPI_ALLREDUCE(MPI_IN_PLACE, counts, 3, MPI_INTEGER8, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    direct_segments = counts(1)
    sampled_segments = counts(2)
    lost_flights = counts(3)
    energy_totals = real([energy_emitted_acc, energy_extinguished_acc, &
       energy_transmitted_acc, energy_scattered_escaped_acc, &
       energy_lost_direct_acc, energy_lost_scattered_acc, &
       energy_path_absorbed_acc], wp)
    call MPI_ALLREDUCE(MPI_IN_PLACE, energy_totals, 7, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    energy_emitted = energy_totals(1)
    energy_extinguished = energy_totals(2)
    energy_transmitted = energy_totals(3)
    energy_scattered_escaped = energy_totals(4)
    energy_lost_direct = energy_totals(5)
    energy_lost_scattered = energy_totals(6)
    energy_path_absorbed = energy_totals(7)
    call interaction_totals_reduce()
    if (reemission_on) call gas_absorbed_power_publish()
    absorbed_field_valid = reemission_on
    if (par%absorption_ledger) call ledger_reduce()
    if (par%use_metals) then
       block
         use species_mod, only : species_rate_path_reduce
         call species_rate_path_reduce()
       end block
    end if
    if (par%hei_metastable .or. par%hei_metastable_coupling) then
       block
         use hei_metastable_mod, only : hei_metastable_path_reduce
         call hei_metastable_path_reduce()
       end block
    end if
  end subroutine ion_score_reduce

  !--- the interaction totals (c), (d) and the root sums over ranks, in
  !--- double precision (accumulated in the extended kind)
  subroutine interaction_totals_reduce()
    use mpi
    integer, parameter :: NSCAL = 6, NBUF = NSCAL + NDIFFUSE_CHANNEL &
                                          + NROOT_TERM*LEDGER_NSRC
    real(kind=wp) :: buf(NBUF)
    integer(kind=int64) :: counts(2)
    integer :: ierr, k
    buf(1:NSCAL) = real([energy_forced_acc, energy_deposited_acc, &
       energy_roulette_acc, energy_event_gas_acc, energy_event_dust_acc, &
       energy_reemission_expected_acc], wp)
    k = NSCAL
    buf(k+1:k+NDIFFUSE_CHANNEL) = real(energy_reemitted_acc, wp)
    k = k + NDIFFUSE_CHANNEL
    buf(k+1:) = reshape(real(root_energy_acc, wp), [NROOT_TERM*LEDGER_NSRC])
    call MPI_ALLREDUCE(MPI_IN_PLACE, buf, NBUF, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    energy_forced = buf(1)
    energy_deposited = buf(2)
    energy_roulette = buf(3)
    energy_event_gas = buf(4)
    energy_event_dust = buf(5)
    energy_reemission_expected = buf(6)
    energy_reemitted = buf(NSCAL+1:NSCAL+NDIFFUSE_CHANNEL)
    root_energy = reshape(buf(k+1:), [NROOT_TERM, LEDGER_NSRC])
    counts = [interaction_count, reemission_count]
    call MPI_ALLREDUCE(MPI_IN_PLACE, counts, 2, MPI_INTEGER8, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
    interaction_count = counts(1)
    reemission_count = counts(2)
    call MPI_ALLREDUCE(MPI_IN_PLACE, max_continuation, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, max_gain_weight, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
    if (reemission_on) then
       call sum_over_ranks(reemitted_leaf(:,1:score_nleaf))
       call sum_over_ranks(event_gas_leaf(1:score_nleaf))
    end if
  end subroutine interaction_totals_reduce

  !--- A of this pass: summed in place over all ranks (sum_over_ranks), then
  !--- published node-shared as gas_absorbed_power, the node's ranks writing
  !--- disjoint leaves.  A diagnostic (GasAbsorbedPower); the split of the
  !--- emission is formed from the photoionization counts
  !--- (previous_photoionizations).
  subroutine gas_absorbed_power_publish()
    integer :: il
    call sum_over_ranks(gas_absorbed_now(1:score_nleaf))
    do il = shared_host_rank() + 1, score_nleaf, shared_host_size()
       gas_absorbed_power(il) = gas_absorbed_now(il)
    end do
    call gas_absorbed_power_sync()
  end subroutine gas_absorbed_power_publish

  subroutine gas_absorbed_power_sync()
    use mpi
    integer :: win, ierr
    call get_window(gas_absorbed_power, win)
    if (win /= MPI_WIN_NULL) call MPI_WIN_FENCE(0, win, ierr)
    call MPI_BARRIER(shared_host_comm(), ierr)
  end subroutine gas_absorbed_power_sync

  !=========================================================================
  ! N_abs,s of the header for leaf il: the photoionizations [1/s] of H0,
  ! He0 and He+ (nabs(1:3)) in the field of the last reduced pass with the
  ! densities of the current gas state, D n_s Sum(path_lum sigma_s/(E erg))
  ! (V Gamma_s, V = (2 half D)^3, Gamma_s = Sum(...)/((2 half)^3 D^2)).
  ! valid is .false. (and nabs = 0) without the path sums of a reduced pass
  ! on the present grid: the first pass, after a restart, after a
  ! re-refinement (absorbed_field_valid).
  !=========================================================================
  subroutine previous_photoionizations(il, nabs, valid)
    use gas_state_mod, only : gas_nH, gas_xHI, gas_xHeI, gas_xHeII
    integer,       intent(in)  :: il
    real(kind=wp), intent(out) :: nabs(3)
    logical,       intent(out) :: valid
    nabs = 0.0_wp
    valid = absorbed_field_valid .and. il >= 1 .and. il <= score_nleaf
    if (.not. valid) return
    nabs(1) = gas_nH(il)*gas_xHI(il)*hhe_rate_path(I_GHI,il)
    nabs(2) = gas_nH(il)*par%abund_He*gas_xHeI(il)*hhe_rate_path(I_GHEI,il)
    nabs(3) = gas_nH(il)*par%abund_He*gas_xHeII(il)*hhe_rate_path(I_GHEII,il)
    nabs = nabs*par%distance2cm
  end subroutine previous_photoionizations


  !--- the ledger sums over ranks, in double precision (the sums of each rank
  !--- were accumulated in the extended kind)
  subroutine ledger_reduce()
    use mpi
    real(kind=wp) :: buf((4 + LEDGER_NABS)*LEDGER_NSRC)
    integer :: ierr, n
    n = LEDGER_NSRC
    buf(1:n)         = real(ledger_emitted, wp)
    buf(n+1:2*n)     = real(ledger_escaped, wp)
    buf(2*n+1:3*n)   = real(ledger_lost, wp)
    buf(3*n+1:4*n)   = real(ledger_roulette, wp)
    buf(4*n+1:)      = reshape(real(ledger_absorbed, wp), [LEDGER_NABS*n])
    call MPI_ALLREDUCE(MPI_IN_PLACE, buf, size(buf), MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    ledger_emitted  = real(buf(1:n), ENERGY_ACC_WP)
    ledger_escaped  = real(buf(n+1:2*n), ENERGY_ACC_WP)
    ledger_lost     = real(buf(2*n+1:3*n), ENERGY_ACC_WP)
    ledger_roulette = real(buf(3*n+1:4*n), ENERGY_ACC_WP)
    ledger_absorbed = real(reshape(buf(4*n+1:), [LEDGER_NABS, n]), ENERGY_ACC_WP)
    call sum_over_ranks(ledger_leaf)
  end subroutine ledger_reduce

  !--- one row for each source that emitted: emitted, absorbed by each
  !--- absorber, escaped [photons/s], and the closure
  !--- (absorbed + escaped)/emitted - 1.  Without dust scattering every
  !--- packet closes exactly (the edge walk scores e^-tau_in - e^-tau_out
  !--- in each leaf, escaped is e^-tau_edge), so the closure is rounding.
  !--- With it the scattered absorption is a path-length estimator and the
  !--- escape is the direct transmission plus the scattered flights that
  !--- leave the grid; the two close in expectation, and a finite sample
  !--- leaves Monte Carlo noise.  A lost flight (FLIGHT_LOST, reported with
  !--- the energy totals) escapes nowhere; it has its own column (lost), and
  !--- the closure is (absorbed + escaped + lost)/emitted - 1.  The roulette
  !--- column is Q = killed - gained in photons, zero in expectation and not
  !--- in the closure.  Under par%diffuse_launch = 'absorption_reemission'
  !--- the rows are the channel a packet carries at the time: a selected
  !--- re-emission is emitted in its new channel (its continued luminosity
  !--- over the new photon energy), and the photons absorbed in a channel
  !--- include those re-emitted later in another, so photon number is not
  !--- conserved across rows, only within each row in expectation.  The
  !--- header names the pass (label).
  subroutine ion_score_report_ledger(label)
    character(len=*), intent(in) :: label
    character(len=8), parameter :: SRC_NAME(0:LEDGER_NSRC-1) = &
       ['star    ', 'dif1 HI ', 'dif2 HeI', 'dif3HeII', 'dif4 HeI', &
        'dif5HeII', 'dif6 met', 'dif7 ff ', 'dif8 val', 'star FUV']
    integer :: is
    real(kind=wp) :: em, ab(LEDGER_NABS), es, ls, rq
    if (.not. par%absorption_ledger .or. mpar%p_rank /= 0) return
    write(*,'(a,5(3x,a6),3x,a,3x,a,3x,a,3x,a,3x,3a)') ' LEDGER: source   emitted      ', &
       LEDGER_ABS_NAME, 'escaped  ', 'lost     ', 'roulette ', &
       '(abs+esc+lost)/em-1', '[', trim(label), ']'
    do is = 0, LEDGER_NSRC-1
       em = real(ledger_emitted(is), wp)
       if (em <= 0.0_wp) cycle
       ab = real(ledger_absorbed(:,is), wp)
       es = real(ledger_escaped(is), wp)
       ls = real(ledger_lost(is), wp)
       rq = real(ledger_roulette(is), wp)
       write(*,'(a,a8,9es12.4e3,es11.2)') ' LEDGER: ', SRC_NAME(is), em, ab, es, &
          ls, rq, (sum(ab) + es + ls)/em - 1.0_wp
    end do
  end subroutine ion_score_report_ledger

  subroutine ion_score_apply_hhe(gHI, gHeI, gHeII, hHI, hHeI, hHeII)
    use octree_mod, only : leaf_half
    implicit none
    !--- inout: each rank writes its share of a node-shared array
    real(kind=wp), intent(inout) :: gHI(:), gHeI(:), gHeII(:)
    real(kind=wp), intent(inout) :: hHI(:), hHeI(:), hHeII(:)
    real(kind=wp) :: fac
    integer :: il

    if (.not. allocated(hhe_rate_path)) &
       call abort_run('H/He path estimators are not allocated')
    !--- node-shared output: the node's ranks take disjoint leaves
    do il = shared_host_rank() + 1, score_nleaf, shared_host_size()
       fac = 1.0_wp / ((2.0_wp*leaf_half(il))**3 * par%distance2cm**2)
       gHI(il)   = hhe_rate_path(I_GHI,il)*fac
       gHeI(il)  = hhe_rate_path(I_GHEI,il)*fac
       gHeII(il) = hhe_rate_path(I_GHEII,il)*fac
       hHI(il)   = hhe_rate_path(I_HHI,il)*fac
       hHeI(il)  = hhe_rate_path(I_HHEI,il)*fac
       hHeII(il) = hhe_rate_path(I_HHEII,il)*fac
    end do
  end subroutine ion_score_apply_hhe

  !--- Sum(path_lum) [erg/s * code length] of the last reduced pass in
  !--- leaf il over every band the gas iteration transports (band_path):
  !--- 4 pi J V integrated over the transported band, up to the factor of
  !--- ion_score_apply_hhe.  Positive exactly when a packet crossed the leaf.
  real(kind=wp) function band_path_luminosity(il) result(path)
    integer, intent(in) :: il
    path = band_path(il)
  end function band_path_luminosity

  subroutine ion_score_apply_dust_heat(heat_dust)
    !--- EUV grain heating [erg/s/cm^3] per leaf.  dust_abs_path
    !--- = Sum(path_lum * sigma_abs,dust(E)) with the factor below is the
    !--- transport form of
    !---   heat_dust = 4 pi J_abs rhokap / distance2cm,
    !--- since jt_ion = Sum path_lum gives 4 pi J dnu after the same factor.
    use octree_mod, only : leaf_half, amr_grid
    implicit none
    real(kind=wp), intent(out) :: heat_dust(:)
    real(kind=wp) :: fac
    integer :: il

    if (.not. allocated(dust_abs_path)) &
       call abort_run('dust-absorption path estimator is not allocated')
    do il = 1, score_nleaf
       fac = 1.0_wp / ((2.0_wp*leaf_half(il))**3 * par%distance2cm**2)
       heat_dust(il) = dust_abs_path(il)*fac &
                       * amr_grid%rhokap(il)/par%distance2cm
    end do
  end subroutine ion_score_apply_dust_heat

  !--- physical energy totals of the last reduced pass [erg/s]: emitted,
  !--- absorbed (in the gas and the grains), escaped (direct transmission
  !--- plus sampled flights leaving the grid) and lost; see the
  !--- definitions at the module variables.  Under absorption_reemission
  !--- (d) the emitted total also carries the sampled re-emission Re, every
  !--- packet emitted into the band as in the default mode, and absorbed is
  !--- the gross absorption emitted - escaped - lost = D + Q + Re.
  subroutine ion_score_energy_totals(emitted, absorbed, escaped, lost)
    real(kind=wp), intent(out) :: emitted, absorbed, escaped, lost
    emitted  = energy_emitted
    escaped  = energy_transmitted - energy_lost_direct + energy_scattered_escaped
    lost     = energy_lost_direct + energy_lost_scattered
    absorbed = energy_extinguished - energy_scattered_escaped - energy_lost_scattered
    if (reemission_on) then
       emitted  = emitted + sum(energy_reemitted)
       absorbed = absorbed + sum(energy_reemitted)
    end if
  end subroutine ion_score_energy_totals

  !--- the totals of absorption_reemission (d) of the last reduced pass
  !--- [erg/s]: reemitted = Re (L_REEMIT); remainder = the explicit
  !--- remainder launched from the leaf emissivity, the launched luminosity
  !--- of the diffuse channels 1-8 (L_LEAFEM); net_absorbed = L_ABS - Re =
  !--- D + Q (L_NETABS).  Optional: Re by channel 1-8.
  subroutine ion_score_reemission_totals(reemitted, remainder, net_absorbed, &
                                         reemitted_by_channel)
    real(kind=wp), intent(out) :: reemitted, remainder, net_absorbed
    real(kind=wp), intent(out), optional :: reemitted_by_channel(NDIFFUSE_CHANNEL)
    real(kind=wp) :: emitted, absorbed, escaped, lost
    call ion_score_energy_totals(emitted, absorbed, escaped, lost)
    reemitted    = sum(energy_reemitted)
    remainder    = sum(root_energy(ROOT_LAUNCHED,1:NDIFFUSE_CHANNEL))
    net_absorbed = absorbed - reemitted
    if (present(reemitted_by_channel)) reemitted_by_channel = energy_reemitted
  end subroutine ion_score_reemission_totals

  !=========================================================================
  ! Accessors for the output of an absorption_reemission pass (the rates
  ! file, gas_rates_mod), valid after ion_score_reduce and until the next
  ! ion_score_reset; every rank holds the reduced values.  [erg/s]
  !=========================================================================
  !--- leaf arrays (absorption_reemission only; size(.) = the leaf count):
  !--- measured = the gas absorption A_i measured in the pass on its own
  !--- state (GasAbsorbedPower).  The a_i the pass split its emission with
  !--- is diffuse_mod::diffuse_absorbed_power_leaf.
  subroutine ion_score_absorbed_power_leaf(measured)
    real(kind=wp), intent(out) :: measured(:)
    if (.not. reemission_on) call abort_run('ion_score_absorbed_power_leaf: '// &
       'par%diffuse_launch = ''absorption_reemission'' only')
    measured = gas_absorbed_now(1:size(measured))
  end subroutine ion_score_absorbed_power_leaf

  !--- reemitted(c, i): sampled re-emission of absorbing leaf i into channel
  !--- c = 1..8 (Sum = Re = L_REEMIT); event_gas(i): Sum G_event of leaf i
  subroutine ion_score_reemission_leaf(reemitted, event_gas)
    real(kind=wp), intent(out) :: reemitted(:,:), event_gas(:)
    if (.not. reemission_on) call abort_run('ion_score_reemission_leaf: '// &
       'par%diffuse_launch = ''absorption_reemission'' only')
    reemitted = reemitted_leaf(:,1:size(reemitted,2))
    event_gas = event_gas_leaf(1:size(event_gas))
  end subroutine ion_score_reemission_leaf

  !--- totals of the pass: sources = launched source luminosity (roots 0
  !--- and 9), launched = sources + explicit remainder (energy_emitted),
  !--- extinguished and transmitted of the direct split, deposited = D,
  !--- roulette = Q (killed - gained), continued_escaped = Cesc,
  !--- continued_lost = Clost, event_gas = Sum G_event, event_dust =
  !--- Sum Dd_event.  Defined in both modes (D, Q, Cesc, Clost are those of
  !--- the dust-scattered histories in the default mode).
  subroutine ion_score_interaction_energy_totals(sources, launched, extinguished, &
       transmitted, deposited, roulette, continued_escaped, continued_lost, &
       event_gas, event_dust)
    real(kind=wp), intent(out) :: sources, launched, extinguished, transmitted
    real(kind=wp), intent(out) :: deposited, roulette, continued_escaped, continued_lost
    real(kind=wp), intent(out) :: event_gas, event_dust
    sources           = root_energy(ROOT_LAUNCHED,0) + root_energy(ROOT_LAUNCHED,LEDGER_NSRC-1)
    launched          = energy_emitted
    extinguished      = energy_extinguished
    transmitted       = energy_transmitted
    deposited         = energy_deposited
    roulette          = energy_roulette
    continued_escaped = energy_scattered_escaped
    continued_lost    = energy_lost_scattered
    event_gas         = energy_event_gas
    event_dust        = energy_event_dust
  end subroutine ion_score_interaction_energy_totals

  !--- the photon ledger by the channel a packet carried (0..LEDGER_NSRC-1)
  !--- [photons/s] (par%absorption_ledger): emitted (launched and, under
  !--- absorption_reemission, re-emitted), absorbed by each absorber (H0,
  !--- He0, He+, metals, dust), escaped, lost
  subroutine ion_score_ledger_channels(emitted, absorbed, escaped, lost)
    real(kind=wp), intent(out) :: emitted(0:LEDGER_NSRC-1)
    real(kind=wp), intent(out) :: absorbed(LEDGER_NABS,0:LEDGER_NSRC-1)
    real(kind=wp), intent(out) :: escaped(0:LEDGER_NSRC-1), lost(0:LEDGER_NSRC-1)
    emitted  = real(ledger_emitted, wp)
    absorbed = real(ledger_absorbed, wp)
    escaped  = real(ledger_escaped, wp)
    lost     = real(ledger_lost, wp)
  end subroutine ion_score_ledger_channels

  !--- Printed totals of an absorption_reemission pass (rank 0), after
  !--- ion_score_reduce.  state_by_channel and remainder_by_channel are the
  !--- frozen-state emission L_c and the explicit remainder of the pass
  !--- (diffuse_mod), passed in since this module is below diffuse_mod.
  !--- The rounding-level identity is the root closure; the gross gas and
  !--- dust absorption of the path estimators against the event sums, and the
  !--- re-emission against its conditional mean, agree in expectation only.
  subroutine ion_score_report_reemission(state_by_channel, remainder_by_channel)
    use octree_mod, only : amr_grid
    real(kind=wp), intent(in) :: state_by_channel(NDIFFUSE_CHANNEL)
    real(kind=wp), intent(in) :: remainder_by_channel(NDIFFUSE_CHANNEL)
    character(len=8), parameter :: ROOT_NAME(0:LEDGER_NSRC-1) = &
       ['star    ', 'dif1 HI ', 'dif2 HeI', 'dif3HeII', 'dif4 HeI', &
        'dif5HeII', 'dif6 met', 'dif7 ff ', 'dif8 val', 'star FUV']
    real(kind=wp) :: emitted, absorbed, escaped, lost, reem, rem, net
    real(kind=wp) :: path_gas, path_dust, closed
    integer :: ir
    if (.not. reemission_on) return
    call ion_score_reemission_totals(reem, rem, net)
    call ion_score_energy_totals(emitted, absorbed, escaped, lost)
    path_gas = 0.0_wp
    if (allocated(gas_absorbed_now)) path_gas = sum(gas_absorbed_now(1:score_nleaf))
    path_dust = 0.0_wp
    if (allocated(dust_abs_path)) &
       path_dust = sum(dust_abs_path(1:score_nleaf)*amr_grid%rhokap(1:score_nleaf))
    if (mpar%p_rank /= 0) return
    write(*,'(a,8es11.3)') ' REEMIT: state emission L by channel 1-8 [erg/s]    :', &
       state_by_channel
    write(*,'(a,8es11.3)') ' REEMIT: explicit remainder X by channel 1-8 [erg/s]:', &
       remainder_by_channel
    write(*,'(a,8es11.3)') ' REEMIT: sampled re-emission Re by channel 1-8      :', &
       energy_reemitted
    write(*,'(a,es14.6,a,es14.6,a,i0,a,i0)') &
       ' REEMIT: Re = ', reem, ',  expected R_exp = ', energy_reemission_expected, &
       ';  interactions ', interaction_count, ', re-emissions ', reemission_count
    write(*,'(a,f12.8,a,es12.4)') ' REEMIT: largest continuation factor t of an interaction = ', &
       max_continuation, ';  largest weight in a history with t > 1 = ', max_gain_weight
    write(*,'(a,es14.6,a,es14.6,a,es11.3)') &
       ' REEMIT: gas absorption, path sum A = ', path_gas, ',  event sum G = ', &
       energy_event_gas, ';  G/A - 1 = ', ratio_dev(energy_event_gas, path_gas)
    !--- the leaf arrays of the output against the totals (rounding)
    write(*,'(a,es11.3,a,es11.3)') &
       ' REEMIT: leaf sums of Re, G_event over the totals - 1 = ', &
       ratio_dev(sum(reemitted_leaf), reem), ', ', &
       ratio_dev(sum(event_gas_leaf), energy_event_gas)
    if (par%ion_add_dust) write(*,'(a,es14.6,a,es14.6,a,es11.3)') &
       ' REEMIT: dust absorption, path sum   = ', path_dust, ',  event sum Dd = ', &
       energy_event_dust, '; Dd/path - 1 = ', ratio_dev(energy_event_dust, path_dust)
    write(*,'(a,es14.6,a,es14.6,a,es14.6)') &
       ' REEMIT: deposited D = ', energy_deposited, ',  roulette Q = ', &
       energy_roulette, ',  net absorbed D+Q = L_NETABS = ', net
    write(*,'(a,es14.6,a,es14.6,a,es14.6)') &
       ' REEMIT: L_EMIT = ', emitted, ' = launched ', energy_emitted, ' + Re ', reem
    write(*,'(a,4es14.6)') &
       ' REEMIT: L_ABS / L_ESC / L_LOST / L_LEAFEM = ', absorbed, escaped, lost, rem
    write(*,'(a)') ' REEMIT: by root   launched      escaped       lost'// &
       '          deposited     roulette      re-emitted    (esc+lost+D+Q)/launched-1'
    do ir = 0, LEDGER_NSRC-1
       if (root_energy(ROOT_LAUNCHED,ir) <= 0.0_wp) cycle
       closed = root_energy(ROOT_ESCAPED,ir) + root_energy(ROOT_LOST,ir) &
                + root_energy(ROOT_DEPOSITED,ir) + root_energy(ROOT_ROULETTE,ir)
       write(*,'(a,a8,6es14.6,es11.2)') ' REEMIT: ', ROOT_NAME(ir), &
          root_energy(:,ir), closed/root_energy(ROOT_LAUNCHED,ir) - 1.0_wp
    end do
  contains
    real(kind=wp) function ratio_dev(a, b)
      real(kind=wp), intent(in) :: a, b
      ratio_dev = 0.0_wp
      if (b > 0.0_wp) ratio_dev = a/b - 1.0_wp
    end function ratio_dev
  end subroutine ion_score_report_reemission

  subroutine ion_score_report_energy_closure()
    !--- Two energy statements of the transported packets.
    !--- direct split: emitted = extinguished + transmitted per packet by
    !---    construction (each packet's luminosity is split by its edge
    !---    transmission), so the residual measures only the accumulation
    !---    error of the reduced sums; past 5e-13 the run aborts.  With dust
    !---    scattering "extinguished" includes the scattered-out energy, and
    !---    without it it is the absorbed energy.
    !--- physical: absorbed / escaped / lost (ion_score_energy_totals), the
    !---    quantities the rates file carries as L_ABS / L_ESC / L_LOST.
    !--- With par%absorption_ledger the absorbed power of the path estimators
    !--- the rates use is compared with the physical absorbed total; the two
    !--- agree exactly without scattering and to Monte Carlo noise with it.
    implicit none
    real(kind=wp) :: closure, emitted, absorbed, escaped, lost

    closure = abs(energy_extinguished+energy_transmitted-energy_emitted) / &
              max(energy_emitted,tiny(1.0_wp))
    call ion_score_energy_totals(emitted, absorbed, escaped, lost)
    if (mpar%p_rank == 0) then
       write(*,'(a,2(i0,1x))') &
          ' ION: continuous scored path segments direct/sampled = ', &
          direct_segments, sampled_segments
       write(*,'(a)') ' ION: continuous H/He direct rates: ACTIVE'
       write(*,'(a,3es14.6)') &
          ' ION: continuous energy direct split emitted/extinguished/transmitted = ', &
          energy_emitted, energy_extinguished, energy_transmitted
       write(*,'(a,es11.3)') &
          ' ION: continuous energy direct split closure = ', closure
       write(*,'(a,3es14.6,a,i0,a)') &
          ' ION: continuous energy absorbed/escaped/lost = ', &
          absorbed, escaped, lost, '  (', lost_flights, ' lost flights)'
       !--- the forced-interaction balance (c), when the pass sampled
       !--- interactions
       if (energy_forced > 0.0_wp) then
          write(*,'(a,5es14.6)') &
             ' ION: interactions forced/deposited/roulette/escaped/lost = ', &
             energy_forced, energy_deposited, energy_roulette, &
             energy_scattered_escaped, energy_lost_scattered
          write(*,'(a,es11.3,a,es11.3)') &
             ' ION: (D + Cesc + Clost + Q)/forced - 1 = ', &
             (energy_deposited + energy_scattered_escaped + energy_lost_scattered &
              + energy_roulette)/energy_forced - 1.0_wp, &
             ';  forced/extinguished - 1 = ', &
             energy_forced/max(energy_extinguished, tiny(1.0_wp)) - 1.0_wp
       end if
       if (par%absorption_ledger) write(*,'(a,es14.6,a,es11.3)') &
          ' ION: continuous energy absorbed by the path estimators = ', &
          energy_path_absorbed, ';  minus absorbed, over emitted = ', &
          (energy_path_absorbed - absorbed)/max(emitted, tiny(1.0_wp))
       if (lost > 0.0_wp) write(*,'(a,es11.3,a)') &
          ' ION: WARNING: ', lost/max(emitted, tiny(1.0_wp)), &
          ' of the emitted energy was lost (MAX_WALK_STEPS or a failed'// &
          ' xy-periodic re-entry); the path estimates of those walks are partial.'
    end if
    if (closure > 5.0e-13_wp) then
       block
         character(len=160) :: reason
         write(reason,'(a,es11.3,a)') 'the continuous energy direct split closure ', &
            closure, ' exceeds 5e-13 (accumulation or reduction error).'
         call abort_run(trim(reason))
       end block
    end if
  end subroutine ion_score_report_energy_closure


  subroutine ion_score_apply_sec(edep)
    !--- Photoelectron energy deposition rate above E_SECONDARY [erg/s per
    !--- absorber particle], edep(k, a, leaf) for energy node k and absorber
    !--- a = 1 (H I), 2 (He I), 3 (He II): sec_ion_path with the factor that
    !--- turns Sum(path_lum coeff) into Int 4 pi J_nu coeff dnu, as for the
    !--- photoheating.  Summed over k it is the part of heat_* whose
    !--- photoelectrons exceed E_SECONDARY.
    use octree_mod, only : leaf_half
    implicit none
    real(kind=wp), intent(inout) :: edep(:,:,:)   ! node-shared, this rank's share
    real(kind=wp) :: fac
    integer :: il

    if (.not. allocated(sec_ion_path)) &
       call abort_run('secondary-ionization path estimators are not allocated')
    !--- node-shared output: the node's ranks take disjoint leaves
    do il = shared_host_rank() + 1, score_nleaf, shared_host_size()
       fac = 1.0_wp / ((2.0_wp*leaf_half(il))**3 * par%distance2cm**2)
       edep(:,:,il) = reshape(sec_ion_path(:,il)*fac, [N_EDEP, 3])
    end do
  end subroutine ion_score_apply_sec


  subroutine ion_score_report_metal_rates()
    if (.not. par%use_metals) return
    if (mpar%p_rank == 0) &
       write(*,'(a)') ' ION: continuous metal direct rates: ACTIVE'
  end subroutine ion_score_report_metal_rates

end module ion_score_mod
