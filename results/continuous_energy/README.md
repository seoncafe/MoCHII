# Continuous-energy validation outputs

This directory holds generated HII20/HII40 output.  Git ignores the outputs
(HDF5, logs, line and continuum tables, binaries); the inputs, drivers,
comparison scripts and their text and PDF summaries are tracked.

- `m2_atomic/`: the 2^25-packet Sobol production outputs of 2026-09-24
  (after the atomic-data update), superseded by `m4_nitrogen/`,
  `m5_converged/`, `m8_defaults/` and `m10_drift_test/`; the paper
  reads `m10_drift_test/` (since 2026-09-26) and `m9_hydrogenic_alpha/` for
  `g2a_thermal` and the G1 figure.  With them are the figure runs of
  that date
  (`pdr_L5_atomic`, `hei_high_density_atomic`, `g2a_thermal_atomic`) and the
  sensitivity runs (`*_nohecx`, `*_lowden`, `*_ots`); `run_lart3*.sh` are the
  drivers.
- `m5_converged/`: the paper set run to the drift stopping test
  (`conv_crit = 'vol'`, `fixed_point_distance_mod`) after the energy-ledger
  fix, with the HII40 inputs on `tests/g2_hii/hii40_L6_rout5p2.fits` (gas to
  5.2 pc); `pdr_L5_converged` uses the Sobol launch (the Mersenne-Twister run
  sat at a He noise floor above gas_tol and is kept as
  `pdr_L5_converged_random.*`); `g2a_thermal` uses the cell test and is not
  rerun here (`m4_nitrogen/`).  The paper read this set until 2026-09-26.
- `m10_drift_test/` (2026-09-26): the `m9_hydrogenic_alpha` set
  (`g2a_thermal` excepted, which uses the cell test) continued to
  convergence under the revised stopping test (`fixed_point_distance_mod`:
  rho fitted over 8 iterations, the Monte Carlo floor accepted only when the
  signed change of the ionized volume has no resolvable mean).  Each input
  is its m9 counterpart with the restart from the m9 state (HII40 from
  `hii40_long`, 16 iterations past `hii40_converged`); physics and
  tolerances unchanged.  `run_m10.sh` is the driver (one job per node from a
  shared queue: lart1 and lart2 at 72 ranks, mocafe at 88).
  `compare_m10.py` writes `compare_m10.txt` and `m10_profiles.pdf`.
  Iterations (m9 -> m10): HII40 4 -> 24, `ots` 4 -> 8, `lowden` 5 -> 36,
  `nohecx` 4 -> 43; HII20 26 -> 15, 4 -> 8, 23 -> 10, 28 -> 14; `pdr_L5`
  4 -> 16, `hei_high_density` 4 -> 28.  What moved: the HII40 H front
  +0.012 pc (L(Hbeta) +0.69%, <T[NpNe]> +0.2 K; `lowden` and `nohecx` the
  same), the on-the-spot runs +0.011 / +0.005 pc (+9.8 / +6.6 K, L(Hbeta)
  +0.54 / +0.52%), `hei_high_density` ionized volume +0.43% and
  sum V n_e +1.0%, `pdr_L5` <T> +5 K; HII20 +0.02% (m9 had converged).
  Paper-table values: <T[NpNe]> 8179.7 / 6940.1 K, L(Hbeta) 2.0226e37 /
  4.8456e36, <He+>/<H+> 0.7524 / 0.0472, R_out 1.4638e19 / 8.9875e18 cm.
  `g1_g4/` (`run_g1_g4.sh`): the level-7 G1 run and the two G4 runs
  (`tests/g4_refine/` inputs, with `helium_charge_exchange = .false.`
  pinned as in G1) with the same binary: R_eff 3.0578 (G1 level 7),
  3.0579 (re-refined, 467272 leaves) and 3.0577 pc (native level 7),
  against the 1D reference 3.0355 pc.
- `m9_hydrogenic_alpha/` (2026-09-26): the m8 set restarted with alpha_B of
  H II and He III from the Storey & Hummer (1991) cascade
  (`hydrogenic_recombination_mod`, -1.13% for H II at 1e4 K), the joint
  H0/He0 solve of the He II on-the-spot cascade and the nebular-continuum
  table gaps filled; `run_m9.sh`, `compare_m9.py` / `compare_m9.txt` /
  `m9_profiles.pdf` (m8 vs m9: HII40 -27.3 K, HII20 -17.8 K); `g1/` the
  G1 Stromgren reruns and `hii40_long` HII40 continued for 16 fixed
  iterations, which showed the stopping test of 2026-09-25 had stopped 8 of
  these runs early.  `g2a_thermal_m9` is the run the G2a figure reads.
- `g2a_current/` and `g2a_transported/` (2026-09-26): G2a with the
  on-the-spot and the transported recombination continuum, from the same
  state; in the optically thin He III core the on-the-spot cascade heat
  pins T at te_max (5e4 K), while transported the core sits at ~3.1e4 K.
- `m8_defaults/` (2026-09-26): the `m5_converged` set (and `g2a_thermal` of
  `m4_nitrogen`) restarted from those states with every 2026-09-25 default
  in force: charge-exchange heating, the He I 2^3S coupling, the band below
  the H edge from 4.78 eV (`add_fuv` following the 2^3S population, FUV
  packets as a separate stratum of no_photons/4), collider-resolved
  two-photon quenching, He I/He II on-the-spot photons, the He II diffuse
  channel with the transported continuum, the Mg + H+ temperature-dependent
  defect, and the diffuse packet normalization.  `run_m8.sh` is the driver
  (24 ranks, three at a time).  **This is the set the paper reads**
  (superseded by `m10_drift_test/`)
  (`g2a_thermal` excepted, below).  `compare_defaults.py` writes
  `compare_defaults.txt` and `defaults_profiles.pdf` (m5 vs m8, transported
  and on-the-spot, Cloudy overlaid); `front_zoom.py` writes `front_zoom.pdf`
  (T_e, x(H0), x(C0), x(O0) across the H front).
  - **Convergence:** HII40 8 more iterations, HII20 18, `*_ots` 17 / 9,
    `*_lowden` 8 / 17, `*_nohecx` 9 / 19, `pdr_L5` 18, `hei_high_density`
    20.  `g2a_thermal_defaults` (cell test) did not converge in 100, as in
    `m4_nitrogen`.
  - **Transported (production), m5 -> m8:** <T[NpNe]> 8206.72 -> 8206.82 K
    (HII40), 6946.07 -> 6957.91 K (HII20); L(Hbeta) +0.05% / +0.04%;
    <He+>/<H+> (Ercolano Eq. 17) 0.7582 -> 0.7584 / 0.0459 -> 0.0476; He
    front (x(He+) = 0.5, 0.01 pc shells) 4.2711 -> 4.2719 pc (Cloudy 4.305)
    / 1.2413 -> 1.2471; R_out 1.4596 -> 1.4599e19 / 8.9567 -> 8.9628e18 cm.
    m6 measured the three switches without the band at -7.65 K / -0.70 K,
    so the band below the H edge is about +8 K / +12 K.
  - **Neutral carbon now matches Cloudy.**  MoCHII/Cloudy x(C0), from
    volume-weighted 0.02 pc shells of the gas leaves against Cloudy's
    `.ele_c` interpolated to the shell centers: in m8, 0.99-1.35 over
    1.0-4.27 pc and 0.93-1.23 from the He front to 4.70 pc (HII40), 0.99-1.19
    over 1.0-1.25 pc and 0.95-1.06 from 1.25 to 2.88 pc (HII20).  In m5,
    without photons below 13.6 eV, it was 1.66-4.65 and 3.03-3.75 over the
    same inner ranges, rising to 257 (HII40) and 389 (HII20) just inside
    the H front.  [C I] 9824+9850 falls by 65x (HII40) and
    90x (HII20).
  - **The partly neutral edge is warmer**, by up to ~1000 K in a thin shell
    just outside the front (photoionization heating of C0, S0, etc.):
    [O I] 6300+6363 +14% / +22%, [N I] 5198+5200 +20% / +25%, [C II] 158 um
    +21% / +6%, [S II] 6716+6731 +4.0% / +2.1%, [S III] unchanged.  The
    other switches move no line by more than 1.6% (m6), so these are the
    band's.
  - **He I 5876 at HII20 +4.0%** with <He+>/<H+> +3.7%: at 20 kK the FUV
    luminosity is 14.6x the ionizing one, and 2^3S photoionization returns
    part of the triplet recombinations to He+ (inferred from the matching
    sizes; not isolated by a run).
  - **On the spot, m5 -> m8:** HII40 7961.65 -> 8149.01 K, L(Hbeta) +7.82%,
    H front 4.551 -> 4.713 pc, now 0.02 pc inside the transported run's
    (was 0.18); HII20 6888.25 -> 6929.31 K.  Mostly `hei_on_the_spot`
    (m7: +178 K alone at HII40).
  - **Zero-density cooling (`*_lowden`):** 7907.40 -> 7906.97 K /
    6294.56 -> 6303.62 K.
  - **`g2a_thermal_defaults`:** the paper's G2a gate against its 1D reference
    fails on this run (median |dev| 5.2%, mean bias +5.5%) and passes on
    `m4_nitrogen/g2a_thermal_nitrogen` (0.61%).  The 1D reference solves
    case B with the He I photons returned to helium, the He III case-B
    recombination (`alphaB_HeIII`) with no share of the He II photons given
    to H0 or He0 (no `heii_on_the_spot`), and neither the 2^3S coupling nor
    the band below the H edge, so the paper keeps the m4 run for that
    figure.
  - **`pdr_L5`:** H front 2.129 -> 2.191 pc, L(Hbeta) +5.9% (He I
    on-the-spot photons); the neutral-zone T_e profile (252 -> 172 K) is
    unchanged to 1 K; [C II] 158 um/Hbeta -0.5%, [O I] 63 um/Hbeta -8.8%.
  - **`hei_high_density`:** max N(2^3S) 3.53e15 -> 2.87e15 cm^-2, max
    tau_0(10833) 4.26e3 -> 3.46e3, L(Hbeta) +6.4%, L(10833)/L(Hbeta)
    1.387 -> 1.341.
- `m7_hei_on_the_spot/` (2026-09-25): the on-the-spot Lexington inputs
  (`m5_converged/*_ots_converged`) restarted with `par%hei_on_the_spot` on
  and off. Both use the new defaults for the other three switches of that
  day. `compare_hei_on_the_spot.py` writes `compare_hei_on_the_spot.txt`
  and `hei_on_the_spot_profiles.pdf`, overlaying the transported-continuum
  runs of `m6_new_physics` and Cloudy c25.00.
  - **HII40, on - off:** <T[NpNe]> +178.4 K (7959.0 -> 8137.3), L(Hbeta)
    +7.77%, H front +0.160 pc (4.550 -> 4.710), He front -0.021 pc.
  - **HII40, closer to the transported run:** the offset shrinks from
    -240 K / -0.183 pc (H front) to -62 K / -0.023 pc. [N II] 6584 is +21%
    and [O I] 63 um -28%, both toward the transported values.
  - **HII20, on - off:** +20.8 K, L(Hbeta) +0.36%. The offset from the
    transported run shrinks from -53 K to -32 K, and He front from -0.010 to
    -0.001 pc. T_e in the He+ zone rises by ~550 K.
  - **What is left is the approximation itself:** the gradient that remains
    in the H+ zone (HII20: -400 to +200 K against the transported run) is
    the same in both runs.
- `m6_new_physics/` (2026-09-25): the HII40/HII20 production inputs
  restarted from the `m5_converged` states (`par%restart_file`) with the three
  switches added that day. There are three kinds of run:
  - `*_new`: all three switches on, which is the new default;
  - `*_ctl`: the same binary with all three off (the control build);
  - `hii40_only_*`: HII40 with one switch on at a time.

  All runs converged within 4-6 iterations. `compare_new_physics.py` writes
  `compare_new_physics.txt` and `new_physics_profiles.pdf`, which shows
  T_e, x(H0) and x(He+) against Cloudy c25.00 and a Delta T panel.
  - **Totals, new - control:** <T[NpNe]> is -7.65 K at HII40
    (8206.87 -> 8199.22) and -0.70 K at HII20 (6946.10 -> 6945.40).
    L(Hbeta) changes by +0.044% / +0.004%. The fronts move by less than
    0.001 pc.
  - **Restart band:** control - m5 is +0.14 K / +0.02 K.
  - **Attribution at HII40, each switch alone:**
    - `hei_metastable_coupling`: -9.02 K, mostly 10830 excitation cooling in
      the He+ zone, a uniform -8 to -10 K inside 4 pc;
    - `charge_exchange_heating`: +1.36 K;
    - `twophoton_quenching`: 0 K (it only changes `_nebcont.txt`).

    The three add up to the total.
  - **Near the H front** the Delta T values of tens of K are inside the
    front's Monte Carlo band (control - m5 reaches +45 K there).
  - **Lines at HII40:** the temperature-sensitive UV and auroral lines drop
    by 1.0-1.6%, e.g. C III] 1909 -1.0% and [O III] 4363 -1.0%.
- `convergence/`: the HII40 production input run for a fixed 80 iterations,
  whose logged history calibrated the stopping test of 2026-09-25 (the
  extracted table was `tests/convergence/hii40_changes.txt`, removed
  2026-09-26 with the revised test, which also needs the signed change the
  log does not carry).
- `m4_nitrogen/`: the whole paper set (the two Lexington production runs,
  their `lowden`/`ots`/`nohecx` variants, and the `pdr_L5`, `g2a_thermal`
  and `hei_high_density` figure runs) rerun after the charge-exchange
  corrections of the night of 2026-09-24 (N0/N+ + H from Kingdon & Ferland
  1996, the form-7 Boltzmann factor at the actual T, He + H+ back to
  detailed balance); `compare_nitrogen.py` / `nitrogen_profiles.pdf`
  compare nitrogen with `m2_atomic/` and Cloudy c25.00, and
  `hii40_niter60_nitrogen.*` is the HII40 production input run for a fixed
  60 iterations (the convergence check recorded in CLAUDE.md).
- `m3_evening/`: the He I 10833 figure run and the two Lexington production
  runs of `m2_atomic/` rerun with the binary of the evening of 2026-09-24
  (He + H+ from Kimura et al. 1993; exact 2^3S collision constant and 2^3S
  ionization); `run_lart3.sh` and `run_bench.sh` are the drivers.
- `m1_oxfix/`: the production outputs before that update (oxygen
  charge-exchange correction).
- `production/`: an earlier production set, no longer read by the paper.
- `bin_independence/`: diagnostic-bin invariance runs.
- `rqmc_ensemble/`: independent scrambled-Sobol production ensemble.
- `mpi_scaling/` and `mpi_multinode/`: MPI agreement and timing runs.
- `smoke/` and `misc/`: short diagnostic outputs and logs.

Run scripts under `tests/continuous_energy/` write to the matching directory.
