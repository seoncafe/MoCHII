# Phase 4 diffuse-field items: attribution by control builds (2026-10-11/12)

Three Phase 4 items of
[ATOMIC_PHYSICS_RESTORATION_PLAN_2026-10-11.md](../../ATOMIC_PHYSICS_RESTORATION_PLAN_2026-10-11.md)
have no input switch, so each was removed by a minimal patch in a source
copy and run against the full tree on the 1e5 K blackbody case of the
Phase 4 record:

- (a) diffuse channel 6 (metal ground-recombination continua) removed;
- (b) diffuse channel 7 (free-free inside the band) removed;
- (c) the Ne III ground-term photoionization table (Nahar 2019, `PHTABLE`)
  removed;
- (d) the full tree.

Work directory: `../MoCHII_v1.00_phase5_attrib/` (beside the tree; deleted
on 2026-10-12, 8.0 GB, after the measurement; the patches, inputs, drivers
and scripts are kept in this folder, see Files). Nothing under `src/` or
`data/` of the tree was changed.

## Verdict

1. **The [S III] +0.78% of the Phase 4 record is channel 6, and it is a
   property of the 2-iteration state, not of the nebula.** Removing channel
   6 alone reproduces it from scratch: [S III] 9532 +0.77%, 33.5 um +0.89%
   (record: +0.66%, +0.78%), and a second realization (`qmc_seed` 777 for
   both builds) gives +0.53% / +0.65%. The carrier is the C III -> C II
   ground continuum (threshold 24.38 eV): in the 2-iteration state the gas
   is still ionized to the grid edge (left over from the initial guess), S+
   is the dominant sulfur stage beyond 2.9 pc, and the stellar field above
   the He0 edge has been absorbed there by x(He0) ~ 0.1. Restarted from a
   common state with a real H front (24 iterations), the same channel-6
   patch moves [S III] by +3e-5 / -4e-4 (9532 / 33.5 um), inside the seed
   noise (+2.7e-4), and the same again at seed 777.
2. **Channel 6 in the restarted state** acts in two places, both reproduced
   at both seeds and 5-20x the seed noise: (i) the O III, Ne III, O IV and
   C IV ground continua (54.9-77.4 eV) raise the He+ and Ne2+
   photoionization rates in the He2+ zone (Gamma(Ne2+) +1.3 to +3.1%),
   giving Ne3+ +0.66 / +0.70%, [Ne IV] 2422 +0.66 / +0.76%, [Ne V] 3426
   +0.56 / +0.61%, He II 4686 +0.13 / +0.18%; (ii) the C II ground continuum
   raises Gamma(S+) (and Gamma(C+), Gamma(O+)) near the H front by
   +0.5 to +0.8%, so the minor-stage lines there fall: [S II] 6716
   -0.31 / -0.37%, C II] 2325 -0.24 / -0.29%, [N II] 6584 -0.17 / -0.19%,
   [O II] 7320 -0.15 / -0.17%. L(Hbeta) does not move (<1e-4).
3. **Channel 7** is below every measurable level: 6.0e31 erg/s, 1.2e-7 of
   the diffuse field (2.5e31 erg/s, 2e-8, after 2 iterations). Removing it changes the lines by a median 5-6e-5,
   which is the response of the Monte Carlo realization to a shift of the
   diffuse emission CDF by 1e-8, not physics.
4. **The Ne III table** reproduces the Phase 4 record in the 2-iteration
   state (Ne3+ +1.73%, [Ne IV] +1.19%, [Ne V] 3426 +0.91%, [Ne III] 3869
   -0.34%; record +1.83 / +1.20 / +0.89 / -0.39%) and is larger in the
   restarted state (Ne3+ +2.19%, [Ne IV] +1.87%, [Ne V] +1.24%, [Ne III]
   -0.27%). It does not move the sulfur lines (<= 1.6e-4 restarted,
   inside the seed noise after 2 iterations).

No code defect was found. Two observations (not defects) are in
"Observations".

## Builds and patches

Source: the tree `src/` (73 `.f90`/`.inc` files and `Makefile.dep`) copied
at 2026-10-11 19:33 KST (`d/SNAPSHOT.txt`, `d/src_md5.txt`); data:
`data/` copied to `data_snapshot/` at the same time (identical to the tree
then). Since then the tree changed `define.f90`, `nlevel_mod.f90`,
`setup.f90` and `data/atomic/README.md` (md5 check 2026-10-12 01:00); these
runs are of the 19:33 source. Build `make MoCHII.x -j 16` in each copy,
mpiifort (oneAPI 2026.1, Intel MPI 2021.18), the Makefile's flags.

(a), `a/src/diffuse_mod.f90` (`patch_a.diff`):

```diff
@@ -186,7 +186,8 @@
        !--- channels 6 and 7: the metal ground continua and the free-free
        !--- emission inside the band
        if (all_channels) then
-          dif_ch(CH_METAL_GROUND,il) = ne*vol*sum(metal_ground_weights(il, T))
+          !--- control build (a): channel 6 removed
+          dif_ch(CH_METAL_GROUND,il) = 0.0_wp
           dif_ch(CH_FREEFREE,il)     = ne*vol*sum(freefree_weights(il, T))
        end if
```

(b), `b/src/diffuse_mod.f90` (`patch_b.diff`):

```diff
@@ -187,7 +187,8 @@
        if (all_channels) then
           dif_ch(CH_METAL_GROUND,il) = ne*vol*sum(metal_ground_weights(il, T))
-          dif_ch(CH_FREEFREE,il)     = ne*vol*sum(freefree_weights(il, T))
+          !--- control build (b): channel 7 removed
+          dif_ch(CH_FREEFREE,il)     = 0.0_wp
        end if
```

`dif_ch` is the only route by which channels 6 and 7 enter a run (the leaf
CDF, the packet count and the channel pick all read it; no other module
reads `metal_ground_weights` or `freefree_weights`), so a zero there
removes the channel and nothing else.

(c), the binary of (d) with `par%atomic_dir = atomic_noNe3/`: symbolic links
to every file of `data_snapshot/atomic/` except `element_ne.txt`, which is
the snapshot's with the record `PHTABLE 25013` and its 25013 rows deleted
(lines 62-25075; `patch_c.diff`, 25023 lines). The comment line that
documents the record is kept. Without the record the reader takes the
VFKY96/VY95 fits for Ne III, as for every other transition.

## Conditions

| item | value |
|---|---|
| input | `../MoCHII_v1.00_phase4_20261011/runs/bb1e5_p4.in` (the Phase 4 record's 1e5 K case), with `atomic_dir`/`data_dir` on the snapshot, `jnu_output` and `emis_output` on (outputs only) |
| geometry | `tests/g2_hii/hii40_L6_rout5p2.fits`, 262144 leaves of 0.1625 pc, n_H = 100 cm^-3 to 5.2 pc |
| source | 1e5 K blackbody, L = 9.3598e38 erg/s, band 13.598-150 eV, `add_fuv` on (default) |
| gas | He/H 0.1; C 2.2e-4, N 4.0e-5, O 3.3e-4, Ne 5.0e-5, S 9.0e-6; thermal balance on (3000-5e4 K); no dust |
| diffuse field | `recombination_continuum = 'transported'` (channels 1-7), current defaults otherwise |
| packets | 1048576 ionizing, Sobol launch, `qmc_seed` 20260728 (noise runs 777), `iseed` 1234 |
| host, ranks | lart3, 16 MPI ranks, one run at a time; `I_MPI_ADJUST_{ALLREDUCE,REDUCE,BCAST,GATHER}=1` |
| set `it2` | each build from the initial guess, 2 iterations + final pass (as the Phase 4 record, which ran 36 ranks) |
| set `rs` | each build restarted from `runs/common2/d` (build d, 12 + 12 iterations from the initial guess), 3 iterations + final pass |
| noise | build d at `qmc_seed` 777 (`noise`); build a at 777 (`a_noise`), so `noise/a_noise` is a second realization of the channel-6 pair |

**Why the restart set.** Two iterations are not a nebula. The state after
iteration 2 is ionized to the grid edge (x(H+) >= 0.98 to 5.2 pc, x(He0)
~ 0.1 and x(He2+) ~ 0.035 beyond 3 pc, T_e ~ 7500 K there): the initial
guess (x(He+) 0.5, x(He2+) 0.4 everywhere) fed a He III recombination
field (channels 3 and 5) of 1.3e39 erg/s into the first pass. Between iterations 3 and 11 the gas
beyond ~3.5 pc recombines (signed volume step of x(H+) down to -0.27 per
iteration), and after 24 iterations the H front (x(H+) = 0.5) stands at
~3.98 pc, still advancing by +1.6% of the ionized volume per iteration
(rho 0.91, distance to the fixed point ~0.26 at iteration 8 of common2). Over
the 3 restarted iterations build d itself moves L(Hbeta) +3.7%, [S II]
+11.5%, [N I] 5200 +15%, [S III] 9532 +0.84% (median line 3.6%). All
builds start from the same state and carry the same drift; the table
differences are taken at equal iteration count. The restarted state is
**not converged**; the item effects of a converged state were not
measured.

## Results

Ratios B/A - 1. For an item, A is the build without it and B the full
build (d); "seed noise" is noise/d; "ch. 6 at seed 777" is noise/a_noise.
N(X) = Sum V n_H x(X). Gamma(S+) is the code's rate summed with weight
V n_H in the shell (state independent; a weight by n(S+) of one run selects
that run's low-rate leaves and is biased by the noise: it gave +8% for the
noise pair). Line ratios are of L/L(Hbeta). Median over the 93 lines above
1e-3 L(Hbeta).

### Set `rs`: restarted from the common state (24 iterations)

| quantity | ch. 6 (d/a) | ch. 7 (d/b) | Ne III table (d/c) | seed noise (noise/d) | ch. 6 at seed 777 |
|---|---|---|---|---|---|
| L(Hbeta) | -4.89e-06 | +6.03e-06 | +1.02e-06 | -1.05e-04 | -5.58e-05 |
| <T_e> (n_e^2 V) | +3.60e-04 | +1.72e-05 | +4.22e-05 | +7.18e-05 | +4.18e-04 |
| <x(He2+)> (n_e^2 V) | +1.55e-03 | -3.26e-05 | +7.61e-05 | +2.18e-04 | +2.01e-03 |
| N(C2+) | -1.41e-05 | -6.67e-07 | +2.62e-05 | +4.39e-05 | +4.36e-05 |
| N(N2+) | -1.35e-04 | -8.18e-07 | -1.68e-05 | +8.74e-05 | -9.55e-05 |
| N(O2+) | +3.03e-04 | +1.81e-05 | +1.10e-05 | -2.58e-05 | +2.89e-04 |
| N(Ne3+) | +6.60e-03 | +5.94e-04 | +2.19e-02 | +4.57e-04 | +6.95e-03 |
| N(S+) | -4.30e-04 | +3.86e-05 | +2.21e-05 | -4.13e-05 | -4.93e-04 |
| N(S2+) | -2.91e-04 | -2.06e-04 | -1.18e-04 | +2.54e-04 | -3.45e-05 |
| Gamma(S+), 2.9-3.6 pc | +4.78e-03 | -1.42e-04 | -4.72e-04 | -8.11e-04 | +4.69e-03 |
| Gamma(S+), 3.6-4.1 pc | +6.92e-03 | -5.21e-06 | +3.28e-04 | +2.36e-04 | +8.19e-03 |
| [S III] 9532 | +2.74e-05 | -1.61e-04 | -7.70e-05 | +2.72e-04 | +2.59e-04 |
| [S III] 33.5 um | -4.30e-04 | -1.66e-04 | -7.64e-05 | +2.71e-04 | -2.34e-04 |
| [S III] 18.7 um | -3.60e-04 | -1.66e-04 | -7.51e-05 | +2.67e-04 | -1.61e-04 |
| [S II] 6716 | -3.05e-03 | +5.47e-04 | +1.60e-04 | -9.05e-04 | -3.73e-03 |
| [O III] 5007 | +9.77e-04 | +2.20e-05 | -1.02e-06 | +1.70e-04 | +1.13e-03 |
| [O II] 7320 | -1.51e-03 | +1.08e-05 | +1.11e-04 | -2.93e-04 | -1.66e-03 |
| [N II] 6584 | -1.70e-03 | +7.56e-05 | +5.23e-05 | -4.48e-04 | -1.89e-03 |
| [N I] 5200 | -3.06e-03 | +1.17e-03 | +7.57e-04 | -3.21e-03 | -5.42e-03 |
| [Ne III] 3869 | +4.27e-04 | -6.89e-05 | -2.68e-03 | +6.60e-05 | +5.37e-04 |
| [Ne IV] 2422 | +6.61e-03 | +4.03e-04 | +1.87e-02 | +1.25e-03 | +7.55e-03 |
| [Ne V] 3426 | +5.62e-03 | +2.01e-04 | +1.24e-02 | +9.64e-04 | +6.08e-03 |
| C II] 2325 | -2.37e-03 | +1.74e-04 | +1.13e-04 | -7.28e-04 | -2.90e-03 |
| He II 4686 | +1.31e-03 | -1.72e-05 | -6.58e-05 | +2.74e-04 | +1.75e-03 |
| He I 5876 | -1.87e-04 | +3.10e-06 | +8.11e-06 | -6.86e-05 | -2.91e-04 |
| median abs, lines | 1.7e-03 | 6.3e-05 | 1.1e-04 | 3.0e-04 | 1.9e-03 |

[N I] 5200 moves by the same amount in the noise pair; it is not
attributed.

### Set `it2`: each build from the initial guess, 2 iterations

| quantity | ch. 6 (d/a) | ch. 7 (d/b) | Ne III table (d/c) | seed noise (noise/d) | ch. 6 at seed 777 |
|---|---|---|---|---|---|
| L(Hbeta) | +1.58e-04 | -3.58e-07 | +1.49e-04 | -2.75e-04 | +1.32e-04 |
| <T_e> (n_e^2 V) | -4.60e-05 | +7.60e-07 | -1.38e-04 | +1.71e-04 | -6.43e-05 |
| <x(He2+)> (n_e^2 V) | -2.45e-04 | -1.48e-05 | -2.87e-04 | -6.16e-03 | +3.17e-04 |
| N(C2+) | +1.11e-02 | -5.06e-06 | +1.39e-03 | -4.56e-03 | +9.17e-03 |
| N(N2+) | +9.04e-03 | +4.49e-04 | +1.61e-03 | -2.58e-03 | +9.92e-03 |
| N(O2+) | +7.31e-03 | +4.35e-04 | +1.62e-03 | -2.57e-03 | +8.07e-03 |
| N(Ne3+) | +5.54e-04 | -1.35e-04 | +1.73e-02 | -5.17e-05 | +1.11e-03 |
| N(S+) | -3.70e-03 | +2.88e-05 | -3.66e-04 | +1.80e-03 | -2.74e-03 |
| N(S2+) | +8.89e-03 | -8.16e-05 | +7.95e-04 | -4.13e-03 | +6.43e-03 |
| Gamma(S+), 2.9-4.0 pc | +1.55e-02 | -1.36e-04 | +1.95e-03 | -4.33e-03 | +9.38e-03 |
| Gamma(S+), 4.0-5.3 pc | +2.17e-02 | -2.18e-04 | +2.79e-03 | -7.54e-03 | +1.87e-02 |
| [S III] 9532 | +7.69e-03 | -9.75e-05 | +3.08e-04 | -3.30e-03 | +5.28e-03 |
| [S III] 33.5 um | +8.94e-03 | -7.60e-05 | +7.21e-04 | -4.00e-03 | +6.50e-03 |
| [S III] 18.7 um | +8.84e-03 | -7.84e-05 | +6.89e-04 | -3.95e-03 | +6.40e-03 |
| [S II] 6716 | -4.14e-03 | +4.43e-05 | -8.90e-04 | +2.43e-03 | -3.19e-03 |
| [O III] 5007 | +4.66e-03 | +3.26e-04 | +8.11e-04 | -2.22e-03 | +5.09e-03 |
| [O II] 7320 | -3.33e-03 | -1.63e-04 | -1.73e-03 | +1.42e-03 | -4.62e-03 |
| [N II] 6584 | -1.84e-03 | -9.71e-05 | -8.94e-04 | +1.05e-03 | -2.20e-03 |
| [N I] 5200 | -6.07e-03 | +4.67e-04 | +1.99e-04 | +7.14e-03 | -4.16e-03 |
| [Ne III] 3869 | +1.63e-03 | +3.96e-04 | -3.39e-03 | -1.49e-03 | +2.62e-03 |
| [Ne IV] 2422 | +7.97e-04 | -3.46e-05 | +1.19e-02 | -5.68e-04 | +9.23e-04 |
| [Ne V] 3426 | +6.49e-04 | +3.58e-07 | +9.13e-03 | -1.29e-03 | +1.18e-03 |
| C II] 2325 | -4.79e-03 | +9.43e-06 | -1.73e-03 | +1.93e-03 | -5.46e-03 |
| He II 4686 | +1.57e-03 | -1.95e-05 | -3.67e-05 | -3.68e-03 | +2.24e-03 |
| He I 5876 | +7.79e-04 | +3.22e-06 | +1.80e-04 | -7.10e-04 | +4.73e-04 |
| median abs, lines | 2.4e-03 | 5.3e-05 | 5.6e-04 | 1.4e-03 | 2.6e-03 |

The seed noise of [S III] here is -0.33 / -0.40% (9532 / 33.5 um); the
Phase 4 record's seed run (36 ranks) gave -0.14% for 33.5 um.

## Mechanism of the [S III] change (set `it2`)

Leaf-level evidence, scripts in `measure/`:

1. *Where.* The [S III] 9532 increase is a sum over the He+ zone (99.5% of
   the line; `compare.py`): x(S2+) accounts for +0.92% of dln(epsilon)
   while n_e (+6e-5) and the rest (T_e, n_e inside the level solve, -2e-5)
   do not move. In radial shells 0.34 pc wide (`n_H V` weighted) x(S2+)
   changes by -0.18 to -0.01% inside 2.7 pc and by +0.3, +0.7, +1.5,
   +1.6, +1.7, +1.1, +2.1, +2.2% in the shells centered at 2.9 to 5.3 pc,
   where S+ is the dominant stage (0.34-0.87) and the gas is at 7500-8500
   K; 68% of the line comes from beyond 2.7 pc in this state.
2. *Rate.* Gamma(S+) in those shells rises +1.55% (2.9-4.0 pc) and +2.17%
   (4.0-5.3 pc) against seed noise -0.43 / -0.75%; at seed 777 +0.94 /
   +1.87%. With S+ dominant, dx(S2+)/x(S2+) ~ x(S+) dGamma/Gamma, which
   is what point 1 shows. T_e (-2e-5) and n_e (+6e-5, emissivity
   weighted) do not move; x(H0) (charge exchange S2+ + H) moves -5e-4.
3. *Energy.* Splitting the rate change over the J_nu bins with the VFKY96
   S+ cross section (`gamma_bins.py`; the estimate reproduces the code's
   rate change to 2-6% of itself): in 2.9-4.0 pc, +1.06e-2 from 23.33-24.58
   eV and +0.35e-2 from 24.58-25.51 eV; in 4.0-5.3 pc, +1.33e-2 and
   +0.57e-2. The 23.33-24.58 eV bin (between the S+ and He0 edges) carries
   +21.5% more mean intensity in 4.0-5.3 pc with channel 6 (`jshell.py`).
4. *Source.* The channel-6 luminosity of the written state by transition
   and window (`channel6/channel6_spectrum.x`, the code's own
   `milne_metal_band_coefficient` and `milne_metal_energy_integral` on the
   written stage fractions): in 2.9-4.0 / 4.0-5.3 pc the 23.34-24.59 eV
   window holds 3.1e34 / 5.1e34 erg/s, of which C III -> C II(ground)
   (threshold 24.383 eV) is 97%, S III -> S II 2.4 / 1.6%, Ne II -> Ne I
   0.9 / 1.0%; the 24.59-34.86 eV window holds 1.0e35 / 1.6e35 erg/s,
   C III -> C II 93%, N III -> N II 7%.
5. *Why so large there.* In this state the stellar field above 24.59 eV is
   absorbed by He0 (x(He0) ~ 0.1, i.e. tau ~ 20 per pc at threshold) and
   the field between 23.34 and 24.59 eV only by H0, so the C II ground
   continuum, made in place where C2+ is abundant, raises the mean
   intensity there by 21% (4.0-5.3 pc) and the S+ photoionization rate by
   1.5-2.2%.

**In the restarted state** the same source acts near the real H front:
Gamma(S+) +0.48% (2.9-3.6 pc) and +0.69% (3.6-4.1 pc), seed 777 +0.47 /
+0.82%, noise -0.08 / +0.02%, again from 23.3-27.5 eV (C II 99% of the
23.34-24.59 eV window, 87-91% of 24.59-34.86 eV). But S2+ is now the
dominant stage there (0.5-0.75, S+ 0.02-0.16 inside 3.9 pc), so the change
lands on S+, not on S2+: N(S2+) -2.9e-4 / -3.5e-5 (noise +2.5e-4), [S III]
inside the noise, [S II] 6716 -0.31 / -0.37%.

## Channel 6 in the He2+ zone (set `rs`)

Gamma(Ne2+) +1.73, +3.11, +1.34% in 0-1.8, 1.8-2.2, 2.2-2.9 pc (seed 777
+1.87, +2.73, +1.33%; noise +0.08, -0.58, +0.25%) and Gamma(He+) +0.45,
+1.58% in 0-1.8, 1.8-2.2 pc (seed 777 +0.46, +1.65%; noise +0.02, +0.08%)
(`gamma_shells.py`). Inside 2.9 pc channel 6 emits 5.3e35 erg/s above
54.42 eV: O IV -> O III (threshold 54.94 eV) 51%, Ne IV -> Ne III (63.42
eV) 26%, O V -> O IV (77.41 eV) 18%, C V -> C IV (64.49 eV) 4%. Result: N(Ne3+) +0.66 / +0.70%, [Ne IV] +0.66 / +0.76%,
[Ne V] +0.56 / +0.61%, <x(He2+)> +0.16 / +0.20%, He II 4686 +0.13 / +0.18%,
<T_e> +3.6e-4 / +4.2e-4.

## Observations (not defects)

- **Written stage fractions and the final pass's channel 6.** The
  channel-6 luminosity recomputed from the written state agrees with the
  printed final-pass value to 0.3% in the restarted state (1.9998e36
  against 1.993e36 erg/s) but not after 2 iterations (1.354e36 against
  2.017e36): the final pass builds channel 6 from the stage fractions of the
  cache refreshed after iteration 2 (`gas_opacity_fill`), while
  `species_write` recomputes the written fractions with
  `species_fractions` on the final pass's photoionization rates and the
  iteration-2 T_e, n_e and hydrogen state. Unconverged, the written metal
  stages therefore belong to a different field than the written n_e (which
  the previous fractions entered); the two coincide as the iteration
  converges. Recorded, not changed; whether the written state should be
  re-solved after the final pass is the author's call.
- **Leaf-level rates are noise dominated at 1M packets.** Gamma(S+) of
  single leaves differs between the two builds by a median +0.4% with 5-95%
  quantiles of -110% / +116% (it2); comparisons must be sums over regions
  with state-independent weights.

## Not checked

- Converged states: the restarted set is 24 + 3 iterations, the H front
  still advancing ~1.5% of the ionized volume per iteration.
- The HII40 benchmark case and PN150: not rerun here.
- Only one rank count (16) and one host (lart3); the Phase 4 record ran 36
  ranks.
- Channel 7 below the realization floor; its physical effect is not
  resolved by these runs (its luminosity, 3e-8 of the diffuse field, bounds
  it).
- The tree changes after 19:33 KST (`define.f90`, `nlevel_mod.f90`,
  `setup.f90`) are not in these builds.

## Files (this folder)

The work directory with the source copies, binaries, data snapshot and run
outputs (8.0 GB, mostly `_rates`/`_emis`) was deleted on 2026-10-12. Kept:

- `SNAPSHOT.txt`, `src_md5.txt`: the source snapshot time and md5 of every
  source of build d.
- `patch_a.diff`, `patch_b.diff`: builds (a) and (b) against d. Build (c)
  is d with `PHTABLE 25013` and its 25013 rows deleted from
  `element_ne.txt` (lines 62-25075 of the snapshot), nothing else.
- `runs/`: the inputs of every run (`it2_*.in`, `rs_*.in`,
  `common2_d.in`), `status.txt` (start and exit times), and the drivers
  `run_all.sh` (stopped after `common/d`), `run_rest.sh`, `run_extra.sh`
  (superseded by `run_rest.sh`). Their paths point into the deleted work
  directory.
- `measure/`: `summarize.py` (tables), `compare.py`, `gamma_bins.py`,
  `gamma_shells.py`, `jshell.py`, `profile.py`, `fronts.py`, `drift.py`,
  `gamma_dist.py`, `quick_old.py`; `channel6/` (diagnostic program
  `channel6_spectrum.f90`, `dump_state.py`, and its outputs `it2_d.txt`,
  `common2_d.txt`, `rs_d.txt`); `rs_table.md`, `it2_table.md`.
