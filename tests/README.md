# tests

Regression cases + baselines, MOCASSIN-suite style.  Each directory is
self-contained: grid (or its builder), input file(s), and a `check_*.py`
gate script.  Gates per stage are defined in md/PLAN.md section 5;
results are recorded in CLAUDE.md and docs/MoCHII_physics.pdf.

**Before comparing two runs bit for bit**, note that bit identity needs the MPI
collective algorithm selection held fixed.  Confirmed here: with no code change
at all, two runs of the same binary on the same input at the same rank count
already differ — at 16 ranks on `cooling_ne/id_new.in` only 5 of the 31
`rates.h5` datasets matched bit for bit, with `J_nu` differing by 5.7e-14
(`_lines.txt` was identical).  Intel MPI choosing its collectives adaptively is
the diagnosis; pinning them with `I_MPI_ADJUST_ALLREDUCE=1` (plus `REDUCE`,
`BCAST`, `GATHER`) is *reported* to restore bit identity, but that has not been
re-confirmed.  A mismatch at the 1e-14 level is therefore not by itself a
regression.

| directory | gate |
|---|---|
| `g0_gamma` | rate integrals vs analytic attenuation |
| `g1_stromgren` | Stromgren sphere vs the 1D Gauss-Seidel reference; AMR = uniform |
| `g2a_thermal` | H/He thermal balance vs the 1D reference |
| `g2_hii` | MOCASSIN Lexington HII20/HII40 + line fluxes; element smokes (`smoke_*.in`), emissivity output (`smoke_emis.in`), He I excited channel (`hii40_hei.in`) |
| `g4_refine` | I-front re-refinement vs native level 7 (+ window recycling) |
| `g4_tng` | Illustris-TNG post-processing demo |
| `io` | HDF5 first-error preservation, failed-create propagation, HDF5-off compilation, and case-insensitive extension detection |
| `thermal_parallel` | Exact `x_HI`, `x_HeI`, `x_HeII`, `n_e`, and `T_e` comparison between serial-thermal and node-local MPI-thermal rate files |
| `continuous_energy` | Continuous-energy gates: C II/He I thresholds; exact-energy H/He Planck/table transport and solver rates; exact metal ionization/heating, opacity, electron, and thermal coupling; threshold-sorted metal scoring benchmark; random/Sobol, 1/3-rank, iterative, and 8/32-bin independence |
| `energy_sampler` | Isolated continuous source samplers: analytic tabulated inverse CDF, adaptive band-limited Planck CDF, full-Planck Barnett--Canfield reference, random/Sobol statistics, source-coordinate mappings, and 1/2/3/5/8-rank QMC identity |
| `d_dusty` | dusty Stromgren, scattering bracket, T_dust/IR, SEDust smokes (`sedust_smoke.in`, `dustemis_smoke.in`, `pahlive_smoke.in`), FUV option (`d_fuv_*.in`, `check_fuv.py`) |
| `peel` | peel-off imaging: optically-thin analytic direct gate (`peel_thin.in`, `check_peel.py`), scattered morphology (`peel_scat.in`), bin cubes (`peel_cube.in`) |
| `uni_dda` | Cartesian (car) grid vs single-level octree: `car_walk='amr'` walk bit-identical; incremental DDA walk (`car_walk='dda'`) to rounding + timing (`check_uni.py`); namelist-built car grid (`car_namelist.in`) + `nH_const`/`rmax` density model reproduce the file-built sphere (`amr` override `car_amr_ovr.in`) |
| `pdr` | PDR physics: metal electrons, photoelectric heating, atomic-hydrogen impact cooling (`pdr_L5.in`, `check_pdr.py`) |
| `grain_kext` | the extinction SEDust serves (the `/kext` curve of `SEDust/data/<model>/sedust_<model>.h5`) against the size integral over the grain model built from the same product, for the four offered models (astrodust, DL07, MRN, Zubko): C_ext, C_abs, albedo and albedo x <cos> on each side of 0.0912 um, plus the band-averaged `s_abs` and `C_ext(lambda_ref)` the transport is built on (`check_kext_table.f90`, `run_check.sh`).  Curve and grid come from one product, so what is left is rounding: one tolerance, 1e-12, for every node |
| `heii_cascade` | the He III case-B cascade diffuse channel (`par%heii_diffuse`, `check_heii_cascade.f90`, `run_check.sh`): the n = 2 exit branching closes and both parts come from SH95; `alpha_eff(2^2P)` reproduces MOCASSIN's independent He II Ly-alpha emissivity fit to 2.3-6.1% over 5-30 kK; the three branch luminosities are positive and ordered Ly-alpha > two-photon > free-bound; every sampled energy lies in [13.598, 54.416] eV; the two-photon in-band moments match independent quadrature of the Nussbaumer & Schmutz (1984) shape to 1e-4, and `<1/E> <E>_n = 1` so equal-energy packets carry the physical photon rate |
| `photo_shells` | metal photoionization cross section and photoelectron heating of every registry ion against Verner's `phfit2.f` on 5-1000 eV: the VFKY96 outer-shell fit below the first inner-shell edge, the sum over occupied subshells (VY95) above it (`check_photo_shells.f90`, `run_check.sh`; `VERNER_DIR` points at the `phfit2.f` tree) |
| `hcoll` | atomic-hydrogen impact on the C I, C II and O I ground terms (`data/atomic/hcoll_*.txt` in `nlevel_mod`): the low-density photon rates reproduce Abrahamsson et al. (2007) Tables 1 and 3 and the Barinovs et al. (2005) C II law, the populations reach Boltzmann at n_HI = 1e14 on both sides of the 1000 K fit edge, and the thermal increment equals the emissivity difference (`check_hcoll.f90`, `run_check.sh`) |
| `convergence` | the `conv_crit = 'vol'` stopping test of `fixed_point_distance_mod`: replayed on the recorded HII40 history (`hii40_changes.txt`, the production input run for a fixed 80 iterations) it stops at iteration 58, where the state equals the 80-iteration state; a geometric decay passes exactly when the drift still to come d rho/(1-rho) is below tol; a slow decay never passes on its step; a noise floor below tol passes (`check_fixed_point.f90`, `run_check.sh`) |
| `restart` | `par%restart_file` (`previous_state_mod`): a short HII20 run writes a state; a restart on the same grid with `gas_niter = 0` returns x_HI, x_HeI, x_HeII, n_e, T_e bit for bit, and a restart on a 32^3 grid gives every gas leaf the values of the previous leaf containing its centre (`run_check.sh [nranks]`, `check_restart.py`) |
| `secondary_ionization` | the fast-photoelectron partition of `photoelectron_degradation_mod` (`par%use_sec_ion`): the 1 keV limit at SvS85's composition reproduces Shull & van Steenberg (1985) Table 2, energy is conserved at any He/H, the heat share is Dalgarno, Yan & Liu (1999) Table 7's eta0 at x -> 0 and 1 at x = 1, and the energy weights interpolate linearly in log E0 (`check_photoelectron_partition.f90`, `run_check.sh`) |
