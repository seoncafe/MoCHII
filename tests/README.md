# tests

Regression cases + baselines, MOCASSIN-suite style.  Each directory is
self-contained: grid (or its builder), input file(s), and a `check_*.py`
gate script.  Gates per stage are defined in docs/PLAN.md section 5;
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
| `pdr` | PDR physics: metal electrons, photoelectric heating, H-impact cooling (`pdr_L5.in`, `check_pdr.py`) |
| `grain_kext` | the extinction SEDust serves from `data/kext_*.dat` against the size integral over the same grain model, on MoCHII's default EUV Q grid (astrodust at `eion_max` 100/150 eV and with no `lam_min`, DL07 at 100 eV, Zubko): the served value against the table value on each side of 0.0912 um, plus the band-averaged `s_abs` and `C_ext(lambda_ref)` the transport is built on (`check_kext_table.f90`, `run_check.sh`).  Every model wavelength is now a table wavelength, so what is left is rounding and the EUV tolerance is 1e-9 |
