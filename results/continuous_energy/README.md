# Continuous-energy validation outputs

This directory holds generated HII20/HII40 output and is intentionally ignored by Git.

- `m2_atomic/`: the 2^25-packet Sobol production outputs the paper uses
  (2026-09-24, after the atomic-data update), with the paper's figure runs
  (`pdr_L5_atomic`, `hei_high_density_atomic`, `g2a_thermal_atomic`) and the
  sensitivity runs (`*_nohecx`, `*_lowden`, `*_ots`); `run_lart3*.sh` are the
  drivers.
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
