# Continuous-energy validation outputs

This directory holds generated HII20/HII40 output and is intentionally ignored by Git.

- `production/`: the 2^25-packet Sobol production outputs used by the paper.
- `bin_independence/`: diagnostic-bin invariance runs.
- `rqmc_ensemble/`: independent scrambled-Sobol production ensemble.
- `mpi_scaling/` and `mpi_multinode/`: MPI agreement and timing runs.
- `smoke/` and `misc/`: short diagnostic outputs and logs.

Run scripts under `tests/continuous_energy/` write to the matching directory.
