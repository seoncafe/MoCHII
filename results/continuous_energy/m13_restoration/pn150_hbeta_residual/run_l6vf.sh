#!/bin/bash
# PN150 level 6 with par%metal_recombination_photons = 'valence_fractions'
# (2026-10-17): from the converged l6c_metals state to convergence, to set
# R_B(H)/Q against l6c_metals (1.0068) and Cloudy (1.1184).
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/pn150_hbeta_residual
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
mpirun -hosts lart1 -n 36 ./MoCHII_vf.x l6vf_metals.in > l6vf_metals.log 2>&1
