cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/rates_iterate_check
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
mpirun -np 16 /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x a2.in > a2.log 2>&1; echo a2 $?
mpirun -np 16 /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x b2.in > b2.log 2>&1; echo b2 $?
