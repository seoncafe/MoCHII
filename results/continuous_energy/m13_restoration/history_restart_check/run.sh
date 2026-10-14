cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/history_restart_check
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
X=/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x
mpirun -np 4 $X nofile.in > nofile.log 2>&1; echo nofile $?
mpirun -np 4 $X old.in > old.log 2>&1; echo old $?
mpirun -np 16 $X a.in > a.log 2>&1; echo a $?
mv a_rates_iterate.h5 a5_rates.h5
mpirun -np 16 /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x b.in > b.log 2>&1; echo b $?
