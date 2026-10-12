cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/refine_rates_check
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
mpirun -np 16 /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x new.in > new.log 2>&1 &
mpirun -np 16 /nfs/mocafe/kiseon/MoCHII/control_refine_20261012/MoCHII.x ctl.in > ctl.log 2>&1 &
wait; echo done
