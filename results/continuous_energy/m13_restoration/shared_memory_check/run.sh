cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/shared_memory_check
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
X=/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x
mpirun -np 16 $X a2.in > a2.log 2>&1 &
mpirun -np 16 $X refine.in > refine.log 2>&1 &
mpirun -np 16 $X emis.in > emis.log 2>&1 &
wait; echo done
