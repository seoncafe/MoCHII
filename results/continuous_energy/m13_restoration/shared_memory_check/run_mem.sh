cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/m13_restoration/shared_memory_check
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
N=$1
( while true; do echo "$(date +%T) $(free -m | awk '/Mem:/{print $3}')"; sleep 15; done ) > mem_np$N.txt &
MON=$!
mpirun -hosts mocafe -n $N /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x pn150_mem.in > pn150_mem_np$N.log 2>&1
echo "exit $?" >> pn150_mem_np$N.log
kill $MON
