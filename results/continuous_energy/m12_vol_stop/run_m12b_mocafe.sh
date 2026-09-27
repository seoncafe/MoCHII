#!/bin/bash
# A third worker (mocafe, 88 ranks) on run_m12b.sh's queue.
R=/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00
cd $R/results/continuous_energy/m12_vol_stop
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
while [ -f queue_b.txt ]; do
  e=$(flock queue_b.lock bash -c 'n=$(head -1 queue_b.txt 2>/dev/null); [ -n "$n" ] && sed -i 1d queue_b.txt; echo $n')
  [ -z "$e" ] && break
  d=${e%%:*}; n=${e##*:}
  echo "${n} start mocafe $(date '+%F %T')" >> driver_b.log
  (cd $R/$d && mpirun -hosts mocafe -n 88 $R/MoCHII_m12run.x ${n}.in > ${n}.log 2>&1)
  echo "${n} exit=$? mocafe $(date '+%F %T')" >> driver_b.log
done
