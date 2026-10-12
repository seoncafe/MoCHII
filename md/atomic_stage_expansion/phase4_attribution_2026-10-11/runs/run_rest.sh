#!/bin/bash
# Continuation after run_all.sh was stopped at the end of common/d: 12 more
# iterations of build d (common2), then the restart set rs and the second
# channel-6 realization (build a at qmc_seed 777).  lart3, 16 ranks.
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
W=/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase5_attrib
S=$W/runs/status.txt
while kill -0 46297 2>/dev/null; do sleep 20; done
echo "common d finished (run_all.sh stopped) $(date +%T)" >> $S
run(){ mkdir -p $W/runs/$1/$2; cd $W/runs/$1/$2; /bin/cp -p ../$2.in run.in
  echo "$1 $2 start $(date +%T)" >> $S
  mpirun -np 16 $W/$3/MoCHII.x run.in > run.log 2>&1
  echo "$1 $2 exit $? $(date +%T)" >> $S; }
run common2 d d
run rs d d; run rs a a; run rs noise d; run rs a_noise a; run rs b b; run rs c d
run it2 a_noise a
