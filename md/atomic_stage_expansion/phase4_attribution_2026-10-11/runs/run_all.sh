#!/bin/bash
# Phase 4 attribution runs on lart3, 16 ranks, one run at a time.
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
W=/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase5_attrib
S=$W/runs/status.txt
run(){ # set build-name binary
  mkdir -p $W/runs/$1/$2; cd $W/runs/$1/$2
  /bin/cp -p ../$2.in run.in
  echo "$1 $2 start $(date +%T)" >> $S
  mpirun -np 16 $W/$3/MoCHII.x run.in > run.log 2>&1
  echo "$1 $2 exit $? $(date +%T)" >> $S
}
run it2 d d; run it2 a a; run it2 b b; run it2 c d; run it2 noise d
run common d d
run rs d d; run rs a a; run rs b b; run rs c d; run rs noise d
