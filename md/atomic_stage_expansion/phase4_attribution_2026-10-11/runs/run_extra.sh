#!/bin/bash
# Second realization of the channel-6 pair: build a with qmc_seed 777, to
# set against the noise run (build d, qmc_seed 777).  Waits for run_all.sh.
export I_MPI_ADJUST_ALLREDUCE=1 I_MPI_ADJUST_REDUCE=1 I_MPI_ADJUST_BCAST=1 I_MPI_ADJUST_GATHER=1
W=/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase5_attrib
S=$W/runs/status.txt
until grep -q "rs noise exit" $S; do sleep 30; done
run(){ mkdir -p $W/runs/$1/$2; cd $W/runs/$1/$2; /bin/cp -p ../$2.in run.in
  echo "$1 $2 start $(date +%T)" >> $S
  mpirun -np 16 $W/$3/MoCHII.x run.in > run.log 2>&1
  echo "$1 $2 exit $? $(date +%T)" >> $S; }
run it2 a_noise a; run rs a_noise a
