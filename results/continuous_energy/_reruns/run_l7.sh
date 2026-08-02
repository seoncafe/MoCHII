#!/bin/bash
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tests/g1_stromgren && mpirun -np 68 /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x g1_stromgren_l7.in > g1_l7.log 2>&1
touch /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/_reruns/l7_done
