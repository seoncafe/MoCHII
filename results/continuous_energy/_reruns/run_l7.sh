#!/bin/bash
cd /nfs/mocafe/kiseon/RT_Codes/MoCHII/tests/g1_stromgren && mpirun -np 68 /nfs/mocafe/kiseon/RT_Codes/MoCHII/MoCHII.x g1_stromgren_l7.in > g1_l7.log 2>&1
touch /nfs/mocafe/kiseon/RT_Codes/MoCHII/results/continuous_energy/_reruns/l7_done
