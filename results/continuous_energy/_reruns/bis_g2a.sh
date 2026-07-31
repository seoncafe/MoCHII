#!/bin/bash
cd /nfs/mocafe/kiseon/RT_Codes/MoCHII_bisect/tests/g2a_thermal && mpirun -np 68 /nfs/mocafe/kiseon/RT_Codes/MoCHII_bisect/MoCHII.x g2a_thermal.in > bis_g2a.log 2>&1
touch /nfs/mocafe/kiseon/RT_Codes/MoCHII/results/continuous_energy/_reruns/bis_done
