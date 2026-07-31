#!/bin/bash
cd /nfs/mocafe/kiseon/RT_Codes/MoCHII/tests/g2a_thermal && mpirun -np 68 /nfs/mocafe/kiseon/RT_Codes/MoCHII/MoCHII.x /nfs/mocafe/kiseon/RT_Codes/MoCHII/results/continuous_energy/_reruns/g2a_hg.in > /nfs/mocafe/kiseon/RT_Codes/MoCHII/results/continuous_energy/_reruns/g2a_hg.log 2>&1
touch /nfs/mocafe/kiseon/RT_Codes/MoCHII/results/continuous_energy/_reruns/hg_done
