#!/bin/bash
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tests/g2a_thermal && mpirun -np 68 /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/_reruns/g2a_hg.in > /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/_reruns/g2a_hg.log 2>&1
touch /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/_reruns/hg_done
