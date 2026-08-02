#!/bin/bash
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tests/d_dusty && mpirun -np 68 /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x d_dust_scat.in > d_dust_scat_new.log 2>&1
touch /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/_reruns/scat_done
