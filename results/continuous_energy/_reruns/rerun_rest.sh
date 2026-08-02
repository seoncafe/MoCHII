#!/bin/bash
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tests/g2a_thermal && mpirun -np 68 /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x g2a_thermal.in > g2a_new.log 2>&1
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tests/d_dusty && mpirun -np 68 /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x d_dust_full.in > d_dust_full_new.log 2>&1
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tests/d_dusty && mpirun -np 68 /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x d_dust_l09.in  > d_dust_l09_new.log 2>&1
touch /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/_reruns/rest_done
