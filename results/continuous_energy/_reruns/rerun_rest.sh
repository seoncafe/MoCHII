#!/bin/bash
cd /nfs/mocafe/kiseon/RT_Codes/MoCHII/tests/g2a_thermal && mpirun -np 68 /nfs/mocafe/kiseon/RT_Codes/MoCHII/MoCHII.x g2a_thermal.in > g2a_new.log 2>&1
cd /nfs/mocafe/kiseon/RT_Codes/MoCHII/tests/d_dusty && mpirun -np 68 /nfs/mocafe/kiseon/RT_Codes/MoCHII/MoCHII.x d_dust_full.in > d_dust_full_new.log 2>&1
cd /nfs/mocafe/kiseon/RT_Codes/MoCHII/tests/d_dusty && mpirun -np 68 /nfs/mocafe/kiseon/RT_Codes/MoCHII/MoCHII.x d_dust_l09.in  > d_dust_l09_new.log 2>&1
touch /nfs/mocafe/kiseon/RT_Codes/MoCHII/results/continuous_energy/_reruns/rest_done
