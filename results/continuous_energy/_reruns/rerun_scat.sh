#!/bin/bash
cd /nfs/mocafe/kiseon/RT_Codes/MoCHII/tests/d_dusty && mpirun -np 68 /nfs/mocafe/kiseon/RT_Codes/MoCHII/MoCHII.x d_dust_scat.in > d_dust_scat_new.log 2>&1
touch /nfs/mocafe/kiseon/RT_Codes/MoCHII/results/continuous_energy/_reruns/scat_done
