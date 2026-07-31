#!/bin/bash
cd /nfs/mocafe/kiseon/RT_Codes/MoCHII/examples/hei_metastable && mpirun -np 68 /nfs/mocafe/kiseon/RT_Codes/MoCHII/MoCHII.x hei_high_density.in > hei_high_density_new.log 2>&1
touch /nfs/mocafe/kiseon/RT_Codes/MoCHII/results/continuous_energy/_reruns/hei_done
