#!/bin/bash
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/examples/hei_metastable && mpirun -np 68 /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x hei_high_density.in > hei_high_density_new.log 2>&1
touch /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/continuous_energy/_reruns/hei_done
