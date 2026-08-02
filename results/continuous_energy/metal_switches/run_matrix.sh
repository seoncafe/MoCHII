#!/bin/bash
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00
R=results/continuous_energy/metal_switches
for tag in hii40_ne_off hii40_heat_off hii40_both_off; do
  echo "=== $tag start $(date +%H:%M:%S) ==="
  mpirun -np 68 ./MoCHII.x $R/$tag.in > $R/$tag.log 2>&1 < /dev/null
  echo "=== $tag exit=$? $(date +%H:%M:%S) ==="
done
echo "ALL DONE $(date +%H:%M:%S)"
