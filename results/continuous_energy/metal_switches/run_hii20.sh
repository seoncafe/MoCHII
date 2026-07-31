#!/bin/bash
cd /nfs/mocafe/kiseon/RT_Codes/MoCHII
R=results/continuous_energy/metal_switches
# wait for the HII40 queue; the frozen binary makes this immune to an M1 rebuild
while ! grep -q "ALL DONE" $R/driver.log 2>/dev/null; do sleep 60; done
for tag in hii20_ne_off hii20_heat_off hii20_both_off; do
  echo "=== $tag start $(date +%H:%M:%S) ==="
  mpirun -np 68 ./MoCHII_preM1.x $R/$tag.in > $R/$tag.log 2>&1 < /dev/null
  echo "=== $tag exit=$? $(date +%H:%M:%S) ==="
done
echo "HII20 ALL DONE $(date +%H:%M:%S)"
