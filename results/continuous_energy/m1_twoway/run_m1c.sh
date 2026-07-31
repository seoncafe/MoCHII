#!/bin/bash
cd /nfs/mocafe/kiseon/RT_Codes/MoCHII
R=results/continuous_energy/m1_twoway
while ! grep -q "HII20 ALL DONE" results/continuous_energy/metal_switches/driver20.log 2>/dev/null; do sleep 60; done
for b in hii40 hii20; do
  echo "=== ${b}_twoway start $(date +%H:%M:%S) ==="
  mpirun -np 68 ./MoCHII.x $R/${b}_twoway.in > $R/${b}_twoway.log 2>&1 < /dev/null
  echo "=== ${b}_twoway exit=$? $(date +%H:%M:%S) ==="
  for e in _rates.h5 _lines.txt _nebcont.txt; do [ -f "${b}_twoway$e" ] && mv "${b}_twoway$e" $R/; done
done
echo "M1C ALL DONE $(date +%H:%M:%S)"
