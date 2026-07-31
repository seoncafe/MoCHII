#!/bin/bash
cd /nfs/mocafe/kiseon/RT_Codes/MoCHII
R=results/continuous_energy/m1_oxfix
for b in hii40 hii20; do
  echo "=== ${b} start $(date +%H:%M:%S) ==="
  mpirun -np 68 ./MoCHII.x $R/${b}_oxfix.in > $R/${b}_oxfix.log 2>&1 < /dev/null
  echo "=== ${b} exit=$? $(date +%H:%M:%S) ==="
  for e in _rates.h5 _lines.txt _nebcont.txt; do [ -f "${b}_oxfix$e" ] && mv "${b}_oxfix$e" $R/; done
done
echo "OXFIX ALL DONE $(date +%H:%M:%S)"
