#!/bin/bash
cd /nfs/mocafe/kiseon/RT_Codes/MoCHII
R=results/continuous_energy/m1_s3
for b in hii40 hii20; do
  echo "=== ${b}_s3 start $(date +%H:%M:%S) ==="
  mpirun -np 68 ./MoCHII.x $R/${b}_s3.in > $R/${b}_s3.log 2>&1 < /dev/null
  echo "=== ${b}_s3 exit=$? $(date +%H:%M:%S) ==="
  for e in _rates.h5 _lines.txt _nebcont.txt; do [ -f "${b}_s3$e" ] && mv "${b}_s3$e" $R/; done
done
echo "S3 ALL DONE $(date +%H:%M:%S)"
