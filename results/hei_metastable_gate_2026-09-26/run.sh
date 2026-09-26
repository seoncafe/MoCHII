#!/bin/bash
# check_hei_meta.py inputs (tests/hei_meta pins) on the current binary, and the
# band-stratified FUV gate, 2026-09-26.
cd /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/hei_metastable_gate_2026-09-26
X=/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/MoCHII.x
for b in hei_lowd hei_highd hei_fuv; do (mpirun -np 16 $X $b.in > $b.log 2>&1; echo "$b exit=$?" >> driver.log) & done
(WORKDIR=/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/fuv_stratified_gate_2026-09-26 \
   /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/tests/fuv_stratified/run_check.sh 16 \
   > /nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00/results/fuv_stratified_gate_2026-09-26/gate.out 2>&1; echo "fuv exit=$?" >> driver.log) &
wait
python3 check_hei_meta.py > gate.out 2>&1; echo "hei_meta gate exit=$?" >> driver.log
echo "ALL DONE" >> driver.log
