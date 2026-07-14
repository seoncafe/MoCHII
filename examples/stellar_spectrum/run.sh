#!/bin/bash
# Photoionization by a young stellar-population SED (par%ion_spectrum).
python3 make_spectrum.py            # spec_young_fsps.dat -> fsps_young_eV.dat
mpirun -np 8 ../../MoCHII.x stellar_spectrum.in
