#!/bin/bash
# Stromgren sphere with a physical-unit stellar spectrum (L_lambda vs lambda);
# the ionizing luminosity is derived from the file, not set in the namelist.
python3 make_physical_spectrum.py
mpirun -np 8 ../../MoCHII.x physical_spectrum.in
