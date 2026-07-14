#!/bin/bash
# Dusty photoionized sphere with EUV dust absorption, scattering, and IR emission
mpirun -np 8 ../../MoCHII.x dusty_nebula.in
