#!/bin/bash
# Neutral cloud in the standard interstellar radiation field (Draine ISRF
# preset): FUV penetrates and photoionizes C/Mg in a skin, dust stops it inside.
mpirun -np 8 ../../MoCHII.x isrf_cloud.in
