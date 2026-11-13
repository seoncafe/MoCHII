#!/bin/bash
# Stromgren sphere: H/He photoionization equilibrium.
mpirun -np 8 ../../MoCHII.x stromgren_sphere.in
# the same sphere with the diffuse field re-emitted at the absorption point
mpirun -np 8 ../../MoCHII.x stromgren_reemission.in
