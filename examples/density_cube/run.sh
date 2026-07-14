#!/bin/bash
# 3D FITS density cube read directly onto a Cartesian grid
mpirun -np 8 ../../MoCHII.x density_cube.in
