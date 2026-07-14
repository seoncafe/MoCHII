#!/bin/bash
# Solution-driven octree re-refinement of the ionization front
mpirun -np 8 ../../MoCHII.x adaptive_refinement.in
