#!/bin/bash
# He I 2^3S metastable diagnostic (par%hei_metastable): the 2^3S population
# n_3, the 10830 A absorption opacity, and the 2^3S photoionization rate Phi_3.
# The switch adds only the diagnostic output; the physics is unchanged.
#   low density   -> radiative limit (n_3/n_HeII ~ alpha_3 n_e/A_31)
#   high density  -> collisional plateau (n_e > n_crit suppresses n_3)
#   low + FUV     -> soft-UV Phi_3 depletes n_3 near the star
mpirun -np 8 ../../MoCHII.x hei_low_density.in
mpirun -np 8 ../../MoCHII.x hei_high_density.in
mpirun -np 8 ../../MoCHII.x hei_low_density_fuv.in
