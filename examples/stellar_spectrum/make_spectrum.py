#!/usr/bin/env python3
"""Convert a young stellar-population SED to the MoCHII ionizing-band format.

Source: spec_young_fsps.dat (FSPS SSP, MIST/MILES, solar Z, Kroupa IMF, 10 Myr;
from MoCafe v2.00 data/), columns lambda[um], L_lambda [shape only].
MoCHII's par%ion_spectrum wants (E[eV], L_E), so:
    E[eV]  = 1.23984 / lambda[um]       (= hc/lambda)
    L_E    ~ L_lambda * lambda**2        (L_E dE = L_lambda dlambda)
The absolute scale is arbitrary -- MoCHII renormalizes the shape so the band
[eion_min, eion_max] carries par%luminosity.
"""
import numpy as np
src = np.loadtxt('spec_young_fsps.dat')
lam_um, Llam = src[:, 0], src[:, 1]
E  = 1.23984 / lam_um
LE = Llam * lam_um**2
order = np.argsort(E);  E, LE = E[order], LE[order]
band = (E >= 10.0) & (E <= 200.0)          # cover 13.6-100 eV with a margin
np.savetxt('fsps_young_eV.dat', np.c_[E[band], LE[band]],
           header='FSPS SSP 10 Myr young population (from spec_young_fsps.dat)\n'
                  'E[eV]   L_E[arb per eV, shape only]', fmt='%.6e')
print('wrote fsps_young_eV.dat: %d points, %.1f-%.1f eV'
      % (band.sum(), E[band].min(), E[band].max()))
