#!/usr/bin/env python3
"""Write a multi-column source-spectrum table for multiple_sources.in.

par%src_spectrum_file lets each point source draw its ionizing-band spectrum
from a column of one shared file (column i = source i), instead of a blackbody
temperature src_tstar(i).  This script writes two_spectra.txt with the same
Planck shape the code's planck_nu evaluates,
    B_E ~ E^3 / (exp(E * ev2erg / (kB * T)) - 1)     [arbitrary normalization],
on 400 log-spaced points over 10-110 eV (the ionizing band 13.598-100 eV sits
strictly inside).  The absolute scale is immaterial: MoCHII normalizes each
source's spectrum to its own src_lum(i).

Columns match multiple_sources.in: source 1 = 5e4 K, source 2 = 3.5e4 K.
Using this file (uncomment par%src_spectrum_file in the .in) reproduces the
src_tstar run to the file-interpolation level.
"""
import numpy as np

EV2ERG = 1.602176634e-12      # define.f90 ev2erg
KB_CGS = 1.380649e-16         # define.f90 kboltz_cgs


def planck_nu(E, T):
    """Code's planck_nu shape (arbitrary normalization), E [eV], T [K]."""
    x = E * EV2ERG / (KB_CGS * T)
    return E**3 / np.expm1(x)


E = np.logspace(np.log10(10.0), np.log10(110.0), 400)
cols = [planck_nu(E, 5.0e4), planck_nu(E, 3.5e4)]
with open('two_spectra.txt', 'w') as f:
    f.write('# tabulated Planck B_E (arbitrary normalization); shape = '
            'E^3/(exp(E*ev2erg/(kB*T))-1)\n')
    f.write('# col1 = E [eV]; col2 = B_E(5e4 K); col3 = B_E(3.5e4 K)\n')
    for i in range(E.size):
        f.write(f'{E[i]:.8e} {cols[0][i]:.8e} {cols[1][i]:.8e}\n')
print('wrote two_spectra.txt: %d points, %.1f-%.1f eV' % (E.size, E[0], E[-1]))
