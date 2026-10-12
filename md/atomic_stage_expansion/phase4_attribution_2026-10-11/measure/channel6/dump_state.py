#!/usr/bin/env python3
"""Write the leaf state of one run as a stream file for channel6_spectrum.x:
nleaf (int32), then float64 arrays T_e, n_e, n_H, V[cm^3], shell index,
and x stages (7 columns) of c, n, o, ne, s.  dump_state.py <run dir> <out>"""
import sys, h5py, numpy as np
d = sys.argv[1]
r = h5py.File(f'{d}/bb1e5_rates.h5', 'r'); e = h5py.File(f'{d}/bb1e5_emis.h5', 'r')
xyz = r['LeafXYZ/data'][:]; rad = np.sqrt((xyz**2).sum(0))
V = (r['LeafSize/data'][:]*3.0856775814913673e18)**3
shell = np.digitize(rad, [2.9, 4.0, 5.3]).astype(np.float64) + 1.0    # 1..4
with open(sys.argv[2], 'wb') as f:
    np.array([len(rad)], dtype=np.int32).tofile(f)
    for a in [e['T_e/data'][:], e['n_e/data'][:], e['n_H/data'][:], V, shell]:
        np.ascontiguousarray(a, dtype=np.float64).tofile(f)
    for el in ['c', 'n', 'o', 'ne', 's']:
        np.ascontiguousarray(e[f'x_{el}_stages/data'][:].T, dtype=np.float64).tofile(f)
