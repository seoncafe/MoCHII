#!/usr/bin/env python3
"""Line drift of build d over the 3 restarted iterations: rs/d against
common2/d (L/Hbeta ratios - 1)."""
import numpy as np
R = '/nfs/mocafe/kiseon/MoCHII/MoCHII_v1.00_phase5_attrib/runs'
def lines(fn):
    L = {}; hb = None
    for l in open(fn):
        if l.startswith('# L(Hbeta) ='): hb = float(l.split('=')[1].split()[0])
        if l.startswith('#') or not l.strip(): continue
        p = l.split()
        try: L[(p[0], p[1], p[2])] = float(p[3])
        except Exception: pass
    return L, hb
a, ha = lines(f'{R}/common2/d/bb1e5_lines.txt'); b, hb = lines(f'{R}/rs/d/bb1e5_lines.txt')
print(f'L(Hbeta) {hb/ha-1:+.3e}')
for k, nm in [(('s','3','9532.25'),'[S III] 9532'), (('s','3','335008.38'),'[S III] 33.5um'), (('s','2','6718.29'),'[S II] 6716'),
              (('n','1','5201.71'),'[N I] 5200'), (('o','3','5008.24'),'[O III] 5007'), (('ne','4','2422.39'),'[Ne IV] 2422'), (('he','2','4687.00'),'He II 4686')]:
    print(f'{nm:16s} {(b[k]/hb)/(a[k]/ha)-1:+.3e}')
big = [(b[k]/hb)/(v/ha) - 1 for k, v in a.items() if k in b and v/ha > 1e-3]
print('median abs', np.median(np.abs(big)))
