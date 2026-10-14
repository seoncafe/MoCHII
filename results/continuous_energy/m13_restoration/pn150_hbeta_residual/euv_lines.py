#!/usr/bin/env python3
"""Photons per second in the MoCHII PN150 lines shortward of the H edge
(911.75 A, vacuum): line photons that can ionize H.  From
../pn150_p27_nodrift_lines.txt (L [erg/s] by element, stage, lambda)."""
import collections
h_c = 6.62607015e-27*2.99792458e10
Q = 5.416e47
tot = 0.0; by = collections.Counter(); top = []
for r in open('../pn150_p27_nodrift_lines.txt'):
    if r.startswith('#'): continue
    f = r.split(); el, st, wl, L = f[0], int(f[1]), float(f[2]), float(f[3])
    if el in ('h', 'he', 'hc') or wl <= 0 or wl >= 911.75: continue
    n = L*wl*1e-8/h_c
    tot += n; by[(el, st)] += n; top.append((n, el, st, wl))
print(f'metal lines below 911.75 A: {tot:.4e} photons/s = {tot/Q:.4f} Q(H)')
for (el, st), n in by.most_common(10): print(f'  {el:3s} {st}: {n/Q:.4f} Q')
for n, el, st, wl in sorted(top, reverse=True)[:10]: print(f'  {el:3s} {st} {wl:9.2f} A  {n/Q:.4f} Q')
