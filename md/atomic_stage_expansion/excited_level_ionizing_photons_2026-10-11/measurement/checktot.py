import gzip, re, sys, numpy as np
import adf, recomb_photons as R
def file_total(path):
    L = gzip.open(path, 'rt', errors='replace').read().splitlines()
    k = max(i for i, l in enumerate(L) if re.match(r'^\s*----\s+--------', l))
    T = []; A = []
    for l in L[k+1:]:
        t = l.split()
        if len(t) < 2: break
        try: T.append(adf.fnum(t[0])); A.append(adf.fnum(t[1]))
        except Exception: break
    return np.array(T), np.array(A)
for f in sys.argv[1:]:
    d = adf.read(f); T, A = file_total(f)
    r = R.nl_rates(d, T)
    s = sum(r.values())
    print(f.split('/')[-1], ' '.join(f'{x:.3f}' for x in (s / A)[4:12]))
