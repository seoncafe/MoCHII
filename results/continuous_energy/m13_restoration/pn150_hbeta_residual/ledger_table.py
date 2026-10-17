#!/usr/bin/env python3
"""Absorption ledgers (par%absorption_ledger, LEDGER: lines of the last
transport pass in each log) per stellar Q(H) = 5.416e47 s^-1, side by side:
for each source the photons emitted and those absorbed by H0, and the H0
total; the sum of the H0 absorptions of every source other than diffuse
channel 1 (the ground captures of H) is the case-B H ionization rate the
radiation supplies.  Usage: ledger_table.py <log> ..."""
import re, sys
Q = 5.416e47
ABS = ['H0', 'He0', 'He+', 'metals', 'dust']
runs = {}
for f in sys.argv[1:]:
    rows = {}
    for line in open(f):
        if line.startswith(' LEDGER:') and 'emitted' not in line:
            name = line[9:17].strip()
            w = re.sub(r'(\d)([+-]\d{3})\b', r'\1E\2', line[17:]).split()
            vals = [float(x) for x in w[:7]]
            rows[name] = dict(em=vals[0], esc=vals[6], **dict(zip(ABS, vals[1:6])))
    runs[f] = rows
names = []
for rows in runs.values():
    names += [n for n in rows if n not in names]
print(f"{'source':10s}" + ''.join(f"{f.replace('ledger_', '').replace('.log', ''):>26s}" for f in runs))
print(f"{'':10s}" + ''.join(f"{'emitted':>13s}{'to H0':>13s}" for f in runs))
for n in names:
    print(f"{n:10s}" + ''.join(
        f"{runs[f][n]['em']/Q:13.4f}{runs[f][n]['H0']/Q:13.4f}" if n in runs[f] else f"{'-':>13s}{'-':>13s}"
        for f in runs))
for lab, sel in (('H0 total', lambda n: True), ('H0 not dif1', lambda n: n != 'dif1')):
    print(f"{lab:10s}" + ''.join(
        f"{'':13s}{sum(r['H0'] for k, r in runs[f].items() if sel(k.split()[0]))/Q:13.4f}" for f in runs))
for a in ABS[1:4]:
    print(f"{a + ' tot':10s}" + ''.join(
        f"{'':13s}{sum(r[a] for r in runs[f].values())/Q:13.4f}" for f in runs))
