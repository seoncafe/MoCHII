#!/usr/bin/env python3
"""Line-flux gate: MoCHII Tier-2 line ratios vs MOCASSIN lineFlux.out.

MoCHII '<case>_lines.txt' rows: elem stage lambda[A] L[erg/s] L/L(Hbeta).
MOCASSIN CEL rows: elem ion lower upper lambda[A] ratio counter.
Lines are matched by (element Z, stage, wavelength within 0.5%).

Both are ratios to Hbeta; MOCASSIN's Hbeta is its full transfer, ours the
case-B effective-recombination fit — a few-% systematic on all ratios.
Report the classic diagnostics and the full matched set.
"""
import sys
import numpy as np
import os

MBASE = os.path.expanduser(
    "~/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/tests/baselines")
ELEM_Z = {"c": 6, "n": 7, "o": 8, "ne": 10, "s": 16}

CASES = [("HII20", "hii20_dif_lines.txt", f"{MBASE}/01_gas_HII20/lineFlux.out"),
         ("HII40", "hii40_dif_lines.txt", f"{MBASE}/02_gas_HII40/lineFlux.out")]

#--- classic diagnostics to headline: (Z, stage, lambda [A], label)
KEY_LINES = [
    (8, 3, 5008.24, "[O III] 5007"),
    (8, 3, 4960.29, "[O III] 4959"),
    (8, 3, 4364.44, "[O III] 4363"),
    (8, 3, 883392.0, "[O III] 88um"),
    (8, 3, 518134.0, "[O III] 52um"),
    (8, 2, 3727.10, "[O II] 3726"),
    (8, 2, 3729.86, "[O II] 3729"),
    (7, 2, 6585.27, "[N II] 6584"),
    (7, 2, 5756.19, "[N II] 5755"),
    (7, 3, 573566.0, "[N III] 57um"),
    (10, 2, 128139.0, "[Ne II] 12.8um"),
    (10, 3, 155545.0, "[Ne III] 15.5um"),
    (10, 3, 3869.86, "[Ne III] 3869"),
    (16, 2, 6718.29, "[S II] 6716"),
    (16, 2, 6732.67, "[S II] 6731"),
    (16, 3, 187056.0, "[S III] 18.7um"),
    (16, 3, 9071.10, "[S III] 9069"),
    (6, 2, 1577371.0, "[C II] 158um"),
    (6, 3, 1908.73, "C III] 1909"),
]


def read_mochii(path):
    out = {}
    hbeta = None
    for ln in open(path):
        if ln.startswith("#"):
            if "L(Hbeta) =" in ln:
                hbeta = float(ln.split("=")[1].split()[0])
            continue
        toks = ln.split()
        if len(toks) != 5:
            continue
        el, st, wl, L, r = toks[0], int(toks[1]), float(toks[2]), \
            float(toks[3]), float(toks[4])
        out[(ELEM_Z[el], st, wl)] = r
    return out, hbeta


def read_mocassin(path):
    """CEL rows: elem ion lower upper lambda ratio counter."""
    out = {}
    hbeta = None
    for ln in open(path):
        if "Hbeta [E36" in ln:
            hbeta = float(ln.split(":")[1].split()[0])*1e36
        toks = ln.split()
        if len(toks) == 7:
            try:
                z, ion, l, u = (int(t) for t in toks[:4])
                wl, r = float(toks[4]), float(toks[5])
            except ValueError:
                continue
            if z in ELEM_Z.values() and wl > 100.0:
                out[(z, ion, wl)] = r
    return out, hbeta


def match(table, z, st, wl, tol=5e-3):
    best, bdev = None, tol
    for (tz, ts, twl), r in table.items():
        if tz == z and ts == st:
            dev = abs(twl/wl - 1.0)
            if dev < bdev:
                best, bdev = r, dev
    return best


for label, ours_f, ref_f in CASES:
    if not os.path.exists(ours_f):
        print(f"[{label}] {ours_f} not found; skipped\n")
        continue
    ours, hb_us = read_mochii(ours_f)
    ref, hb_moc = read_mocassin(ref_f)
    print(f"=== {label} ===  L(Hbeta): MoCHII {hb_us:.3e}, "
          f"MOCASSIN {hb_moc:.3e} erg/s ({hb_us/hb_moc-1:+.1%})")
    print(f"{'line':16s} {'MoCHII':>10s} {'MOCASSIN':>10s} {'ratio':>8s}")
    for z, st, wl, name in KEY_LINES:
        r_us = match(ours, z, st, wl)
        r_mc = match(ref, z, st, wl)
        if r_us is None and r_mc is None:
            continue
        su = f"{r_us:10.4f}" if r_us is not None else "         -"
        sm = f"{r_mc:10.4f}" if r_mc is not None else "         -"
        rr = f"{r_us/r_mc:8.2f}" if (r_us and r_mc) else "       -"
        print(f"{name:16s} {su} {sm} {rr}")
    print()
