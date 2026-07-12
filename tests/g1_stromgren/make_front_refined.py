#!/usr/bin/env python3
"""Front-refined octree sphere for the G1 AMR gate.

Same generic-AMR schema as MoCafe's make_amr_sphere.py (columns x, y, z,
level, nH, T, vx, vy, vz; keys BOXLEN, ORIGINX/Y/Z), but the refinement
targets the I-front: full level inside a shell around the expected
Stromgren radius, coarse elsewhere.  This is the refinement that matters
for the problem (docs/PLAN.md section 4) — the center-biased builder puts
its COARSEST cells at the front.

    python3 make_front_refined.py --out amr_sphere_front46.fits
"""
import argparse
import numpy as np


def build_leaves(half, lmin, lmax, r_front, w_front):
    cx, cy, cz, lv = [], [], [], []

    def target_level(x, y, z, h):
        # refine to lmax when the CELL (center +/- h) can overlap the front
        # shell [r_front - w_front, r_front + w_front]
        r = np.sqrt(x*x + y*y + z*z)
        rad = np.sqrt(3.0)*h
        if (r + rad >= r_front - w_front) and (r - rad <= r_front + w_front):
            return lmax
        return lmin

    def rec(x, y, z, h, level):
        if level < target_level(x, y, z, h):
            hc = 0.5*h
            for b in range(8):
                ix = b & 1; iy = (b >> 1) & 1; iz = (b >> 2) & 1
                rec(x + (2*ix - 1)*hc, y + (2*iy - 1)*hc, z + (2*iz - 1)*hc,
                    hc, level + 1)
        else:
            cx.append(x); cy.append(y); cz.append(z); lv.append(level)

    rec(0.0, 0.0, 0.0, half, 0)
    return (np.array(cx), np.array(cy), np.array(cz),
            np.array(lv, dtype=np.int32))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--boxlen", type=float, default=8.0)
    ap.add_argument("--level-min", type=int, default=4)
    ap.add_argument("--level-max", type=int, default=6)
    ap.add_argument("--n0", type=float, default=100.0)
    ap.add_argument("--r-front", type=float, default=3.04, help="front radius")
    ap.add_argument("--w-front", type=float, default=0.35, help="shell half width")
    args = ap.parse_args()

    half = 0.5*args.boxlen
    x, y, z, lv = build_leaves(half, args.level_min, args.level_max,
                               args.r_front, args.w_front)
    n = x.size
    r = np.sqrt(x*x + y*y + z*z)
    nH = np.where(r < half, args.n0, 0.0)
    T = np.full(n, 1.0e4)
    v0 = np.zeros(n)

    from astropy.io import fits
    cols = [fits.Column(name="x", format="1D", array=x),
            fits.Column(name="y", format="1D", array=y),
            fits.Column(name="z", format="1D", array=z),
            fits.Column(name="level", format="1J", array=lv.astype("i4")),
            fits.Column(name="nH", format="1D", array=nH),
            fits.Column(name="T", format="1D", array=T),
            fits.Column(name="vx", format="1D", array=v0),
            fits.Column(name="vy", format="1D", array=v0),
            fits.Column(name="vz", format="1D", array=v0)]
    tab = fits.BinTableHDU.from_columns(cols)
    for k, v in dict(BOXLEN=args.boxlen, ORIGINX=-half, ORIGINY=-half,
                     ORIGINZ=-half).items():
        tab.header[k] = v
    fits.HDUList([fits.PrimaryHDU(), tab]).writeto(args.out, overwrite=True)
    print(f"make_front_refined: wrote {n} leaves to {args.out} "
          f"(levels {lv.min()}-{lv.max()}, front shell "
          f"{args.r_front}+-{args.w_front})")


if __name__ == "__main__":
    main()
