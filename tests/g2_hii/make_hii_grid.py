#!/usr/bin/env python3
"""Hollow uniform sphere for the MOCASSIN HII20/HII40 gate (G2c).

Same generic-AMR schema as the other builders: nH = n0 between r_in and
r_out (the MOCASSIN benchmarks start the gas at Rin = 3e18 cm), zero
elsewhere.  Uniform octree level.

    python3 make_hii_grid.py --rin 0.97223 --rout 3.07874 --out hii20.fits
    python3 make_hii_grid.py --rin 0.97223 --rout 4.73154 --out hii40.fits
"""
import argparse
import numpy as np


def build_uniform(half, level):
    n = 2**level
    step = 2*half/n
    c = -half + (np.arange(n) + 0.5)*step
    x, y, z = np.meshgrid(c, c, c, indexing="ij")
    lv = np.full(x.size, level, dtype=np.int32)
    return x.ravel(), y.ravel(), z.ravel(), lv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--rin", type=float, required=True, help="inner radius [pc]")
    ap.add_argument("--rout", type=float, required=True, help="outer radius [pc]")
    ap.add_argument("--level", type=int, default=6)
    ap.add_argument("--n0", type=float, default=100.0)
    args = ap.parse_args()

    half = args.rout
    x, y, z, lv = build_uniform(half, args.level)
    r = np.sqrt(x*x + y*y + z*z)
    nH = np.where((r >= args.rin) & (r < args.rout), args.n0, 0.0)
    n = x.size
    T = np.full(n, 1.0e4)
    v0 = np.zeros(n)

    from astropy.io import fits
    cols = [fits.Column(name="x", format="1D", array=x),
            fits.Column(name="y", format="1D", array=y),
            fits.Column(name="z", format="1D", array=z),
            fits.Column(name="level", format="1J", array=lv),
            fits.Column(name="nH", format="1D", array=nH),
            fits.Column(name="T", format="1D", array=T),
            fits.Column(name="vx", format="1D", array=v0),
            fits.Column(name="vy", format="1D", array=v0),
            fits.Column(name="vz", format="1D", array=v0)]
    tab = fits.BinTableHDU.from_columns(cols)
    for k, v in dict(BOXLEN=2*half, ORIGINX=-half, ORIGINY=-half,
                     ORIGINZ=-half).items():
        tab.header[k] = v
    fits.HDUList([fits.PrimaryHDU(), tab]).writeto(args.out, overwrite=True)
    print(f"make_hii_grid: wrote {n} leaves to {args.out} "
          f"(boxlen {2*half:.4f} pc, shell {args.rin}-{args.rout} pc, "
          f"f_gas={np.mean(nH>0):.3f})")


if __name__ == "__main__":
    main()
