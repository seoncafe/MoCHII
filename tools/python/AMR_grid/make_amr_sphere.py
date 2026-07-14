#!/usr/bin/env python3
"""Build a uniform-density nH sphere (or shell) octree for MoCHII.

Writes the generic AMR file MoCHII reads via ``par%amr_file`` with
``grid_type='amr'`` (src/read_generic_amr.f90): columns x, y, z, level, nH
[, xHI] and header keywords BOXLEN, ORIGINX/Y/Z, NLEAF.

The medium is a uniform-density sphere (or shell): every leaf whose center
lies in rmin <= r < rmax gets nH = n0, the rest get nH = 0.  The octree can be
uniform (every leaf at --level) or radially refined (finer toward the center),
which exercises the octree's level-crossing ray traversal while the physical
medium stays a uniform sphere.

Uses the AMRGrid octree builder from AMR_grid.py.

Examples
--------
    python make_amr_sphere.py --level 5 --n0 100 --out sphere_uniform.fits
    python make_amr_sphere.py --level-min 3 --level-max 6 --n0 100 --rmax 0.8 \
        --out sphere_refined.fits
    python make_amr_sphere.py --level 5 --n0 100 --rmin 0.4 --rmax 0.8 \
        --xhi 1.0 --out shell.h5
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from AMR_grid import AMRGrid   # noqa: E402


def _target_level(r, half, lmin, lmax):
    """Refinement level as a function of center radius: lmax toward the center,
    ramping down to lmin at the box corner (finer toward the center)."""
    if lmax == lmin:
        return lmin
    t = np.minimum(1.0, r / (0.75 * half))
    return int(round(lmax - t * (lmax - lmin)))


def build_grid(boxlen, lmin, lmax, n0, rmin, rmax, xhi=None):
    half = 0.5 * boxlen
    grid = AMRGrid(boxlen=boxlen)          # origin centered at 0, cubic box

    def radius(c):
        return np.sqrt(c.cx*c.cx + c.cy*c.cy + c.cz*c.cz)

    if lmax == lmin:
        grid.refine_uniform(lmin)
    else:
        grid.refine_uniform(lmin)
        grid.refine(lambda c: c.level < _target_level(radius(c), half, lmin, lmax),
                    lmax)

    def density(x, y, z):
        r = np.sqrt(x*x + y*y + z*z)
        return np.where((r >= rmin) & (r < rmax), n0, 0.0)

    grid.set_density(density)
    if xhi is not None:
        grid.set_neutral_fraction(xhi)     # uniform initial neutral fraction
    return grid


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build a uniform-density nH sphere/shell octree for MoCHII.")
    ap.add_argument("--out", required=True,
                    help="output file (.fits / .fits.gz / .h5)")
    ap.add_argument("--boxlen", type=float, default=2.0,
                    help="box side length [code units] (default 2.0)")
    ap.add_argument("--level", type=int, default=None,
                    help="uniform refinement level (overrides --level-min/max)")
    ap.add_argument("--level-min", type=int, default=3,
                    help="coarsest level for radial refinement (default 3)")
    ap.add_argument("--level-max", type=int, default=6,
                    help="finest level at the center for radial refinement "
                         "(default 6)")
    ap.add_argument("--n0", type=float, default=1.0,
                    help="nH inside the sphere/shell [cm^-3] (default 1.0)")
    ap.add_argument("--rmax", type=float, default=None,
                    help="outer radius [code units] (default = box half-size)")
    ap.add_argument("--rmin", type=float, default=0.0,
                    help="inner radius [code units]; > 0 makes a shell "
                         "(default 0)")
    ap.add_argument("--xhi", type=float, default=None,
                    help="optional uniform initial neutral fraction column")
    args = ap.parse_args(argv)

    half = 0.5 * args.boxlen
    rmax = args.rmax if args.rmax is not None else half
    if args.level is not None:
        lmin = lmax = args.level
    else:
        lmin, lmax = args.level_min, args.level_max

    grid = build_grid(args.boxlen, lmin, lmax, args.n0, args.rmin, rmax,
                      xhi=args.xhi)
    grid.write(args.out)
    print(grid.info())
    print(f"make_amr_sphere: rmin={args.rmin:g}, rmax={rmax:g}, "
          f"n0={args.n0:g}"
          + (f", xHI={args.xhi:g}" if args.xhi is not None else ""))


if __name__ == "__main__":
    main()
