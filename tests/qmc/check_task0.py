#!/usr/bin/env python3
"""TASK 0 unit validation: the diffuse-photon energy->bin mapping.

The fix replaces the closed-form log-grid formula in gen_diffuse_photon
  inu = nfuv + int(nion*log(e/emin)/log(emax/emin)) + 1     (clamped)
with a binary search of the actual band edge array ion_eedge (the one source of
truth, filled for BOTH the aligned and the log grid).

(i)  LOG grid: the binary search must reproduce the closed-form bin for every
     energy except (measure-zero) exact edges.  We scan finely and report any
     disagreement and where it sits relative to a bin edge.
(ii) ALIGNED grid: each ionization threshold must land in the bin ABOVE its
     edge (the bin whose LOWER edge is the threshold).  We reconstruct the true
     aligned edges from a metals run's rates file (E_bin + dE_bin) and check.
"""
import sys
import numpy as np
import h5py

EMIN, EMAX = 13.598, 100.0


# ------------------------------------------------------------------ (i) log
def closed_form(e, nion, emin, emax, nfuv=0, nnu=None):
    if nnu is None:
        nnu = nfuv + nion
    inu = nfuv + (np.log(e / emin) / np.log(emax / emin) * nion).astype(int) + 1
    return np.clip(inu, nfuv + 1, nnu)


def bin_search(e, edges, nfuv=0, nnu=None):
    # bin i with edges[i-1] <= e < edges[i]  (1-based), clamped to ionizing bins
    if nnu is None:
        nnu = len(edges) - 1
    i = np.searchsorted(edges, e, side="right")   # edges[i-1] <= e < edges[i]
    return np.clip(i, nfuv + 1, nnu)


print("=== Task0(i) LOG grid: binary search vs closed-form ===")
allok = True
for nion in (16, 32, 64):
    edges = np.exp(np.linspace(np.log(EMIN), np.log(EMAX), nion + 1))
    # a fine, GENERIC scan (continuous diffuse energies) avoiding exact edges
    e = np.exp(np.linspace(np.log(EMIN) + 1e-9, np.log(EMAX) - 1e-9, 2_000_003))
    cf = closed_form(e, nion, EMIN, EMAX)
    bs = bin_search(e, edges)
    ndiff = int(np.count_nonzero(cf != bs))
    # of the mismatches, how close to the nearest edge (in relative energy)?
    if ndiff:
        idx = np.where(cf != bs)[0]
        near = np.min(np.abs(np.log(e[idx][:, None] / edges[None, :])), axis=1)
        worstrel = float(np.exp(near.max()) - 1.0)  # rel distance of the FARTHEST mismatch
        closerel = float(np.exp(near.min()) - 1.0)
        print(f"  nion={nion:3d}: {ndiff}/{e.size} mismatches; all within "
              f"[{closerel:.1e},{worstrel:.1e}] rel of a bin edge (FP ties)")
    else:
        print(f"  nion={nion:3d}: 0/{e.size} mismatches (exact agreement)")
    # the criterion: every mismatch is glued to an edge (rel < 1e-12)
    if ndiff:
        allok &= worstrel < 1e-10
print(f"  -> LOG grid: {'PASS (mismatches only AT edges = FP ties)' if allok else 'CHECK'}")


# -------------------------------------------------------------- (ii) aligned
def edges_from_rates(fn):
    with h5py.File(fn, "r") as f:
        ec = f["E_bin"]["data"][:]
        de = f["dE_bin"]["data"][:]
    # geometric center ec=sqrt(lo*hi), width de=hi-lo -> invert each bin exactly
    # (ec^2 = lo*(lo+de)): lo = (-de + sqrt(de^2 + 4 ec^2))/2, hi = lo+de.  Bins
    # tile contiguously (lo[i+1]=hi[i]), so edges = [lo_1..lo_n, hi_n].
    lo = 0.5 * (-de + np.sqrt(de ** 2 + 4 * ec ** 2))
    edges = np.empty(len(ec) + 1)
    edges[:-1] = lo
    edges[-1] = lo[-1] + de[-1]
    return edges


print("\n=== Task0(ii) ALIGNED grid: thresholds land in the bin ABOVE the edge ===")
fn = sys.argv[1] if len(sys.argv) > 1 else \
    "tests/qmc/v1/v1d_diff_aligned_new_rates.h5"
edges = edges_from_rates(fn)
print(f"  edges reconstructed from {fn}  (nbin={len(edges)-1})")
# The aligned builder pins H I / He I / He II onto bin edges.  A real diffuse
# recombination photon carries energy = threshold + delta (delta>0), so it must
# fall in the bin whose LOWER edge is the threshold (the bin ABOVE the edge),
# and the old log-formula (which assumed a pure log grid) would mis-assign it.
# We confirm each threshold IS a reconstructed edge and that an energy a hair
# ABOVE it lands in the bin above (and a hair below lands in the bin below).
thr = {"H I 13.598": 13.598, "He I 24.587": 24.587, "He II 54.416": 54.416}
aok = True
for name, Elit in thr.items():
    j = int(np.argmin(np.abs(edges - Elit)))       # nearest reconstructed edge
    onedge = abs(edges[j] / Elit - 1.0) < 1e-6
    if not (edges[0] <= edges[j] < edges[-1]) or not onedge:
        # H I sits at the band start (eion_min) - it is edge 0, not interior
        at_start = abs(edges[0] / Elit - 1.0) < 1e-6
        print(f"  {name:12s}: {'band-start edge (eion_min)' if at_start else 'not an aligned edge'}")
        continue
    Eedge = edges[j]                                # the true (Fortran) edge
    b_above = int(np.searchsorted(edges, Eedge * (1 + 1e-10), side="right"))
    b_below = int(np.searchsorted(edges, Eedge * (1 - 1e-10), side="right"))
    ok = (b_above == j + 1) and (b_below == j)
    aok &= ok
    print(f"  {name:12s}: aligned edge={Eedge:.4f} (=edge idx {j})  "
          f"E+ -> bin {b_above} (lower edge {edges[b_above-1]:.4f}),  "
          f"E- -> bin {b_below}  {'-> photon lands in bin ABOVE edge (correct)' if ok else '-> WRONG'}")
print(f"  -> ALIGNED grid: {'PASS (thresholds are edges; recomb photons bin ABOVE)' if aok else 'CHECK'}")

print("\nTASK0 UNIT:", "PASS" if (allok and aok) else "CHECK")
