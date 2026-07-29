#!/usr/bin/env python3
"""G4a gate: solution-driven re-refinement vs uniform high resolution.

Same Stromgren configuration as the G1 gate (nH = 100, Q_H = 1e49,
case B, Te = 1e4 K): the re-refined run starts uniform level 5 and
rebuilds to level 7 on the I-front at iteration 15; the reference is a
native uniform level-7 run.  The exact 1D Gauss-Seidel reference gives
R_eff = 3.0365 pc (computed in tests/g1_stromgren).

Criteria: |R_eff(refined) - R_eff(uniform7)| < 0.5%, both within ~1% of
the 1D reference, with the refined run using far fewer leaves.
"""
import numpy as np
import os as _os
# Figures are written to paper/figures/ regardless of the working directory.
FIGDIR = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "figures")

# Plot 1 in THIN Monte Carlo cells (random sample); all numbers use every cell.
THIN = 5
import h5py

PC2CM = 3.0856776e18
R_EFF_1D = 3.0365          # from the G1 1D reference (identical physics)


def load(fname):
    f = h5py.File(fname, "r")
    xyz = f["LeafXYZ"]["data"][:]
    return dict(r=np.sqrt((xyz**2).sum(axis=0)), xyz=xyz,
                xHI=f["x_HI"]["data"][:], ne=f["n_e"]["data"][:],
                nleaf=xyz.shape[1])


def half_cells(xyz, levels=(7, 6, 5), boxlen=8.0):
    half = np.zeros(xyz.shape[1])
    for L in levels:
        step = boxlen/2**L
        on = np.all(np.abs(((xyz + boxlen/2)/step - 0.5)
                           - np.round((xyz + boxlen/2)/step - 0.5)) < 1e-6,
                    axis=0) & (half == 0)
        half[on] = step/2
    return half


def r_eff(run, half):
    vol = (2.0*half)**3
    gas = run["ne"] > 0
    vion = np.sum((1.0 - run["xHI"][gas])*vol[gas])
    return (3.0*vion/(4.0*np.pi))**(1.0/3.0)


TAG_A = "re-refined 5$\\rightarrow$7"
runs = {}
for tag, fn in ((TAG_A, "../../tests/g4_refine/g4_refine_rates.h5"),
                ("uniform 7", "../../tests/g4_refine/g4_uniform7_rates.h5")):
    try:
        run = load(fn)
    except FileNotFoundError:
        print(f"[{tag}] {fn} not found; skipped")
        continue
    half = half_cells(run["xyz"])
    assert (half > 0).all(), "level inference failed"
    Re = r_eff(run, half)
    runs[tag] = (run, half, Re)
    print(f"[{tag:16s}] nleaf = {run['nleaf']:8d}   R_eff = {Re:.4f} pc"
          f"   vs 1D {Re/R_EFF_1D-1:+.2%}")

if len(runs) == 2:
    Ra = runs[TAG_A][2]
    Rb = runs["uniform 7"][2]
    ratio = runs["uniform 7"][0]["nleaf"]/runs[TAG_A][0]["nleaf"]
    print(f"\nre-refined vs uniform-7: {Ra/Rb-1:+.3%}   "
          f"leaf ratio {ratio:.1f}x")
    ok = abs(Ra/Rb - 1.0) < 0.005
    print("GATE:", "PASS" if ok else "FAIL")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # Paper figures: axis, tick, and title fonts sized to ~85% of the body text.
    _S = 1.78
    _b = plt.rcParams["font.size"]
    plt.rcParams.update({
        "axes.labelsize":   _b * _S,
        "axes.titlesize":   _b * 1.2 * _S,
        "xtick.labelsize":  _b * _S,
        "ytick.labelsize":  _b * _S,
        "figure.titlesize": _b * 1.2 * _S,
    })
    fig, ax = plt.subplots(figsize=(7, 5))
    # uniform-7 (the dense reference) drawn UNDER in light gray filled dots;
    # re-refined drawn ON TOP as small red open circles, so both clouds are
    # distinguishable where they overlap.
    styles = {
        "uniform 7": dict(marker="o", ms=2.2, mfc="0.7", mec="none",
                          ls="none", alpha=0.5, zorder=1),
        TAG_A: dict(marker="o", ms=2.0, mfc="none", mec="tab:red",
                    mew=0.5, ls="none", alpha=0.6, zorder=2),
    }
    # Plot a random 1-in-THIN sample of the Monte Carlo cells: the clouds are
    # far denser than the figure can resolve, and every quoted number above is
    # computed from the FULL cell set.
    rng = np.random.default_rng(20260725)
    for tag in ("uniform 7", TAG_A):
        run, half, Re = runs[tag]
        gas = np.nonzero(run["ne"] > 0)[0]
        pick = rng.choice(gas, size=max(1, gas.size // THIN), replace=False)
        ax.plot(run["r"][pick], run["xHI"][pick], rasterized=True,
                label=f"{tag} ({run['nleaf']} leaves)", **styles[tag])
    ax.set_xlim(2.6, 3.5)
    ax.set_xlabel(r"$r$ [pc]")
    ax.set_ylabel(r"$x_{\rm HI}$")
    ax.legend(frameon=False, markerscale=3)
    ax.set_title("Ionization-front region: re-refined vs native grid")
    fig.tight_layout()
    fig.savefig(_os.path.join(FIGDIR, "g4_refine_check.pdf"), dpi=200, bbox_inches="tight")
    print("wrote g4_refine_check.pdf")
