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
for tag, fn in ((TAG_A, "g4_refine_rates.h5"),
                ("uniform 7", "g4_uniform7_rates.h5")):
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
    fig, ax = plt.subplots(figsize=(7, 5))
    for (tag, (run, half, Re)), c in zip(runs.items(), ("C2", "C0")):
        gas = run["ne"] > 0
        ax.plot(run["r"][gas], run["xHI"][gas], ".", ms=1.2, color=c,
                alpha=0.3, label=f"{tag} ({run['nleaf']} leaves)")
    ax.set_xlim(2.6, 3.5)
    ax.set_xlabel(r"$r$ [pc]")
    ax.set_ylabel(r"$x_{\rm HI}$")
    ax.legend(frameon=False, markerscale=8)
    ax.set_title("G4 gate: I-front region")
    fig.tight_layout()
    fig.savefig("g4_refine_check.png", dpi=140)
    print("wrote g4_refine_check.png")
