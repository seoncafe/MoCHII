#!/usr/bin/env python3
"""RQMC variance reduction: Owen-scrambled Sobol launch vs ordinary Monte Carlo.

Reproduces the V3 variance study (tests/qmc/check_v3.py) and plots it for the
paper.  On the fixed-state G0 attenuation test the direct field is an analytic
function of the launch variables, so a low-discrepancy (Owen-scrambled Sobol)
launch set reduces the noise of Gamma_HI directly.  For each photon count N and
each launch mode this script:

  * reads the replicate rates files (four MC replicates with different iseed,
    four RQMC replicates with different qmc_seed);
  * forms the cell-wise relative error of Gamma_HI against the analytic
    cell-averaged reference (the same reference check_v3.py / the g0 checker
    use), and its root-mean-square over the shell 0.15 < r < 0.85;
  * takes the RMS across replicates.

The figure plots RMS error vs N on log-log for MC and RQMC, with reference
slopes N^-1/2 (Monte Carlo) and N^-1 (RQMC), and annotates the effective
photon gain (the variance ratio) at each N.

Output: figures/qmc_variance.pdf   (vector)

Data (regenerable, NOT stored in git -- run tests/qmc/run_v3.sh once):
  ../../tests/qmc/v3/{mc,qmc}_{p17,p20}_*_rates.h5

Run:  python3 make_qmc.py
"""

import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"

# Axis, tick, and title fonts sized to ~85% of the body text once the single-
# column figure is scaled down by \includegraphics in the two-column paper.
_S = 1.75
_b = plt.rcParams["font.size"]
plt.rcParams.update({
    "axes.labelsize":   _b * _S,
    "axes.titlesize":   _b * 1.2 * _S,
    "xtick.labelsize":  _b * _S,
    "ytick.labelsize":  _b * _S,
    "figure.titlesize": _b * 1.2 * _S,
})

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "..", "figures")
V3 = os.path.join(HERE, "..", "..", "tests", "qmc", "v3")
OUT = os.path.join(FIGDIR, "qmc_variance.pdf")

# --- geometry and band of the V3 fixed-state test (from check_v3.py) ---------
EV2ERG = 1.602176634e-12
KB = 1.380649e-16
PC2CM = 3.0856776e18
ETH_HI, ETH_HeI = 13.598, 24.587
NNU, EMIN, EMAX = 16, 13.598, 100.0
TSTAR, LBAND = 1.0e5, 1.0e38
NH, XHI, XHEI, YHE = 1.0, 0.1, 1.0, 0.1
RSPH_CM = 1.0 * PC2CM
HALF_PC = 1.0 / 32.0

# photon-count tags produced by run_v3.sh (tag -> N)
PTAGS = [("p17", 2 ** 17), ("p18", 2 ** 18), ("p19", 2 ** 19), ("p20", 2 ** 20)]


def sigma_vfky96(E, Eth, E0, s0, ya, P, yw, y0, y1):
    E = np.asarray(E, float)
    x = E / E0 - y0
    z = np.sqrt(x * x + y1 * y1)
    Q = 5.5 - 0.5 * P
    s = s0 * ((x - 1.0) ** 2 + yw ** 2) * z ** (-Q) \
        * (1.0 + np.sqrt(z / ya)) ** (-P) * 1e-18
    return np.where(E >= Eth, s, 0.0)


def sig_HI(E):
    return sigma_vfky96(E, ETH_HI, 4.298e-1, 5.475e4, 3.288e1, 2.963, 0, 0, 0)


def sig_HeI(E):
    return sigma_vfky96(E, ETH_HeI, 1.361e1, 9.492e2, 1.469, 3.188, 2.039,
                        4.434e-1, 2.136)


def band_bins(nnu, emin, emax, nsub=32):
    edge = np.exp(np.linspace(np.log(emin), np.log(emax), nnu + 1))
    ec = np.sqrt(edge[:-1] * edge[1:])
    lum = np.zeros(nnu)
    for i in range(nnu):
        es = edge[i] + (edge[i + 1] - edge[i]) * (np.arange(nsub) + 0.5) / nsub
        x = es * EV2ERG / (KB * TSTAR)
        lum[i] = np.sum(es ** 3 / np.expm1(x)) * (edge[i + 1] - edge[i]) / nsub
    lum *= LBAND / lum.sum()
    return ec, lum


def gamma_ana(r_cm, ec, lum):
    kap = NH * (XHI * sig_HI(ec) + YHE * XHEI * sig_HeI(ec))
    tau = np.outer(np.minimum(r_cm, RSPH_CM), kap)
    return np.sum(lum * sig_HI(ec) / (ec * EV2ERG) * np.exp(-tau), axis=1) \
        / (4 * np.pi * r_cm ** 2)


def gamma_ana_cellavg(xyz_pc, half_pc, ec, lum, nsub=4):
    off = (np.arange(nsub) + 0.5) / nsub - 0.5
    ox, oy, oz = np.meshgrid(off, off, off, indexing="ij")
    pts = np.stack([o.ravel() for o in (ox, oy, oz)])
    acc = np.zeros(xyz_pc.shape[1])
    for k in range(pts.shape[1]):
        p = xyz_pc + 2.0 * half_pc * pts[:, k:k + 1]
        r = np.sqrt((p ** 2).sum(axis=0)) * PC2CM
        acc += gamma_ana(r, ec, lum)
    return acc / pts.shape[1]


def build_reference():
    """Analytic cell-averaged Gamma_HI on the run geometry (needs one h5)."""
    import h5py
    any_h5 = sorted(glob.glob(os.path.join(V3, "*_rates.h5")))
    if not any_h5:
        return None
    with h5py.File(any_h5[0], "r") as f0:
        xyz = f0["LeafXYZ"]["data"][:]
        e_bin = f0["E_bin"]["data"][:]
    r_pc = np.sqrt((xyz ** 2).sum(axis=0))
    ec, lum = band_bins(NNU, EMIN, EMAX)
    if not np.allclose(ec, e_bin, rtol=1e-9):
        raise SystemExit("bins disagree with the analytic log grid "
                         "(is ion_align_edges off in run_v3.sh?)")
    ana = gamma_ana_cellavg(xyz, HALF_PC, ec, lum)
    sel = (r_pc > 0.15) & (r_pc < 0.85)
    return ana, sel


def rms_cell(fn, ana, sel):
    import h5py
    with h5py.File(fn, "r") as f:
        g = f["Gamma_HI"]["data"][:]
    rel = g[sel] / ana[sel] - 1.0
    return np.sqrt(np.mean(rel ** 2))


def collect():
    """{mode: {N: mean RMS over replicates}} from the V3 outputs, or None."""
    ref = build_reference()
    if ref is None:
        return None
    ana, sel = ref
    res = {"mc": {}, "qmc": {}}
    for mode in ("mc", "qmc"):
        for tag, N in PTAGS:
            fns = sorted(glob.glob(os.path.join(V3, f"{mode}_{tag}_*_rates.h5")))
            if not fns:
                continue
            rc = np.array([rms_cell(fn, ana, sel) for fn in fns])
            # RMS across replicates (the check_v3.py aggregate)
            res[mode][N] = float(np.sqrt(np.mean(rc ** 2)))
    return res


def main():
    res = collect()
    if not res or not res["mc"] or not res["qmc"]:
        print("MISSING DATA: no V3 rates files under", V3)
        print("  The tests/qmc/ tree stores only inputs, logs and checkers;")
        print("  the *_rates.h5 outputs are regenerable but not in git.")
        print("  Regenerate them once (single point source, fixed-state G0):")
        print("      bash tests/qmc/run_v3.sh")
        print("  then re-run:  python3 make_qmc.py")
        sys.exit(1)

    Ns = np.array([N for _, N in PTAGS if N in res["mc"] and N in res["qmc"]])
    mc = np.array([res["mc"][N] for N in Ns])
    qmc = np.array([res["qmc"][N] for N in Ns])
    gain = (mc / qmc) ** 2                    # variance ratio = eff. photon gain

    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    ax.loglog(Ns, 100 * mc, "o-", color="#1f77b4", lw=2.0, ms=8,
              label=r'Monte Carlo (\texttt{random})')
    ax.loglog(Ns, 100 * qmc, "s-", color="#d62728", lw=2.0, ms=8,
              label=r'RQMC (Owen--Sobol)')

    # reference slopes anchored at the smallest-N measured point of each method
    Nr = np.array([Ns[0], Ns[-1]], float)
    ax.loglog(Nr, 100 * mc[0] * (Nr / Ns[0]) ** -0.5, ":", color="#1f77b4",
              lw=1.4, label=r'$N^{-1/2}$')
    ax.loglog(Nr, 100 * qmc[0] * (Nr / Ns[0]) ** -1.0, ":", color="#d62728",
              lw=1.4, label=r'$N^{-1}$')

    # annotate the effective photon gain at each N
    for N, m, q, g in zip(Ns, mc, qmc, gain):
        ax.annotate(r'$%.0f\times$' % g,
                    xy=(N, 100 * np.sqrt(m * q)),
                    ha="center", va="center", fontsize=_b * _S * 0.95,
                    color="black",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white",
                              ec="0.6", lw=0.6))

    ax.set_xlabel(r'launched photons $N$')
    ax.set_ylabel(r'cell-wise RMS error in $\Gamma_{\rm HI}$ [\%]')
    ax.set_title(r'Quasi-random launch: variance reduction')
    ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.5)
    ax.legend(loc="lower left", frameon=False, fontsize=_b * _S * 0.85)

    fig.tight_layout()
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)

    print("wrote", OUT)
    for N, m, q, g in zip(Ns, mc, qmc, gain):
        print("  N = 2^%d: RMS  MC = %6.3f%%  RQMC = %6.3f%%  "
              "eff. photon gain = %6.1fx" % (int(np.log2(N)), 100 * m,
                                             100 * q, g))


if __name__ == "__main__":
    main()
