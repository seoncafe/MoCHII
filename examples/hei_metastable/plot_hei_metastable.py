#!/usr/bin/env python3
"""Figures for the He I 2^3S metastable diagnostic example.

Reads the three runs' <base>_heimeta.h5 (n_2s3, Phi_3, n_e, T_e) and
<base>_rates.h5 (x_HeII, x_HI) through tools/python/mochii_output.py and makes:

  hei_profiles.png  radial n_3/n_HeII: low-density run in the radiative limit
                    alpha_3 n_e/A_31, high-density run pushed below it toward the
                    collisional plateau alpha_3/(q_31a+q_31b).
  hei_maps.png      high-density run: the 2^3S column density and the 10830 A
                    doublet line-center optical depth (tau_0 >> 1 in the core --
                    the optically-thin warning is physically correct).
  hei_phi3.png      soft-UV depletion: n_3^fuv/n_3^noFUV and Phi_3/A_31 vs radius.

The metastable balance and the coefficients are from OH18 (Oklopcic & Hirata
2018) as transcribed in data/atomic/hei_metastable.txt (see the file header for
provenance); this script reads the rate coefficients from that same file.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "tools", "python"))
from mochii_output import HeiMetaData, RatesData          # noqa: E402

ATOM = os.path.join(HERE, "..", "..", "data", "atomic")
YHE = 0.1                                     # par%He_abund in all three runs
NH = {"hei_low_density": 100.0, "hei_high_density": 1.0e4,
      "hei_low_density_fuv": 100.0}
RCUT = {"hei_low_density": 1.0, "hei_high_density": 0.038}   # core radius [pc]


# --- rate coefficients from the data file (single source of truth) -----------
def parse_coeffs():
    c = {}
    with open(os.path.join(ATOM, "hei_metastable.txt")) as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            t = ln.split()
            c[t[0]] = [float(x) for x in t[1:]]
    return c


C = parse_coeffs()
A31 = C["A31"][0]
KB_EV = 1.380649e-16 / 1.602176634e-12


def alpha3(T):
    Cc, p = C["ALPHA3"]
    return Cc * (T / 1e4) ** p


def qexc(T, key):
    pref, sc, Ee, div, a, b, cc, d = C[key]
    kT = KB_EV * T
    ups = a * np.exp(b * T) + cc * np.exp(d * T)
    return pref * np.sqrt(sc / kT) * np.exp(-Ee / kT) * ups / div


T0 = 1.0e4                                    # par%te_fixed in all runs
ALPHA3_0 = alpha3(T0)                         # 2.10e-13 cm^3/s
Q31_0 = qexc(T0, "QEXC31A") + qexc(T0, "QEXC31B")   # 3.224e-8 cm^3/s
PLATEAU = ALPHA3_0 / Q31_0                    # alpha_3/(q_31a+q_31b) ~ 6.51e-6
NCRIT = A31 / Q31_0                           # A_31/(q_31a+q_31b) ~ 3945


# --- helpers -----------------------------------------------------------------
def load(base):
    """Return radius [pc], n_2s3, Phi_3, n_e, T_e, and n_HeII for a run.

    n_HeII = nH * He_abund * x_HeII with nH the run's uniform density and
    x_HeII from the rates file (leaf order matches the heimeta file)."""
    hm = HeiMetaData(os.path.join(HERE, base + "_heimeta.h5"), atomic_dir=ATOM)
    rt = RatesData(os.path.join(HERE, base + "_rates.h5"))
    r = np.sqrt((hm.xyz ** 2).sum(axis=0))          # code units = pc
    nH = np.full(hm.nleaf, NH[base])
    ne = hm.fields["n_e"]
    nH[ne <= 0.0] = 0.0                              # outside the rmax sphere
    xHeII = rt.fields["x_HeII"]
    nHeII = nH * YHE * xHeII
    return dict(hm=hm, r=r, n3=hm.fields["n_2s3"], Phi3=hm.fields["Phi_3"],
                ne=ne, Te=hm.fields["T_e"], nHeII=nHeII, nH=nH, xHeII=xHeII)


def radial_median(r, y, mask, nbin=24, rmax=None):
    """Median of y in nbin radial bins over the masked leaves."""
    rm = r[mask]
    ym = y[mask]
    if rmax is None:
        rmax = rm.max()
    edges = np.linspace(0.0, rmax, nbin + 1)
    cen, med = [], []
    for i in range(nbin):
        s = (rm >= edges[i]) & (rm < edges[i + 1])
        if s.sum() >= 3:
            cen.append(0.5 * (edges[i] + edges[i + 1]))
            med.append(np.median(ym[s]))
    return np.array(cen), np.array(med)


# --- figure 1: radial n_3/n_HeII profiles ------------------------------------
def fig_profiles(dl, dh):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

    for ax, d, base, title in (
            (axes[0], dl, "hei_low_density", "low density (radiative limit)"),
            (axes[1], dh, "hei_high_density",
             "high density (collisional plateau)")):
        # restrict to the ionized region (He+ present): where n_3/n_HeII is the
        # meaningful ratio.  Beyond the front n_HeII -> 0 and the ratio diverges.
        m = (d["nHeII"] > 0.0) & (d["ne"] > 0.0) & (d["xHeII"] > 0.5)
        ratio = np.full(d["n3"].shape, np.nan)
        ratio[m] = d["n3"][m] / d["nHeII"][m]
        radlim = alpha3(T0) * d["ne"] / A31        # alpha_3 n_e / A_31 per leaf

        ax.scatter(d["r"][m], ratio[m], s=2, alpha=0.15, color="0.6",
                   rasterized=True, label="leaves")
        cen, med = radial_median(d["r"], ratio, m)
        ax.plot(cen, med, "o-", color="C0", ms=4, lw=1.5,
                label=r"$n_3/n_{\rm He^+}$ (median)")
        # radiative limit alpha_3 n_e/A_31 vs radius (median of the leaf-by-leaf line)
        cenr, medr = radial_median(d["r"], radlim, m)
        ax.plot(cenr, medr, "--", color="C3", lw=1.5,
                label=r"radiative $\alpha_3 n_e/A_{31}$")
        ax.set_xlabel(r"$r$ [pc]")
        ax.set_title(title)
        ax.set_yscale("log")

    axes[0].set_ylabel(r"$n_3/n_{\rm He^+}$")
    # the high-density panel also shows the density-independent plateau
    axes[1].axhline(PLATEAU, color="C2", lw=1.5, ls=":",
                    label=r"plateau $\alpha_3/(q_{31a}+q_{31b})$")
    axes[0].legend(fontsize=8, loc="lower left")
    axes[1].legend(fontsize=8, loc="lower left")
    fig.suptitle(r"He I 2$^3$S metastable population "
                 r"($T_e = 10^4$ K, $n_{\rm crit} = " + f"{NCRIT:.0f}" +
                 r"$ cm$^{-3}$)")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(HERE, "hei_profiles.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


# --- figure 2: high-density maps ---------------------------------------------
def fig_maps(dh):
    hm = dh["hm"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    col, ext = hm.column_map("n_2s3", npix=256)
    tau, _ = hm.tau10830_map(component="doublet", npix=256, warn=True)

    for ax, img, lbl, title in (
            (axes[0], col, r"$\log_{10} N(2^3{\rm S})$ [cm$^{-2}$]",
             r"2$^3$S column density"),
            (axes[1], tau, r"$\log_{10}\tau_0$(10830 doublet)",
             r"10830 \AA\ line-center $\tau_0$")):
        show = np.log10(np.where(img > 0, img, np.nan))
        im = ax.imshow(show.T, origin="lower", extent=ext, cmap="inferno")
        fig.colorbar(im, ax=ax, label=lbl)
        ax.set_xlabel(r"$x$ [pc]")
        ax.set_ylabel(r"$y$ [pc]")
        ax.set_title(title)
    # mark the tau=1 contour: the optically-thin diagnostic breaks down inside
    axes[1].contour(np.log10(np.where(tau > 0, tau, np.nan)).T, levels=[0.0],
                    colors="cyan", linewidths=1.2, extent=ext)
    fig.suptitle(r"High-density run: the core is optically THICK at 10830 "
                 r"($\tau_0 \gg 1$, cyan = $\tau_0 = 1$)")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(HERE, "hei_maps.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out, float(np.nanmax(tau))


# --- figure 3: soft-UV Phi_3 depletion ---------------------------------------
def fig_phi3(dl, df):
    # identical grids (same nx/box) -> matched leaf order; restrict to the
    # ionized region where He+ is present (the front-cell ratio is noise).
    m = (dl["nH"] > 0.0) & (dl["n3"] > 0.0) & (dl["xHeII"] > 0.5)
    dep = np.full(dl["n3"].shape, np.nan)
    dep[m] = df["n3"][m] / dl["n3"][m]              # n_3^fuv / n_3^noFUV
    phi_ratio = df["Phi3"] / A31
    rmax = dl["r"][m].max()

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    cen, med = radial_median(dl["r"], dep, m, rmax=rmax)
    ax.plot(cen, med, "o-", color="C0", ms=4, lw=1.6,
            label=r"$n_3^{\rm FUV}/n_3^{\rm no\,FUV}$")
    ax.set_xlabel(r"$r$ [pc]")
    ax.set_ylabel(r"$n_3^{\rm FUV}/n_3^{\rm no\,FUV}$", color="C0")
    ax.tick_params(axis="y", labelcolor="C0")
    ax.set_ylim(0.0, 1.05)

    ax2 = ax.twinx()
    cen2, med2 = radial_median(df["r"], phi_ratio, m, rmax=rmax)
    ax2.plot(cen2, med2, "s--", color="C3", ms=4, lw=1.6,
             label=r"$\Phi_3/A_{31}$")
    ax2.set_ylabel(r"$\Phi_3/A_{31}$", color="C3")
    ax2.tick_params(axis="y", labelcolor="C3")

    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [l.get_label() for l in lines], fontsize=9,
              loc="center right")
    ax.set_title(r"Soft-UV photoionization depletes the 2$^3$S reservoir "
                 r"near the star")
    fig.tight_layout()
    out = os.path.join(HERE, "hei_phi3.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out, float(np.nanmax(phi_ratio)), np.nanmin(med)


# --- main --------------------------------------------------------------------
def main():
    dl = load("hei_low_density")
    dh = load("hei_high_density")
    df = load("hei_low_density_fuv")

    p1 = fig_profiles(dl, dh)
    p2, tau_max = fig_maps(dh)
    p3, phi_max, dep_min = fig_phi3(dl, df)

    # --- console summary -----------------------------------------------------
    print("\n=== He I 2^3S metastable example summary ===")
    print("coefficients @ 1e4 K (data/atomic/hei_metastable.txt):")
    print("  alpha_3 = %.3e cm^3/s,  A_31 = %.3e s^-1" % (ALPHA3_0, A31))
    print("  q_31a+q_31b = %.3e cm^3/s,  plateau alpha_3/(q31a+q31b) = %.3e"
          % (Q31_0, PLATEAU))
    print("  n_crit = A_31/(q31a+q31b) = %.0f cm^-3" % NCRIT)

    # low-density: collisional suppression vs the radiative limit in the core
    mc = (dl["r"] < RCUT["hei_low_density"]) & (dl["nHeII"] > 0.0) \
        & (dl["ne"] > 0.0)
    ratio_c = dl["n3"][mc] / dl["nHeII"][mc]
    radlim_c = alpha3(T0) * dl["ne"][mc] / A31
    suppr = 1.0 - ratio_c / radlim_c
    print("\nlow density (radiative limit):")
    print("  core <n_e> = %.1f cm^-3 (<< n_crit)" % dl["ne"][mc].mean())
    print("  collisional suppression below alpha_3 n_e/A_31: "
          "median %.1f%%, max %.1f%%"
          % (100 * np.median(suppr), 100 * suppr.max()))

    # high-density: fraction of the plateau in the core
    mh = (dh["r"] < RCUT["hei_high_density"]) & (dh["nHeII"] > 0.0) \
        & (dh["ne"] > 0.0)
    ratio_h = dh["n3"][mh] / dh["nHeII"][mh]
    frac = np.median(ratio_h) / PLATEAU
    radlim_h = alpha3(T0) * dh["ne"][mh] / A31
    below_rad = np.median(ratio_h / radlim_h)          # ~ 1/(1+n_e/n_crit)
    print("\nhigh density (collisional plateau):")
    print("  core <n_e> = %.2e cm^-3 (~%.1f n_crit)"
          % (dh["ne"][mh].mean(), dh["ne"][mh].mean() / NCRIT))
    print("  core n_3/n_HeII = %.3f of the radiative limit alpha_3 n_e/A_31 "
          "(factor %.1f suppression)" % (below_rad, 1.0 / below_rad))
    print("  core n_3/n_HeII = %.3f of the plateau alpha_3/(q31a+q31b)" % frac)
    print("  max 10830 doublet tau_0 through the core = %.2e (optically THICK)"
          % tau_max)

    # FUV: soft-UV depletion near the star
    print("\nlow density + soft-UV FUV (Phi_3 sink):")
    print("  max Phi_3/A_31 = %.3f" % phi_max)
    print("  strongest n_3 depletion (n_3^FUV/n_3^noFUV min) = %.2f"
          % dep_min)

    print("\nfigures written:")
    for p in (p1, p2, p3):
        print("  " + p)


if __name__ == "__main__":
    main()
