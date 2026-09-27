#!/usr/bin/env python3
"""He+ + H -> He + H+ + photon (radiative charge transfer): the two power-law
fits in use against their common source and against the recent calculation.

  GJ07 R26   1.25e-15 (T/300)^0.25   Glover & Jappsen (2007) Table 1, citing
                                     Zygelman et al. (1989); MoCHII until 2026-09-26
  SLD98 (19) 1.20e-15 (T/300)^0.25   Stancil, Lepp & Dalgarno (1998) Table 1,
                                     "Zygelman et al. 1989, multiplied by 0.25 to
                                     account for approach probability factor"; EXHALE
  Zygelman et al. (1989) Table I "Direct" column x 1/4 (1-1000 K).
  Courtney et al. (2021) Figure 3, RCT in the NLTE zero-density limit without
  stimulated processes (red solid), digitized from the published figure
  (courtney2021_fig3_rct_digitized.txt: the lowest red pixel in each column of
  a 400-dpi rendering; axes calibrated on the major ticks, 194.2 px/dex in T
  and 228.8 px/dex in k; line width ~3 px, so ~3% per point).  Their
  tabulated data (UGA Molecular Opacity Project website) could not be
  reached on 2026-09-26.  Above ~3e4 K the solid and dashed (stimulated)
  curves separate and the lowest red pixel may not follow the solid one.

Run:  python3 radiative_channel_comparison.py   (writes radiative_channel_comparison.png)
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
c = np.loadtxt(os.path.join(HERE, "courtney2021_fig3_rct_digitized.txt"))
T = np.logspace(0, 6, 400)
zT = np.array([1, 10, 100, 200, 400, 1000.0])
zk = np.array([5.36, 4.21, 4.41, 4.50, 4.83, 5.99])*1e-15/4

fig, ax = plt.subplots(2, 1, figsize=(6.5, 7), sharex=True,
                       gridspec_kw=dict(height_ratios=[2, 1]))
ax[0].loglog(c[:, 0], c[:, 1], "r-", lw=1.5, label="Courtney et al. (2021), Fig. 3 (digitized)")
ax[0].loglog(T, 1.25e-15*(T/300)**0.25, "b--", label=r"GJ07 R26, $1.25\times10^{-15}(T/300)^{0.25}$")
ax[0].loglog(T, 1.20e-15*(T/300)**0.25, "k-", lw=1, label=r"SLD98 (19), $1.20\times10^{-15}(T/300)^{0.25}$")
ax[0].loglog(zT, zk, "go", label=r"Zygelman et al. (1989) Table I $\times$ 1/4")
ax[0].set_ylabel(r"$k$ [cm$^3$ s$^{-1}$]")
ax[0].legend(fontsize=7)
ci = np.interp(np.log10(T), np.log10(c[:, 0]), c[:, 1])
ax[1].semilogx(T, 1.25e-15*(T/300)**0.25/ci, "b--", label="GJ07 / Courtney")
ax[1].semilogx(T, 1.20e-15*(T/300)**0.25/ci, "k-", lw=1, label="SLD98 / Courtney")
ax[1].axhline(1, c="0.5", lw=0.5)
ax[1].axvspan(5e3, 2e4, color="0.9")
ax[1].set_ylim(0.3, 1.4)
ax[1].set_ylabel("ratio")
ax[1].set_xlabel(r"$T$ [K]")
ax[1].legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "radiative_channel_comparison.png"), dpi=110)
for t in (1e3, 3e3, 5e3, 8e3, 1e4, 1.5e4, 2e4, 3e4):
    k = np.interp(np.log10(t), np.log10(c[:, 0]), c[:, 1])
    print(f"{t:8.0f}  Courtney {k:.3e}  GJ07/C {1.25e-15*(t/300)**.25/k:.3f}  SLD98/C {1.20e-15*(t/300)**.25/k:.3f}")
