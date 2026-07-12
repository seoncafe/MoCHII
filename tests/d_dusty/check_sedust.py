"""Check the SEDust smoke output: energy consistency against Heat_dust*V
in the rates file, PAH features, and a figure for the docs.

Run:  python3 check_sedust.py
"""

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "sedust_smoke"

# --- dust SED ------------------------------------------------------------
lam, Llam = np.loadtxt(f"{BASE}_dustsed.txt", unpack=True)
Ltot = np.trapz(Llam, lam)

# --- absorbed power from the rates file ----------------------------------
with h5py.File(f"{BASE}_rates.h5", "r") as f:
    hd = f["Heat_dust"][:]
    xyz = f["LeafXYZ"][:]
    # leaf volume from the level structure is not stored; recover from
    # the dustsed header instead (written as sum Heat_dust V).
Labs = None
with open(f"{BASE}_dustsed.txt") as fh:
    for ln in fh:
        if "integral L_lambda" in ln:
            Labs = float(ln.split("=")[1].split()[0])
            break

print(f"total dust luminosity  int L_lambda dlambda = {Ltot:.4e} erg/s")
print(f"header value (= sum Heat_dust V)            = {Labs:.4e} erg/s")
print(f"trapezoid/header ratio = {Ltot/Labs:.4f}")

# PAH features: local maxima near 6.2, 7.7, 11.3, 17 um
nuLnu = Llam * lam
for w in (6.2, 7.7, 11.3, 17.0):
    i = np.argmin(np.abs(lam - w))
    print(f"lambda*L_lambda at {w:5.1f} um : {nuLnu[i]:.3e} erg/s")
ipk = np.argmax(nuLnu)
print(f"peak of lambda*L_lambda at {lam[ipk]:.1f} um")

# --- figure ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.2, 3.6))
ax.loglog(lam, nuLnu, color="k", lw=1.0)
ax.set_xlim(1, 3e3)
sel = (lam > 1) & (lam < 3e3)
ax.set_ylim(nuLnu[sel].max() * 1e-4, nuLnu[sel].max() * 3)
ax.set_xlabel(r"$\lambda\ [\mu{\rm m}]$")
ax.set_ylabel(r"$\lambda L_\lambda\ [{\rm erg\,s^{-1}}]$")
for w in (6.2, 7.7, 11.3, 17.0):
    ax.axvline(w, color="0.7", lw=0.5, zorder=0)
ax.set_title(r"MoCHII + SEDust: stochastic dust SED (smoke test)")
fig.tight_layout()
fig.savefig("sedust_smoke_sed.png", dpi=200)
print("figure written: sedust_smoke_sed.png")
