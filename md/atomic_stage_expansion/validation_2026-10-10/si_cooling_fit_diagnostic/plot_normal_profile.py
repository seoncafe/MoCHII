"""Plot the accepted-accuracy S I staged fit and full residual profile."""
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as C
import fit_cooling as F

table = HERE / "normal_generator/cooling_s_1.txt"
rows = [line.split() for line in table.read_text().splitlines()
        if line.strip() and not line.startswith("#")]
assert int(rows[0][0]) == len(rows) - 1
amplitude = np.array([float(row[0]) for row in rows[1:]])
excitation = np.array([float(row[1]) for row in rows[1:]])
temperature = np.logspace(3, 5, 1601)
source = C.cooling_low_density_limit("s", 1, temperature)
fitted = F.evaluate(amplitude, excitation, temperature)
residual = 100 * (fitted / source - 1)

figure, axes = plt.subplots(2, 1, figsize=(6.3, 6.4), sharex=True,
                           gridspec_kw={"height_ratios": [2, 1]})
axes[0].loglog(temperature, source, "k-", label="CHIANTI source")
axes[0].loglog(temperature, fitted, "r--", label="21-term fit")
axes[0].set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
axes[0].legend(frameon=False)
axes[1].semilogx(temperature, residual, "b-", linewidth=1)
axes[1].axhline(0, color="0.4", linewidth=0.7)
axes[1].axhline(1, color="0.6", linewidth=0.7, linestyle=":")
axes[1].axhline(-1, color="0.6", linewidth=0.7, linestyle=":")
axes[1].set_ylim(-1.1, 1.1)
axes[1].set_xlabel(r"$T$ [K]")
axes[1].set_ylabel("Relative error [%]")
for axis in axes:
    axis.set_xlim(1e3, 1e5)
figure.tight_layout()
figure.savefig(HERE / "normal_generator/si_cooling_profile.pdf")
plt.close(figure)
