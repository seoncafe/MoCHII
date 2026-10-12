"""Plot all saved Ne I populations without recalculating physical products."""

from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
native_path = HERE / "native_profiles.txt"
direct_path = HERE / "direct_profiles.npz"
manifest = json.loads((HERE / "validation_manifest.json").read_text())
for path, key in ((native_path, "native_profile_sha256"),
                  (direct_path, "direct_profile_sha256")):
    if hashlib.sha256(path.read_bytes()).hexdigest() != manifest[key]:
        raise ValueError("Saved validated population profile changed")
native = np.loadtxt(native_path).reshape(401, 4, 34)
direct = np.load(direct_path)
fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
colors = plt.get_cmap("viridis")(np.linspace(0, 1, 31))
positive_populations = native[:, :, 3:][native[:, :, 3:] > 0]
population_floor = 0.1 * float(positive_populations.min())
for density_index, ax in enumerate(axes.flat):
    temperature = direct["temperature_K"]
    population = direct[f"populations_{density_index}"]
    for level in range(31):
        ax.loglog(temperature, population[:, level], color=colors[level],
                  linewidth=0.7)
        ax.loglog(temperature, native[:, density_index, level + 3],
                  color=colors[level], linestyle="--", linewidth=0.7)
    density = float(direct["density_cm3"][density_index])
    ax.set_title(r"$n_e = " + f"{density:g}" + r"\ {\rm cm}^{-3}$")
    ax.set_xlabel(r"$T$ [K]")
    ax.set_ylabel(r"Level population fraction")
    ax.set_ylim(population_floor, 1.5)
fig.suptitle("All 31 states: raw-source calculation (solid), compiled solver (dashed)")
fig.tight_layout(rect=(0, 0, 0.92, 0.96))
bar_ax = fig.add_axes([0.93, 0.14, 0.02, 0.68])
bar = fig.colorbar(matplotlib.cm.ScalarMappable(
    norm=matplotlib.colors.Normalize(1, 31), cmap="viridis"), cax=bar_ax)
bar.set_label("Source level index")
bar.set_ticks([1, 10, 20, 31])
output = HERE / "all_level_populations.pdf"
if output.exists():
    original = HERE / "all_level_populations_initial_axis.pdf"
    if original.exists():
        raise FileExistsError("Preserve existing figure versions")
    output.rename(original)
fig.savefig(output)
plt.close(fig)
