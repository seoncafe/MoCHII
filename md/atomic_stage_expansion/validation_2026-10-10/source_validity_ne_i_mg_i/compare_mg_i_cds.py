"""Compare CHIANTI Mg I fine-state strengths with published BSR LS tables.

The CDS inputs are external scientific source data. Supply a local private copy
with --cds-dir; this script does not redistribute their numerical contents.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
from chianti_cooling import read_excitation_atom, upsilon  # noqa: E402


LS_LEVELS = {
    1: (1,),
    2: (2, 3, 4),
    3: (5,),
    4: (6,),
    5: (7,),
    6: (8,),
    7: (9, 10, 11),
    8: (12, 13, 14),
    9: (15,),
}
TEMPERATURES = (500, 1000, 2000, 3000, 4000, 5000, 6000, 7000,
                8000, 9000, 10000, 15000, 20000)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cds-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    args = parser.parse_args()
    _, _, collisions = read_excitation_atom("mg", 1)
    measurements = []
    hashes = {}
    for temperature in TEMPERATURES:
        path = args.cds_dir / f"{temperature}k_b.dat"
        raw = path.read_bytes()
        hashes[path.name] = hashlib.sha256(raw).hexdigest()
        bsr = np.array([[float(value) for value in line.split()]
                        for line in raw.decode("ascii").splitlines()])
        if bsr.shape != (25, 25) or not np.all(np.isfinite(bsr)):
            raise ValueError(f"Invalid CDS matrix: {path}")
        ratios = []
        pairs = []
        for lower in LS_LEVELS:
            for upper in LS_LEVELS:
                if upper <= lower:
                    continue
                selected = sum(float(upsilon(row, temperature)[0])
                               for row in collisions
                               if (row["ll"] in LS_LEVELS[lower]
                                   and row["ul"] in LS_LEVELS[upper])
                               or (row["ll"] in LS_LEVELS[upper]
                                   and row["ul"] in LS_LEVELS[lower]))
                published = bsr[lower - 1, upper - 1]
                if not selected > 0 or not published > 0:
                    raise ValueError(f"Nonpositive strength at {temperature} K")
                ratio = selected / published
                ratios.append(ratio)
                pairs.append([lower, upper, ratio])
        measurements.append({
            "temperature_K": temperature,
            "number_of_LS_pairs": len(ratios),
            "max_absolute_fractional_difference": float(np.max(np.abs(np.array(ratios) - 1))),
            "mean_absolute_fractional_difference": float(np.mean(np.abs(np.array(ratios) - 1))),
            "smallest_ratio": float(min(ratios)),
            "largest_ratio": float(max(ratios)),
            "worst_pair": max(pairs, key=lambda pair: abs(pair[2] - 1)),
            "LS_pair_ratios": pairs,
        })
    result = {"source": "CDS J/A+A/606/A11 BSR matrices", "sha256": hashes,
              "temperatures": measurements}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    temperatures = np.array(TEMPERATURES, dtype=float)
    ratios_by_temperature = np.array([[pair[2] for pair in item["LS_pair_ratios"]]
                                      for item in measurements])
    figure, axes = plt.subplots(1, 2, figsize=(10, 4))
    for column in range(ratios_by_temperature.shape[1]):
        axes[0].plot(temperatures, ratios_by_temperature[:, column],
                     linewidth=0.7, alpha=0.7)
    axes[0].axhline(1, color="black", linewidth=0.8)
    axes[0].set(xlabel="Electron temperature [K]",
                ylabel=r"$\sum_{J_i,J_j}\Upsilon_{\rm CHIANTI}/\Upsilon_{\rm BSR}$",
                xscale="log", yscale="log")
    axes[1].plot(temperatures,
                 [item["max_absolute_fractional_difference"] for item in measurements],
                 label="Largest pair difference")
    axes[1].plot(temperatures,
                 [item["mean_absolute_fractional_difference"] for item in measurements],
                 label="Mean over 36 pairs")
    axes[1].set(xlabel="Electron temperature [K]", ylabel="Absolute fractional difference",
                xscale="log", yscale="log")
    axes[1].legend()
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure)
    plt.close(figure)
    for item in measurements:
        print(item["temperature_K"], item["max_absolute_fractional_difference"],
              item["mean_absolute_fractional_difference"], item["worst_pair"])


if __name__ == "__main__":
    main()
