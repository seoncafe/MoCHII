"""Plot source-fit Fe V--VII recombination curves for audit only.

Usage: MPLCONFIGDIR=/tmp/mochii-mpl python tools/fitting/plot_fe_recombination_expansion.py
The source fits have unverified temperature validity and are not production rates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fe_recombination_expansion import (
    DEFAULT_CHIANTI,
    DEFAULT_CLOUDY,
    chianti_rates,
    cloudy_gu_fe_vii_dr,
    read_chianti_fe,
    read_cloudy_gu_fe_vii,
)


DEFAULT_OUTPUT = Path(__file__).resolve().parents[2] / "md/atomic_stage_expansion/fe_recombination_expansion"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_curves(chianti_fe: Path, cloudy_source: Path, count: int) -> dict:
    if count < 2 or count > 10001:
        raise ValueError("Grid size must be from 2 to 10001")
    temperatures = np.geomspace(1e3, 1e6, count)
    gu = read_cloudy_gu_fe_vii(cloudy_source)
    records = {stage: read_chianti_fe(chianti_fe, stage) for stage in (5, 6, 7)}
    data = {
        "temperature_K": temperatures.tolist(),
        "rate_unit": "cm^3 s^-1",
        "temperature_validity": "unverified for these source fits",
        "source_fits_only": True,
        "chianti": {},
        "cloudy_gu_fe_vii": {
            "source": str(cloudy_source),
            "sha256": sha256(cloudy_source),
            "provenance": gu["provenance"],
            "dr_cm3_s": [cloudy_gu_fe_vii_dr(gu, float(t)) for t in temperatures],
        },
    }
    for stage, record in records.items():
        rates = [chianti_rates(record, float(t)) for t in temperatures]
        sources = [Path(record["rr_source"]), Path(record["dr_source"])]
        data["chianti"][str(stage)] = {
            "transition": record["transition"],
            "sources": [{"path": str(p), "sha256": sha256(p)} for p in sources],
            "rr_cm3_s": [r[0] for r in rates],
            "dr_cm3_s": [r[1] for r in rates],
        }
    return data


def plot_curves(data: dict, output: Path) -> None:
    temperatures = np.asarray(data["temperature_K"])
    colors = {5: "#0072B2", 6: "#009E73", 7: "#D55E00"}
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 8.2), sharex=True,
                             constrained_layout=True)
    for stage in (5, 6, 7):
        rates = data["chianti"][str(stage)]
        axes[0].loglog(temperatures, rates["rr_cm3_s"], color=colors[stage],
                       linewidth=1.8, label=f"CHIANTI Fe {stage} RR")
        axes[1].loglog(temperatures, rates["dr_cm3_s"], color=colors[stage],
                       linewidth=1.8, label=f"CHIANTI Fe {stage} DR")
    axes[1].loglog(temperatures, data["cloudy_gu_fe_vii"]["dr_cm3_s"],
                   color="#7B3294", linestyle="--", linewidth=2,
                   label="Cloudy GuPC Fe 7 DR")
    axes[0].set_ylabel(r"RR coefficient (cm$^3$ s$^{-1}$)")
    axes[1].set_ylabel(r"DR coefficient (cm$^3$ s$^{-1}$)")
    axes[1].set_xlabel("Electron temperature (K)")
    axes[0].set_ylim(1e-13, 3e-10)
    # CHIANTI's tiny low-temperature values are displayed, not clipped.
    axes[1].set_ylim(1e-300, 1e-9)
    axes[1].set_xlim(1e3, 1e6)
    for ax in axes:
        ax.grid(True, which="major", alpha=0.25)
        ax.legend(loc="best", fontsize=8)
    fig.suptitle("Fe V-VII source-fit recombination curves\n"
                 "Temperature validity of source fits is unverified.", fontsize=10)
    fig.savefig(output, format="pdf")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chianti-fe", type=Path, default=DEFAULT_CHIANTI)
    parser.add_argument("--cloudy-source", type=Path, default=DEFAULT_CLOUDY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--count", type=int, default=301)
    args = parser.parse_args()
    data = make_curves(args.chianti_fe, args.cloudy_source, args.count)
    if any(not math.isfinite(x) or x <= 0 for block in data["chianti"].values()
           for key in ("rr_cm3_s", "dr_cm3_s") for x in block[key]):
        raise ValueError("Nonpositive or nonfinite CHIANTI rate on grid")
    if any(not math.isfinite(x) or x <= 0 for x in data["cloudy_gu_fe_vii"]["dr_cm3_s"]):
        raise ValueError("Nonpositive or nonfinite Cloudy rate on grid")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "fe_recombination_curves.json").write_text(
        json.dumps(data, indent=2) + "\n")
    plot_curves(data, args.output_dir / "fe_recombination_curves.pdf")


if __name__ == "__main__":
    main()
