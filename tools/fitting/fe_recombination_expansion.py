"""Audit candidate Fe V--VII recombination rates without changing production data.

CHIANTI stage s denotes the *recombining* Fe ion: Fe s -> Fe (s-1).
The Cloudy GuPC fit begins at Fe VII -> Fe VI. Neither source establishes a
validated low-temperature DR rate for all three transitions.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


DEFAULT_CHIANTI = Path.home() / "RT_Codes/CHIANTI/dbase/fe"
DEFAULT_CLOUDY = Path.home() / "CLOUDY/c25.00/source/ion_recomb_Badnell.cpp"
EV_TO_K = 11604.518121550082


def read_chianti_fe(fe_dir: Path, stage: int) -> dict:
    """Read CHIANTI type-3 RR and type-1 DR for Fe V, VI, or VII."""
    if stage not in (5, 6, 7):
        raise ValueError("Only Fe V--VII are in this audit")
    prefix = fe_dir / f"fe_{stage}" / f"fe_{stage}"
    rr_lines = Path(f"{prefix}.rrparams").read_text().splitlines()
    dr_lines = Path(f"{prefix}.drparams").read_text().splitlines()
    if int(rr_lines[0]) != 3 or int(dr_lines[0]) != 1:
        raise ValueError(f"Unexpected CHIANTI fit type for Fe {stage}")
    rr = rr_lines[1].split()
    energy = dr_lines[1].split()
    coef = dr_lines[2].split()
    if [(int(row[0]), int(row[1])) for row in (rr, energy, coef)] != [(26, stage)] * 3:
        raise ValueError(f"Incorrect charge indexing in Fe {stage} data")
    energies = tuple(map(float, energy[2:]))
    coefficients = tuple(map(float, coef[2:]))
    if len(energies) != len(coefficients) or not energies:
        raise ValueError(f"Malformed CHIANTI DR terms for Fe {stage}")
    if any(e < 0 or c < 0 for e, c in zip(energies, coefficients)):
        raise ValueError("Negative CHIANTI DR energy or coefficient")
    if any((e == 0) != (c == 0) for e, c in zip(energies, coefficients)):
        raise ValueError("Unpaired CHIANTI DR term")
    return {
        "stage": stage,
        "transition": f"Fe {stage} -> Fe {stage - 1}",
        "rr_A": float(rr[2]),
        "rr_eta": float(rr[3]),
        "dr_E_K": energies,
        "dr_c_cm3_s_K3_2": coefficients,
        "rr_source": str(prefix) + ".rrparams",
        "dr_source": str(prefix) + ".drparams",
    }


def chianti_rates(record: dict, temperature_K: float) -> tuple[float, float]:
    """Return RR and DR [cm^3 s^-1] using the CHIANTI file fit conventions."""
    if not math.isfinite(temperature_K) or temperature_K <= 0:
        raise ValueError("Temperature must be finite and positive")
    rr = record["rr_A"] * (temperature_K / 1.0e4) ** (-record["rr_eta"])
    dr = math.fsum(
        c * math.exp(-e / temperature_K)
        for e, c in zip(record["dr_E_K"], record["dr_c_cm3_s_K3_2"])
    ) / temperature_K ** 1.5
    return rr, dr


def read_cloudy_gu_fe_vii(path: Path) -> dict:
    """Extract only the first GuPC row, labeled Fe+6 in Cloudy."""
    source = path.read_text()
    arrays = []
    for name in ("Fe_Gu_c", "Fe_Gu_E"):
        match = re.search(rf"static const double {name}\[9\]\[6\]\s*=\s*\{{\s*\{{([^}}]+)\}}", source)
        if match is None:
            raise ValueError(f"Cloudy {name} table not found")
        arrays.append(tuple(float(x.strip()) for x in match.group(1).split(",") if x.strip()))
    c, e = arrays
    if len(c) != 6 or len(e) != 6 or any(x < 0 for x in c + e):
        raise ValueError("Malformed Cloudy GuPC row")
    if "GuPC: Gu, M. private communication" not in source:
        raise ValueError("Cloudy GuPC provenance marker missing")
    return {"stage": 7, "transition": "Fe 7 -> Fe 6", "c_cm3_s_eV3_2": c,
            "E_eV": e, "source": str(path), "provenance": "Gu private communication (Cloudy label)"}


def cloudy_gu_fe_vii_dr(record: dict, temperature_K: float) -> float:
    """Evaluate Cloudy's GuPC Fe+6 DR expression in eV units [cm^3 s^-1]."""
    if not math.isfinite(temperature_K) or temperature_K <= 0:
        raise ValueError("Temperature must be finite and positive")
    temperature_eV = temperature_K / EV_TO_K
    return math.fsum(
        c * math.exp(-e / temperature_eV)
        for c, e in zip(record["c_cm3_s_eV3_2"], record["E_eV"])
    ) / temperature_eV ** 1.5


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chianti-fe", type=Path, default=DEFAULT_CHIANTI)
    parser.add_argument("--cloudy-source", type=Path, default=DEFAULT_CLOUDY)
    parser.add_argument("--temperatures", type=float, nargs="+", default=[1e3, 1e4, 1e5, 1e6])
    args = parser.parse_args()
    gu = read_cloudy_gu_fe_vii(args.cloudy_source)
    report = {"chianti": {}, "cloudy_gu_fe_vii": {}}
    for stage in (5, 6, 7):
        record = read_chianti_fe(args.chianti_fe, stage)
        report["chianti"][str(stage)] = {
            "source": [record["rr_source"], record["dr_source"]],
            "rates_cm3_s": {str(t): dict(zip(("rr", "dr"), chianti_rates(record, t)))
                            for t in args.temperatures},
        }
    report["cloudy_gu_fe_vii"] = {
        "source": gu["source"], "provenance": gu["provenance"],
        "dr_cm3_s": {str(t): cloudy_gu_fe_vii_dr(gu, t) for t in args.temperatures},
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
