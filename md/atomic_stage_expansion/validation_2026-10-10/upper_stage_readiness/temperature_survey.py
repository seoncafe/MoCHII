"""Read original collision metadata and infer interior BT node temperatures."""

from __future__ import annotations

import json
import math
from pathlib import Path

from audit import ROOT, report
from chianti_cooling import KB_OVER_RY, read_scups


def first_last_finite_temperatures(row):
    xs = [float(x) for x in row["xs"] if 0.0 < x < 1.0]
    if not xs:
        return None
    c = float(row["cups"])
    gap = float(row["de"])
    def inverse(x):
        et = (math.exp(math.log(c) / (1.0 - x)) - c
              if row["ttype"] in (1, 4) else c * x / (1.0 - x))
        return et * gap / KB_OVER_RY
    try:
        lo, hi = inverse(xs[0]), inverse(xs[-1])
    except (OverflowError, ValueError, ZeroDivisionError):
        return None
    return (lo, hi) if 0 < lo <= hi and math.isfinite(hi) else None


def survey():
    entries = []
    for atom in report()["ions"]:
        ion = atom["ion"]
        sources = [Path(x["path"]) for x in atom["source_files"]]
        scups = next((x for x in sources if x.suffix == ".scups"), None)
        if scups:
            rows = read_scups(scups)
            temperatures = [first_last_finite_temperatures(row) for row in rows]
            temperatures = [row for row in temperatures if row is not None]
            lines = scups.read_text().splitlines()
            marker = next((j for j, row in enumerate(lines) if row.strip() == "-1"), None)
            footer = lines[marker + 1:] if marker is not None else []
            entries.append({
                "ion": ion, "source_path": str(scups), "source_sha256":
                next(x["sha256"] for x in atom["source_files"] if x["path"] == str(scups)),
                "source_footer": [line.strip() for line in footer if line.strip() and line.strip() != "-1"],
                "interior_BT_rows": len(temperatures),
                "inferred_first_interior_T_K_range": [min(x[0] for x in temperatures),
                                                      max(x[0] for x in temperatures)],
                "inferred_last_interior_T_K_range": [min(x[1] for x in temperatures),
                                                     max(x[1] for x in temperatures)],
                "inference_limit": "Rounded BT spline coordinates, not necessarily original source grid endpoints or published accuracy bounds",
            })
        else:
            from chianti_cooling import bound_excitation_atom
            element, stage = ion.split("_")
            _, _, collisions = bound_excitation_atom(element, int(stage))
            grids = [row["temperature_grid"] for row in collisions]
            entries.append({"ion": ion, "source_path": str(next(x for x in sources if x.suffix == ".coll")),
                            "source_sha256": next(x["sha256"] for x in atom["source_files"] if x["path"].endswith(".coll")),
                            "Stout_first_T_K_range": [min(float(x[0]) for x in grids), max(float(x[0]) for x in grids)],
                            "Stout_last_T_K_range": [min(float(x[-1]) for x in grids), max(float(x[-1]) for x in grids)]})
    return entries


if __name__ == "__main__":
    path = Path(__file__).with_name("temperature_result.json")
    path.write_text(json.dumps(survey(), indent=2) + "\n")
    print(f"wrote {path}")
