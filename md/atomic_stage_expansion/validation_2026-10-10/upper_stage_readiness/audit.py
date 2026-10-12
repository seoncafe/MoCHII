"""Read-only inventory of newly staged oxygen-to-calcium excitation ions."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
from chianti_cooling import BT_SOURCE_TMIN_K, bound_excitation_atom  # noqa: E402
from write_excitation_atom import BT_SOURCE_TMAX_K, _source_files  # noqa: E402

STAGED = ROOT / "md/atomic_stage_expansion/validation_2026-10-10/langevin64_exit_heat/products/atomic"
PRODUCTION = ROOT / "data/atomic"
VALIDATION = ROOT / "md/atomic_stage_expansion/validation_2026-10-10"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stage_count(path):
    return int(next(row for row in path.read_text().splitlines()
                    if row.startswith("ELEMENT ")).split()[3])


def roman(number):
    digits = ((10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I"))
    name = ""
    for value, symbol in digits:
        while number >= value:
            name += symbol
            number -= value
    return name


def report():
    products = list((ROOT / "md/atomic_stage_expansion/validation_2026-10-09").rglob("*.txt"))
    products += list(VALIDATION.rglob("*.txt"))
    validations = list((ROOT / "md/atomic_stage_expansion/validation_2026-10-09").rglob("*.json"))
    validations += list(VALIDATION.rglob("*.json"))
    records = []
    for installed in sorted(PRODUCTION.glob("element_*.txt")):
        element = installed.stem[8:]
        staged = STAGED / installed.name
        for stage in range(stage_count(installed) + 1, stage_count(staged) + 1):
            levels, radiation, collisions = bound_excitation_atom(element, stage)
            sources = _source_files(element, stage)
            source_records = [{"path": str(p), "sha256": digest(p)} for p in sources]
            stem = f"{element}_{stage}"
            atom_products = sorted(str(p.relative_to(ROOT)) for p in products
                                   if p.name == f"nlevel_{stem}.txt")
            cooling_products = sorted(str(p.relative_to(ROOT)) for p in products
                                      if p.name == f"cooling_{stem}.txt")
            source_header = []
            scups = next((p for p in sources if p.suffix == ".scups"), None)
            if scups:
                for line in scups.read_text().splitlines():
                    if re.search(r"temperature range|temperatures? provided in the range|collision strengths provided in the range", line, re.I):
                        source_header.append(line.strip())
            node_count = sum(len(c["xs"]) for c in collisions)
            bt_floor = BT_SOURCE_TMIN_K.get((element, stage))
            ion_name = element.capitalize() + roman(stage)
            direct_checks = []
            for check in validations:
                if ion_name not in check.parts and check.stem != ion_name:
                    continue
                if check.name not in ("result.json", "native_result.json", f"{ion_name}.json"):
                    continue
                try:
                    data = json.loads(check.read_text())
                except (OSError, ValueError):
                    continue
                direct_checks.append({"path": str(check.relative_to(ROOT)),
                                      "passed": data.get("passed"),
                                      "native_accepted": data.get("native_accepted"),
                                      "model_accepted": data.get("model_accepted")})
            record = {
                "ion": stem, "production_stages": stage_count(installed),
                "staged_stages": stage_count(staged), "levels": len(levels),
                "radiative_rows": len(radiation), "collision_rows": len(collisions),
                "collision_nodes": node_count, "source_files": source_records,
                "source_header_temperature_notes": source_header,
                "explicit_BT_floor_K": bt_floor,
                "explicit_BT_upper_K": BT_SOURCE_TMAX_K.get((element, stage)),
                "installed_atom": (PRODUCTION / f"nlevel_{stem}.txt").exists(),
                "installed_fit": (PRODUCTION / f"cooling_{stem}.txt").exists(),
                "private_atom_products": atom_products,
                "private_cooling_products": cooling_products,
                "direct_check_results": direct_checks,
            }
            records.append(record)
    return {
        "installed_element_files": {p.name: digest(p) for p in sorted(PRODUCTION.glob("element_*.txt"))},
        "staged_element_files": {p.name: digest(p) for p in sorted(STAGED.glob("element_*.txt"))},
        "ions": records,
    }


if __name__ == "__main__":
    output = Path(__file__).with_name("result.json")
    if output.exists():
        raise FileExistsError("Preserve the existing source-readiness snapshot")
    output.write_text(json.dumps(report(), indent=2) + "\n")
    print(f"wrote {output}")
