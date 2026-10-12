"""Construct isolated full-bound oxygen atoms with identified E1 strengths.

These files are physical candidates for validation, not installed atomic data.
Same-parity M1/E2 alternatives test conditional frequency sensitivity only.
"""

import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools" / "fitting"))
from atomic_excitation_energies import selected_energies  # noqa: E402
from chianti_cooling import read_excitation_atom  # noqa: E402
from write_excitation_atom import write_excitation_atom  # noqa: E402

EINSTEIN = 6.6703e15
HERE = Path(__file__).parent
OUTPUT = HERE / "candidate_atoms"
SZ_OPTICAL = {(2, 4): 6.791e-3, (3, 4): 2.046e-2}
result = {}

for stage, protected_prefix in ((2, 5), (3, 6)):
    raw, radiation, collisions = read_excitation_atom("o", stage)
    selected = selected_energies("o", stage, raw, radiation, collisions)
    source_path = Path.home() / f"RT_Codes/CHIANTI/dbase/o/o_{stage}/o_{stage}.elvlc"
    configs = {}
    for line in source_path.read_text().splitlines():
        if line.strip() == "-1":
            break
        configs[int(line[:7])] = line[7:37]

    def parity(index):
        occupations = re.findall(r"(\d+)([spdfgh])(\d*)", configs[index])
        if not occupations:
            raise ValueError(f"Missing O {stage} level {index} configuration")
        return sum((int(n) if n else 1) * "spdfgh".index(orbital)
                   for _, orbital, n in occupations) % 2

    records = []
    rates = {name: {} for name in ("base", "M1", "E2")}
    for lower, upper, wavelength, gf, source_rate in radiation:
        if source_rate <= 0:
            continue
        physical_lower, physical_upper = sorted((lower, upper), key=lambda i: selected[i][1])
        gap = selected[physical_upper][1] - selected[physical_lower][1]
        if gap <= 0:
            raise ValueError(f"O {stage} degenerate radiative pair {lower}-{upper}")
        opposite = parity(lower) != parity(upper)
        theory_frequency = raw[upper][2] - raw[lower][2]
        wavelength_frequency = 1e8 / abs(wavelength)
        if opposite:
            if gf <= 0:
                raise ValueError(f"O {stage} unclassified opposite-parity rate {lower}-{upper}")
            jl = (raw[lower][0] - 1) / 2
            ju = (raw[upper][0] - 1) / 2
            if abs(jl - ju) > 1 or (jl == 0 and ju == 0):
                raise ValueError(f"O {stage} non-E1 angular momentum {lower}-{upper}")
            inferred = 1e8 * (source_rate * raw[upper][0] / (EINSTEIN * gf)) ** .5
            possibilities = [("wavelength", wavelength_frequency)]
            if theory_frequency > 0:
                possibilities.append(("theoretical", theory_frequency))
            origin, source_frequency = min(possibilities, key=lambda item: abs(item[1] / inferred - 1))
            frequency_error = abs(source_frequency / inferred - 1)
            if frequency_error > .01:
                raise ValueError(f"O {stage} E1 frequency identity fails {lower}-{upper}")
            base = (source_rate * (gap / source_frequency) ** 3
                    * raw[upper][0] / raw[physical_upper][0])
        else:
            inferred = None
            origin, source_frequency = "wavelength", wavelength_frequency
            frequency_error = None
            base = source_rate
        pair = (physical_lower, physical_upper)
        values = {"base": base, "M1": base, "E2": base}
        if not opposite and max(pair) > protected_prefix:
            values["M1"] = base * (gap / source_frequency) ** 3
            values["E2"] = base * (gap / source_frequency) ** 5
        for name, value in values.items():
            rates[name][pair] = rates[name].get(pair, 0.0) + value
        records.append({"source_lower": lower, "source_upper": upper,
                        "physical_lower": physical_lower, "physical_upper": physical_upper,
                        "opposite_parity": opposite, "source_A_s-1": source_rate,
                        "source_frequency_kind": origin,
                        "source_frequency_cm-1": source_frequency,
                        "inferred_frequency_cm-1": inferred,
                        "inferred_frequency_relative_difference": frequency_error,
                        "physical_gap_cm-1": gap, "base_A_s-1": base})

    if stage == 3:
        for name in rates:
            for pair, absolute_rate in SZ_OPTICAL.items():
                rates[name][pair] = absolute_rate

    artifacts = {}
    base_dir = OUTPUT / f"O{stage}_base"
    base_dir.mkdir(parents=True, exist_ok=True)
    writer = write_excitation_atom("o", stage, output_dir=base_dir)
    source_lines = Path(writer["path"]).read_text().splitlines()
    for name, radiative in rates.items():
        destination = OUTPUT / f"O{stage}_{name}"
        destination.mkdir(parents=True, exist_ok=True)
        target = destination / f"nlevel_o_{stage}.txt"
        lines = source_lines.copy()
        marker = next(i for i, line in enumerate(lines) if line.startswith("NRAD "))
        count = int(lines[marker].split()[1])
        if count != len(radiative):
            raise ValueError(f"O {stage} changed radiative pair count: {count} != {len(radiative)}")
        source_pairs = {tuple(map(int, line.split()[:2])) for line in lines[marker+1:marker+1+count]}
        if source_pairs != set(radiative):
            raise ValueError(f"O {stage} changed radiative pair identities")
        replacements = [f"{l} {u} {radiative[l,u]:.17g}" for l,u in sorted(radiative)]
        lines[marker+1:marker+1+count] = replacements
        header = [f"# Private O {stage} {name} radiative candidate; no production installation",
                  "# Opposite-parity E1 strengths transformed to the selected physical photon gap",
                  "# Same-parity rates retained as estimates or conditionally scaled for sensitivity",
                  "# O III ground optical A: Storey and Zeippen (2000) Table 3 absolute totals"]
        target.write_text("\n".join(header + lines) + "\n")
        artifacts[name] = {"path": str(target), "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
                           "nlev": writer["nlev"], "nrad": writer["nrad"],
                           "nups": writer["nups"], "nodes": writer["nodes"]}
    result[f"O{stage}"] = {"artifacts": artifacts,
                           "opposite_parity_rows": sum(r["opposite_parity"] for r in records),
                           "reversed_opposite_parity_rows": sum(r["opposite_parity"] and r["physical_upper"] == r["source_lower"] for r in records),
                           "max_source_frequency_identification_error": max(r["inferred_frequency_relative_difference"] for r in records if r["opposite_parity"]),
                           "rows": records}
(HERE / "candidate_manifest.json").write_text(json.dumps({"passed": True, "ions": result}, indent=2) + "\n")
print(json.dumps({ion: {k: v for k, v in report.items() if k != "rows"}
                  for ion, report in result.items()}, indent=2))
