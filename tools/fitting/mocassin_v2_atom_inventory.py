"""Inspect the atoms actually selected by the MOCASSIN v2 data reader.

This is a read-only source inventory. Missing files are unavailable line
models, not missing charge states. External files remain unchanged.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


ELEMENTS = (("c", 6), ("n", 7), ("o", 8), ("ne", 10), ("s", 16),
            ("ar", 18), ("mg", 12), ("fe", 26), ("si", 14),
            ("cl", 17), ("ca", 20))
MOCASSIN = Path.home() / "RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def number(value):
    return float(value.rstrip(",").replace("D", "E").replace("d", "e"))


def read_atom(path, include_tables=False):
    """Follow hydro_mod.makeElements, including its two table formats.

    emission_mod.equilibrium applies the Upsilon-to-rate formula in both
    formats; the historical iRats comment is not the invoked behavior.
    """
    lines = path.read_text().splitlines()
    comments = int(lines[0])
    cursor = comments + 1
    nlev, ntemp = map(int, lines[cursor].split()[:2])
    cursor += 1
    labels = lines[cursor:cursor + nlev]
    cursor += nlev
    temperatures = [number(line.split()[0]) for line in lines[cursor:cursor + ntemp]]
    cursor += ntemp
    mode = number(lines[cursor].split()[0])
    cursor += 1
    collision = {}
    if mode == 0:
        previous_lower = previous_upper = None
        while True:
            lower, upper, strength = lines[cursor].split()[:3]
            cursor += 1
            lower, upper, strength = int(lower), int(upper), number(strength)
            if strength == 0:
                break
            if lower == 0:
                lower = previous_lower
            if upper == 0:
                upper = previous_upper
            if lower is None or upper is None:
                raise ValueError(f"Unresolved collision continuation: {path}")
            previous_lower, previous_upper = lower, upper
            collision.setdefault((lower, upper), []).append(strength)
        radiation = []
        for _ in range(nlev * (nlev - 1) // 2):
            lower, upper, rate = lines[cursor].split()[:3]
            cursor += 1
            radiation.append((int(lower), int(upper), number(rate)))
    elif mode == 1:
        while True:
            fields = lines[cursor].split()
            cursor += 1
            while len(fields) < 2 + ntemp:
                fields.extend(lines[cursor].split())
                cursor += 1
            lower, upper = map(int, fields[:2])
            values = list(map(number, fields[2:2 + ntemp]))
            if len(values) != ntemp:
                raise ValueError(f"Collision-rate row length: {path}")
            collision[(lower, upper)] = values
            if lower >= nlev - 1:
                break
        radiation = []
        for _ in range(2502):
            fields = lines[cursor].split()
            cursor += 1
            if len(fields) != 12:
                raise ValueError(f"Packed Fe II radiation row: {path}")
            for offset in range(0, 12, 3):
                # The file stores lower,upper; the v2 packed reader stores
                # A(lower,upper), unlike the emitting-first consumer.
                radiation.append((int(fields[offset]), int(fields[offset + 1]),
                                  number(fields[offset + 2])))
        fields = lines[cursor].split()
        cursor += 1
        if len(fields) != 9:
            raise ValueError(f"Final packed Fe II radiation row: {path}")
        for offset in range(0, 9, 3):
            radiation.append((int(fields[offset]), int(fields[offset + 1]),
                              number(fields[offset + 2])))
    else:
        raise ValueError(f"Unsupported reader mode {mode}: {path}")
    levels = []
    for index in range(nlev):
        fields = lines[cursor].split()
        cursor += 1
        level, weight, energy = int(fields[0]), number(fields[1]), number(fields[2])
        levels.append(dict(index=level, label=labels[index].strip(),
                           statistical_weight=weight, energy_cm_inverse=energy))
    level_map = {row["index"]: row for row in levels}
    if sorted(level_map) != list(range(1, nlev + 1)):
        raise ValueError(f"Missing or duplicate level identity: {path}")
    if any(len(values) != ntemp for values in collision.values()):
        raise ValueError(f"Incomplete collision-temperature pair: {path}")
    if any(not 1 <= lower < upper <= nlev for lower, upper in collision):
        raise ValueError(f"Invalid collision pair: {path}")
    if any(not math.isfinite(value) or value < 0
           for values in collision.values() for value in values):
        raise ValueError(f"Invalid source collision coefficient: {path}")
    if any(not math.isfinite(rate) or rate < 0 for _, _, rate in radiation):
        raise ValueError(f"Invalid source spontaneous probability: {path}")
    positive_radiation = [row for row in radiation if row[2] > 0]
    inverted = [[lower, upper] for lower, upper, _ in positive_radiation
                if level_map[upper]["energy_cm_inverse"] <= level_map[lower]["energy_cm_inverse"]]
    trailing = [line.split() for line in lines[cursor:] if line.strip()]
    result = dict(nlev=nlev, temperature_points=ntemp, temperatures_K=temperatures,
                collision_format="effective_collision_strength", reader_format=int(mode),
                collision_pair_count=len(collision),
                collision_value_count=sum(map(len, collision.values())),
                positive_collision_pair_count=sum(any(v > 0 for v in row) for row in collision.values()),
                radiation_row_count=len(radiation),
                positive_radiation_pair_count=len(positive_radiation),
                reader_downward_positive_A_count=0 if mode == 1 else len(positive_radiation),
                packed_radiation_reader_defect=mode == 1,
                nonpositive_photon_energy_pairs=inverted, levels=levels,
                trailing_recombination_rows=len(trailing),
                trailing_recombination_field_counts=sorted(set(map(len, trailing))),
                source_comments=lines[1:comments + 1])
    if include_tables:
        result['collision_tables'] = [dict(lower=lower, upper=upper, upsilon=values)
                                      for (lower, upper), values in collision.items()]
        result['spontaneous_probabilities'] = [dict(lower=lower, upper=upper, A_s_inverse=rate)
                                             for lower, upper, rate in radiation]
    return result


def inventory(root=MOCASSIN):
    filenames = (root / "data/fileNames.dat").read_text().splitlines()
    selected = {}
    cursor = 0
    for charge in range(3, 31):
        for stage in range(1, min(charge + 1, 10) + 1):
            selected[(charge, stage)] = filenames[cursor][:20].strip()
            cursor += 1
    records = []
    for element, charge in ELEMENTS:
        for stage in range(1, 8):
            relative = selected[(charge, stage)]
            path = root / relative
            record = dict(element=element, Z=charge, stage=stage,
                          relative_file=relative, line_file_available=path.is_file(),
                          fully_stripped=stage == charge + 1)
            if path.is_file():
                record.update(source_sha256=digest(path), atom=read_atom(path))
            records.append(record)
    evidence = [root / "data/fileNames.dat", root / "source/set_input_mod.f90",
                root / "source/hydro_mod.f90", root / "source/emission_mod.f90",
                root / "source/update_mod.f90", root / "source/constants_mod.f90"]
    return dict(scope="Eleven selected metal elements; nstages=7; actual v2 file mapping and reader formats",
                source_root=str(root.resolve()), metal_stage_count=len(records),
                adjacent_transition_count=66, hydrogen_helium_stage_count=5,
                total_stage_count=82, available_line_atoms=sum(r["line_file_available"] for r in records),
                missing_line_files=sum(not r["line_file_available"] for r in records),
                source_sha256={str(path.resolve()): digest(path) for path in evidence}, ions=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, default=MOCASSIN)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Preserve the existing MOCASSIN source inventory")
    result = inventory(args.source_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print({key: result[key] for key in ("metal_stage_count", "adjacent_transition_count",
                                       "available_line_atoms", "missing_line_files")})


if __name__ == "__main__":
    main()
