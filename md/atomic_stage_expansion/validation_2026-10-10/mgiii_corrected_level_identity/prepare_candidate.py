"""Write and structurally verify the corrected-level Mg III source-node atom."""

from __future__ import annotations

from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OLD = HERE.parent / "mgiii_radiative_candidate"
sys.path.insert(0, str(ROOT / "tools/fitting"))
sys.path.insert(0, str(HERE))
import write_excitation_atom as W  # noqa: E402
from run_profiles import source_atom  # noqa: E402


OLD_ATOM = OLD / "nlevel_mg_3_line_strength_candidate.txt"
FIRST27 = OLD / "first27_original_A/nlevel_mg_3.txt"
RESULT = HERE / "candidate_result.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(path: Path):
    lines = path.read_text().splitlines()
    mark = lines.index("NLEVEL_VERSION 2") + 1
    nlev = int(lines[mark].split()[1])
    levels = {}
    for line in lines[mark + 1:mark + 1 + nlev]:
        index, weight, energy = line.split()
        levels[int(index)] = (float(weight), float(energy))
    mark += 1 + nlev
    nrad = int(lines[mark].split()[1])
    radiative = {}
    for line in lines[mark + 1:mark + 1 + nrad]:
        lower, upper, value = line.split()
        radiative[(int(lower), int(upper))] = float(value)
    mark += 1 + nrad
    assert lines[mark] == "NODES 207915"
    assert lines[mark + 1] == "NUPS 13861"
    mark += 2
    collision = {}
    for _ in range(13861):
        fields = lines[mark].split()
        lower, upper, count = int(fields[0]), int(fields[1]), int(fields[-1])
        pair = tuple(sorted((lower, upper)))
        nodes = tuple(tuple(map(float, line.split()))
                      for line in lines[mark + 1:mark + 1 + count])
        if pair in collision:
            raise ValueError(f"Duplicate collision pair: {pair}")
        collision[pair] = (tuple(fields[2:-1]), nodes)
        mark += count + 1
    assert lines[mark] == "BT_TMIN_K 1800"
    return levels, radiative, collision


def main() -> None:
    if RESULT.exists():
        raise FileExistsError(f"Preserve existing candidate result: {RESULT}")
    atom_path = HERE / "atom/nlevel_mg_3.txt"
    if atom_path.exists():
        raise FileExistsError(f"Preserve existing candidate atom: {atom_path}")
    if digest(OLD_ATOM) != "46e71f59c1c885a46ec4304f4a26bc1549e716a5490c6118d4cd242b16711873":
        raise ValueError("The preserved old candidate changed")
    if digest(FIRST27) != "ff9347360eeaa58bf04f008d89dbf344810f085fce3f7013100471df38ec8095":
        raise ValueError("The preserved first-27 atom changed")
    source_paths = W._source_files("mg", 3)
    source_hashes = {str(path): digest(path) for path in source_paths}
    metadata = W.write_excitation_atom("mg", 3, levels=167,
                                       output_dir=HERE / "atom")
    if source_hashes != {str(path): digest(path) for path in source_paths}:
        raise ValueError("Raw source changed while writing the candidate")
    if (not metadata["passed"] or (metadata["nlev"], metadata["nrad"],
            metadata["nups"], metadata["nodes"]) != (167, 3364, 13861, 207915)
            or metadata["selected_ids"] != list(range(1, 168))
            or metadata["bt_tmin_K"] != 1800):
        raise ValueError("Corrected full-atom metadata changed")
    old_levels, old_radiation, old_collision = records(OLD_ATOM)
    new_levels, new_radiation, new_collision = records(atom_path)
    if old_levels[167][1] != 628148.375 or new_levels[167][1] != 634430:
        raise ValueError("The ID 167 observed-energy defect was not corrected")
    if {key: value for key, value in old_levels.items() if key != 167} != {
            key: value for key, value in new_levels.items() if key != 167}:
        raise ValueError("A level other than ID 167 changed")
    if old_collision != new_collision:
        raise ValueError("Collision types, BT scaling, or source nodes changed")
    if len(new_radiation) != 3364:
        raise ValueError("The radiative row count changed")
    for pair, rate in new_radiation.items():
        if 167 not in pair and old_radiation.get(pair) != rate:
            raise ValueError(f"Unrelated radiative rate changed: {pair}")
    levels, direct_radiation, direct_collision = source_atom()
    if len(direct_collision) != 13861 or levels != new_levels:
        raise ValueError("Direct source-level identities differ")
    expected = defaultdict(float)
    for lower, upper, _, _, rate in direct_radiation:
        expected[(lower, upper)] += rate
    if set(expected) != set(new_radiation):
        raise ValueError("Serialized and direct radiative identities differ")
    max_error = max(abs(new_radiation[pair] / rate - 1)
                    for pair, rate in expected.items())
    if max_error > 5e-15:
        raise ValueError("Serialized corrected A differs from independent raw-source transform")
    changed_167 = [pair for pair in new_radiation
                   if 167 in pair and old_radiation.get(pair) != new_radiation[pair]]
    result = {
        "passed": True,
        "candidate_path": str(atom_path.resolve()),
        "candidate_sha256": digest(atom_path),
        "old_candidate_sha256": digest(OLD_ATOM),
        "first27_sha256": digest(FIRST27),
        "level_167_identity_checked": True,
        "old_level_167_cm1": old_levels[167][1],
        "selected_level_167_cm1": new_levels[167][1],
        "other_selected_levels_unchanged": True,
        "source_collision_records_unchanged": True,
        "first27_atom_unchanged": True,
        "changed_radiative_rows_involving_167": len(changed_167),
        "maximum_relative_serialized_vs_independent_A_error": max_error,
        "nlev": 167, "nrad": 3364, "nups": 13861,
        "nodes": 207915, "bt_tmin_K": 1800,
        "raw_source_sha256": source_hashes,
        "code_sha256": {str(path): digest(path) for path in
                        (ROOT / "tools/fitting/chianti_cooling.py",
                         ROOT / "tools/fitting/write_excitation_atom.py",
                         ROOT / "tools/fitting/electric_dipole_probabilities.py",
                         ROOT / "tools/fitting/atomic_excitation_energies.py")},
    }
    RESULT.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
