"""Write and structurally check Ne I with two measured Stout decay rates."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OLD = HERE.parent / "source_validity_ne_i_mg_i/ne_i_candidate/nlevel_ne_1.txt"
STOUT = Path.home() / "CLOUDY/c25.00/data/stout/ne/ne_1/ne_1"
sys.path.insert(0, str(ROOT / "tools/fitting"))
import write_excitation_atom as W  # noqa: E402


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
    mark += nlev + 1
    nrad = int(lines[mark].split()[1])
    radiation = {}
    for line in lines[mark + 1:mark + 1 + nrad]:
        lower, upper, value = line.split()
        pair = (int(lower), int(upper))
        if pair in radiation:
            raise ValueError(f"Duplicate radiative pair {pair}")
        radiation[pair] = float(value)
    mark += nrad + 1
    assert lines[mark] == "NODES 6975"
    assert lines[mark + 1] == "NUPS 465"
    mark += 2
    collisions = {}
    for _ in range(465):
        fields = lines[mark].split()
        count = int(fields[-1])
        pair = tuple(sorted((int(fields[0]), int(fields[1]))))
        nodes = tuple(tuple(map(float, line.split()))
                      for line in lines[mark + 1:mark + count + 1])
        if pair in collisions:
            raise ValueError(f"Duplicate collision pair {pair}")
        collisions[pair] = (tuple(fields[2:-1]), nodes)
        mark += count + 1
    floor = [line for line in lines if line.startswith("BT_TMIN_K ")]
    return levels, radiation, collisions, floor


def source_rates():
    rows = {}
    for line in Path(str(STOUT) + ".tp").read_text().splitlines():
        fields = line.split()
        if len(fields) == 4 and fields[0] == "A":
            pair = (int(fields[1]), int(fields[2]))
            rows[pair] = float(fields[3])
    return rows


def main():
    output = HERE / "atom/nlevel_ne_1.txt"
    result_path = HERE / "candidate_result.json"
    if output.exists() or result_path.exists():
        raise FileExistsError("Preserve existing corrected Ne I candidate")
    old_hash = digest(OLD)
    source_paths = W._source_files("ne", 1)
    additional_sources = [Path(str(STOUT) + suffix) for suffix in (".nrg", ".tp")]
    source_hashes = {str(path): digest(path) for path in
                     [*source_paths, *additional_sources]}
    metadata = W.write_excitation_atom("ne", 1, levels=31,
                                       output_dir=HERE / "atom")
    if source_hashes != {str(path): digest(path) for path in
                         [*source_paths, *additional_sources]}:
        raise ValueError("Raw source changed during candidate generation")
    if (not metadata["passed"] or (metadata["nlev"], metadata["nrad"],
            metadata["nups"], metadata["nodes"]) != (31, 144, 465, 6975)
            or metadata["selected_ids"] != list(range(1, 32))
            or metadata["bt_tmin_K"] != 1000.36):
        raise ValueError("Corrected Ne I source counts or floor changed")
    old_levels, old_radiation, old_collision, old_floor = records(OLD)
    levels, radiation, collision, floor = records(output)
    if old_levels != levels or old_collision != collision:
        raise ValueError("Existing level identity or collision source changed")
    if any(radiation.get(pair) != rate for pair, rate in old_radiation.items()):
        raise ValueError("An original CHIANTI radiative rate changed")
    added = {pair: rate for pair, rate in radiation.items()
             if pair not in old_radiation}
    expected_added = {(3, 15): 9.00e5, (5, 15): 6.82e7}
    if added != expected_added or any(pair in old_radiation for pair in expected_added):
        raise ValueError("Expected two Stout transitions were not added exactly")
    stout = source_rates()
    if any(stout.get(pair) != rate for pair, rate in expected_added.items()):
        raise ValueError("Stout source A values changed")
    if floor != ["BT_TMIN_K 1000.36"]:
        raise ValueError("Ne I inferred BT floor record changed")
    if old_floor:
        raise ValueError("Historical Ne I atom unexpectedly had BT metadata")
    result = {
        "passed": True,
        "candidate_path": str(output.resolve()),
        "candidate_sha256": digest(output),
        "historical_atom_sha256": old_hash,
        "stout_added_radiative_pairs_checked": True,
        "added_pairs_s-1": {f"{l},{u}": rate for (l, u), rate in added.items()},
        "historical_radiative_rows": len(old_radiation),
        "nlev": 31, "nrad": 144, "nups": 465, "nodes": 6975,
        "bt_tmin_K": 1000.36,
        "raw_source_sha256": source_hashes,
        "code_sha256": {str(path): digest(path) for path in
                        (ROOT / "tools/fitting/chianti_cooling.py",
                         ROOT / "tools/fitting/write_excitation_atom.py",
                         ROOT / "tools/fitting/atomic_excitation_energies.py")},
        "previous_atom_physically_stale": True,
    }
    result_path.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
