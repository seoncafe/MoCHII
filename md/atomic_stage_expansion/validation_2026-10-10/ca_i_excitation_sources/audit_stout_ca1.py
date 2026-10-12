"""Count the local Ca I energy, radiation, and collision records."""

import hashlib
import json
from pathlib import Path


source = Path("/home/kiseon/CLOUDY/c25.00/data/stout/ca/ca_1")
levels = {}
for row in (source / "ca_1.nrg").read_text().splitlines():
    fields = row.split()
    if len(fields) >= 4 and fields[0].isdigit() and fields[1].replace(".", "", 1).isdigit():
        levels[int(fields[0])] = {
            "energy_cm1": float(fields[1]),
            "g": int(fields[2]),
            "label": row.split('"', 1)[1].rsplit('"', 1)[0],
        }

radiation = []
for row in (source / "ca_1.tp").read_text().splitlines():
    fields = row.split()
    if len(fields) >= 4 and fields[0] == "A":
        radiation.append((int(fields[1]), int(fields[2]), float(fields[3])))

collisions = []
for row in (source / "ca_1.coll").read_text().splitlines():
    fields = row.split()
    if fields and fields[0] in {"CE", "CS", "RATE", "BT"}:
        collisions.append(row)

decays = {}
for lower, upper, rate in radiation:
    assert lower in levels and upper in levels and rate > 0
    decays.setdefault(upper, set()).add(lower)


def reaches_ground(level, visited=None):
    if level == 1:
        return True
    if visited is None:
        visited = set()
    if level in visited:
        return False
    visited.add(level)
    return any(reaches_ground(lower, visited.copy()) for lower in decays.get(level, ()))


report = {
    "source_files": {
        name: {
            "path": str(source / name),
            "sha256": hashlib.sha256((source / name).read_bytes()).hexdigest(),
        }
        for name in ("ca_1.nrg", "ca_1.tp", "ca_1.coll")
    },
    "energy_level_count": len(levels),
    "radiative_A_record_count": len(radiation),
    "radiative_A_records_with_both_indices_at_most_39": sum(
        max(lower, upper) <= 39 for lower, upper, _ in radiation
    ),
    "electron_collision_record_count": len(collisions),
    "levels_2_to_39_without_radiative_exit": [
        level for level in range(2, 40) if level not in decays
    ],
    "levels_2_to_39_without_radiative_path_to_ground": [
        level for level in range(2, 40) if not reaches_ground(level)
    ],
    "first_nine_levels": {str(level): levels[level] for level in range(1, 10)},
}

output = Path(__file__).with_name("stout_ca1_inventory.json")
output.write_text(json.dumps(report, indent=2) + "\n")
print(output)
