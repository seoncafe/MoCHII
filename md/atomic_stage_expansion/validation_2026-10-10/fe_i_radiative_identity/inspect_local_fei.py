"""Measure the local Fe I radiation inventory against the final article excerpt."""

import hashlib
import json
from pathlib import Path
import re
import subprocess


STOUT = Path.home() / "CLOUDY/c25.00/data/stout/fe/fe_1/fe_1"
PAPER = Path.home() / "MoCHII/references/Conroy_2020_ApJ_902_68.pdf"
RY_CM1 = 109737.3156816


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def levels(path):
    result = {}
    for line in path.read_text().splitlines()[1:]:
        match = re.match(r'^(\d+)\s+([\d.]+)\s+(\d+)\s+"(.*)"', line)
        if match:
            index, energy, weight, label = match.groups()
            result[int(index)] = (float(energy), int(weight), label)
    if sorted(result) != list(range(1, 847)):
        raise ValueError("Unexpected Stout Fe I level identities")
    return result


def radiation(path, energies):
    rates = {}
    records = 0
    multipole_sums = 0
    for line in path.read_text().splitlines()[1:]:
        fields = line.split()
        if not fields or fields[0] != "A":
            continue
        if len(fields) < 4:
            raise ValueError("Malformed Stout A record")
        first, second = int(fields[1]), int(fields[2])
        rate = float(fields[3])
        if first not in energies or second not in energies:
            raise ValueError("Out-of-bounds Stout A record")
        if energies[first][0] >= energies[second][0] or rate <= 0:
            raise ValueError("Stout A direction, energy, or rate is invalid")
        if (first, second) in rates:
            multipole_sums += 1
        rates[(first, second)] = rates.get((first, second), 0.0) + rate
        records += 1
    return rates, records, multipole_sums


def paper_excerpt(path):
    content = subprocess.check_output(["pdftotext", "-layout", str(path), "-"], text=True)
    start = content.index("Index       Level")
    table = content[start:content.index("portion of Figure 2", start)]
    result = {}
    for line in table.splitlines():
        match = re.match(r'^\s*(\d{1,2})\s+(.+?)\s+(\d\.\d{6})\s+(\d\.\d{6})\s+', line)
        if match and 1 <= int(match.group(1)) <= 50:
            index, label, _, observed = match.groups()
            result[int(index)] = (label, float(observed))
    # Table 2 wraps row 25 around Equation (2); the value is printed in its continuation.
    if "25                         0.175326         0.181456" not in table:
        raise ValueError("Published Table 2 row 25 changed")
    result[25] = ("3d64s4p 7D2", 0.181456)
    if sorted(result) != list(range(1, 51)):
        raise ValueError("Published Table 2 excerpt does not contain 50 rows")
    return result


def coarse_configuration(label):
    compact = re.sub(r'[^a-zA-Z0-9]', '', label)
    if "3d6" in compact and "4s2" in compact:
        return "3d6_4s2"
    if "3d7" in compact and "4s" in compact:
        return "3d7_4s"
    if "3d6" in compact and "4s" in compact and "4p" in compact:
        return "3d6_4s4p"
    raise ValueError(f"Unrecognized published low-state configuration: {label}")


def term_and_j(label, *, paper):
    compact = re.sub(r'[^A-Za-z0-9]', '', label)
    terms = re.findall(r'([1357][SPDFGHIK])(\d)', compact)
    if not terms:
        raise ValueError(f"Unrecognized published low-state term: {label}")
    term, j = terms[-1]
    if paper:
        return term, int(j)
    # Stout sometimes appends a term sequence number after the J value.
    return term, int(j)


def main():
    paths = {ext: Path(f"{STOUT}.{ext}") for ext in ("nrg", "tp", "coll")}
    energies = levels(paths["nrg"])
    rates, rate_records, multipole_sums = radiation(paths["tp"], energies)
    excerpt = paper_excerpt(PAPER)
    residuals = {str(index): energies[index][0] - value[1] * RY_CM1
                 for index, value in excerpt.items()}
    identity_mismatches = []
    for index, (paper_label, _) in excerpt.items():
        stout_energy, stout_weight, stout_label = energies[index]
        paper_term, paper_j = term_and_j(paper_label, paper=True)
        stout_term, _ = term_and_j(stout_label, paper=False)
        if (coarse_configuration(paper_label) != coarse_configuration(stout_label)
                or paper_term != stout_term or stout_weight != 2 * paper_j + 1
                or ("4p" in paper_label) != ("*" in stout_label)
                or abs(stout_energy - excerpt[index][1] * RY_CM1) > 0.06):
            identity_mismatches.append(index)
    counts = {}
    for size in (10, 50, 150, 485, 846):
        selected = {(low, upper): value for (low, upper), value in rates.items()
                    if low <= size and upper <= size}
        emitting = {upper for _, upper in selected}
        counts[str(size)] = {
            "positive_spontaneous_pairs": len(selected),
            "emitting_levels": len(emitting),
            "radiatively_closed_excited_levels": [index for index in range(2, size + 1)
                                                    if index not in emitting],
        }
    coll_lines = paths["coll"].read_text().splitlines()
    colliders = {line.split()[1] for line in coll_lines if line.startswith("RATE ")}
    result = {
        "stout_source_sha256": {extension: sha256(path) for extension, path in paths.items()},
        "published_article_sha256": sha256(PAPER),
        "stout_energy_levels": len(energies),
        "stout_spontaneous_pairs": len(rates),
        "stout_spontaneous_records": rate_records,
        "stout_repeated_pair_multipole_records": multipole_sums,
        "stout_collision_colliders": sorted(colliders),
        "stout_electron_collision_rows": sum(line.startswith("RATE ELECTRON")
                                              for line in coll_lines),
        "published_first_50_observed_energy_max_abs_difference_cm1": max(
            abs(value) for value in residuals.values()),
        "published_first_50_configuration_term_parity_J_energy_mismatch_indices": identity_mismatches,
        "published_first_50_observed_energy_difference_cm1": residuals,
        "local_stout_radiation_within_first_levels": counts,
        "published_table3_stout_2_to_1_s1": rates[(1, 2)],
        "published_table3_stout_3_to_2_s1": rates[(2, 3)],
    }
    output = Path(__file__).with_name("local_inventory.json")
    if output.exists():
        raise FileExistsError("Preserve existing Fe I inventory")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
