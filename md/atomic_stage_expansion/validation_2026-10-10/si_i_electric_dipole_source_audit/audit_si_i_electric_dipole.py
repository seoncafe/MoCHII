"""Compare original NIST MCHF Si I records with the local CHIANTI radiation table.

The NIST .lsj file is external scientific source preserved byte-for-byte.
This script reads it and the unchanged CHIANTI files and writes only the
specified diagnostic JSON. It does not select replacement probabilities.
"""

import argparse
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path


NIST_SHA256 = "438150a064b12c202433f590122da7eaf819287e4656d0a6c5ba88586744d904"
CHI_ELVLC_SHA256 = "1538e98411f8c2f03414aebd36f901f7b88bc79109b1f0aea7254af070bd4f57"
CHI_WGFA_SHA256 = "be54f54d2e91695b766c0cb545e744d0d540165add1231775b7361f82ab21c42"
CHI_DIR = Path("/home/kiseon/RT_Codes/CHIANTI/dbase/si/si_1")
DEFAULT_NIST = Path(__file__).with_name("Si_14.16.cff-mchf.lsj.db")
DEFAULT_DISPLAY = Path(__file__).with_name("Si_14.16.cff-mchf-lin.dat.mp")
DISPLAY_SHA256 = "049b438f9513a045f9cc2c36dba052f88734a2710d2015da88174b6713b592ee"
LEVEL_RE = re.compile(r"^\s*(\d+)\s+(-\d+\.\d+)\s+(\S+)\s*$")
ENERGY_RE = re.compile(r"^\s*([\d.]+) CM-1\s+([\d.]+) ANGS\(VAC\)")
RATE_RE = re.compile(
    r"^\s*(E1|E2|M1|M2)\s+S\s*=\s*([\d.eE+-]+)\s+"
    r"GF\s*=\s*([\d.eE+-]+)\s+AKI\s*=\s*([\d.eE+-]+)"
)
VELOCITY_RE = re.compile(r"^\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s*$")
ORBITAL_RE = re.compile(r"(\d[spdfg])(?:(\d+)|\((\d+)\))?")
DISPLAY_E1_RE = re.compile(
    r"\b(\d+)\s+(\d+)\s+E1\s+([\d.]+)\s+([\d.]+)\s+"
    r"([\d.]+)\s+[\d.]+\s+([\d.eE+-]+)"
)
ORBITAL_ANGULAR_MOMENTUM = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def level_identity(configuration, term, twice_j):
    orbitals = ORBITAL_RE.findall(configuration)
    if not orbitals:
        raise ValueError(f"Unparsed orbital configuration: {configuration}")
    return (tuple((orbital, int(bare_count or parenthesized_count or 1))
                  for orbital, bare_count, parenthesized_count in orbitals),
            term.rstrip("*"), int(twice_j))


def nist_identity(label, twice_j):
    configuration, term = label.rsplit("_", 1)
    return level_identity(configuration, term, twice_j)


def chianti_levels(path):
    levels = {}
    identities = {}
    for line in path.read_text().splitlines():
        if line.lstrip().startswith("-1"):
            break
        fields = line.split()
        if not fields:
            continue
        index = int(fields[0])
        configuration = " ".join(fields[1:-5])
        term = fields[-5] + fields[-4]
        twice_j = round(2 * float(fields[-3]))
        identity = level_identity(configuration, term, twice_j)
        if identity in identities:
            raise ValueError(f"Ambiguous CHIANTI identity: {identity}")
        identities[identity] = index
        levels[index] = {"identity": identity, "observed_cm_inverse": float(fields[-2]),
                         "theoretical_cm_inverse": float(fields[-1])}
    return levels, identities


def nist_records(path, identities):
    lines = path.read_text().splitlines()
    records = []
    unmapped = []
    for line_number, line in enumerate(lines):
        rate = RATE_RE.match(line)
        if rate is None:
            continue
        if line_number < 3:
            raise ValueError("Incomplete NIST transition record")
        first = LEVEL_RE.match(lines[line_number - 3])
        second = LEVEL_RE.match(lines[line_number - 2])
        energy = ENERGY_RE.match(lines[line_number - 1])
        if first is None or second is None or energy is None:
            raise ValueError(f"Unexpected NIST record at line {line_number + 1}")
        first_identity = nist_identity(first[3], int(first[1]))
        second_identity = nist_identity(second[3], int(second[1]))
        first_id, second_id = identities.get(first_identity), identities.get(second_identity)
        record = {
            "nist_line": line_number + 1,
            "first_label": first[3], "second_label": second[3],
            "first_chianti_id": first_id, "second_chianti_id": second_id,
            "first_g": int(first[1]) + 1, "second_g": int(second[1]) + 1,
            "first_total_energy_hartree": float(first[2]),
            "second_total_energy_hartree": float(second[2]),
            "multipole": rate[1], "source_gap_cm_inverse": float(energy[1]),
            "source_wavelength_angstrom": float(energy[2]),
            "length_S_atomic_units": float(rate[2]),
            "length_GF": float(rate[3]), "length_A_per_second": float(rate[4]),
        }
        if line_number + 1 < len(lines):
            velocity = VELOCITY_RE.match(lines[line_number + 1])
            if velocity is not None:
                record.update({"velocity_S_atomic_units": float(velocity[1]),
                               "velocity_GF": float(velocity[2]),
                               "velocity_A_per_second": float(velocity[3])})
        records.append(record)
        if first_id is None or second_id is None:
            unmapped.append(record)
    return records, unmapped


def display_multiplicity(path, records):
    display = []
    for line in path.read_text().splitlines():
        match = DISPLAY_E1_RE.search(line)
        if match is not None:
            display.append((int(match[1]), int(match[2]), float(match[5]), float(match[6])))
    for record in records:
        if record["multipole"] != "E1":
            continue
        record["formatted_display_count"] = sum(
            gi == record["first_g"] and gk == record["second_g"]
            and abs(gap - record["source_gap_cm_inverse"]) <= 0.2
            and abs(strength / record["length_S_atomic_units"] - 1) <= 0.001
            for gi, gk, gap, strength in display)


def parity(identity):
    return sum(ORBITAL_ANGULAR_MOMENTUM[orbital[-1]] * count
               for orbital, count in identity[0]) % 2


def chianti_radiation(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.lstrip().startswith("-1"):
            break
        fields = line.split()
        if not fields:
            continue
        first, second = map(int, fields[:2])
        rows.append({"first_id": first, "second_id": second,
                     "wavelength_angstrom": float(fields[2]),
                     "gf": float(fields[3]), "A_per_second": float(fields[4])})
    return rows


def compare(nist_path, display_path, chi_dir):
    elvlc_path = chi_dir / "si_1.elvlc"
    wgfa_path = chi_dir / "si_1.wgfa"
    hashes = {"nist_original_lsj": sha256(nist_path),
              "nist_formatted_transition_table": sha256(display_path),
              "chianti_elvlc": sha256(elvlc_path),
              "chianti_wgfa": sha256(wgfa_path)}
    expected = {"nist_original_lsj": NIST_SHA256,
                "nist_formatted_transition_table": DISPLAY_SHA256,
                "chianti_elvlc": CHI_ELVLC_SHA256,
                "chianti_wgfa": CHI_WGFA_SHA256}
    if hashes != expected:
        raise ValueError(f"Source hash mismatch: {hashes}")
    levels, identities = chianti_levels(elvlc_path)
    nist, unmapped = nist_records(nist_path, identities)
    display_multiplicity(display_path, nist)
    radiation = chianti_radiation(wgfa_path)
    keyed = defaultdict(list)
    for record in nist:
        a, b = record["first_chianti_id"], record["second_chianti_id"]
        if a is not None and b is not None:
            keyed[tuple(sorted((a, b)))].append(record)
    comparison = []
    for row in radiation:
        if row["gf"] <= 0 or row["A_per_second"] <= 0:
            continue
        a, b = row["first_id"], row["second_id"]
        matches = keyed[tuple(sorted((a, b)))]
        positive = [r for r in matches if r["multipole"] == "E1"]
        if len(positive) > 1:
            raise ValueError(f"Multiple original E1 records for CHIANTI pair {a},{b}")
        physical_lower, physical_upper = sorted((a, b), key=lambda i: levels[i]["observed_cm_inverse"])
        opposite_parity = parity(levels[a]["identity"]) != parity(levels[b]["identity"])
        delta_j = abs(levels[a]["identity"][2] - levels[b]["identity"][2]) / 2
        e1_allowed = (opposite_parity and delta_j <= 1 and
                      not (levels[a]["identity"][2] == 0 and levels[b]["identity"][2] == 0))
        entry = {"chianti": row, "physical_lower_id": physical_lower,
                 "physical_upper_id": physical_upper,
                 "chianti_orientation_is_physical": (a, b) == (physical_lower, physical_upper),
                 "opposite_parity": opposite_parity, "delta_J": delta_j,
                 "E1_selection_allowed": e1_allowed,
                 "matched_multipoles": [r["multipole"] for r in matches],
                 "nist_E1": positive[0] if positive else None}
        if positive:
            if not e1_allowed:
                raise ValueError(f"NIST E1 selection inconsistent with CHIANTI pair {a},{b}")
            source = positive[0]
            entry["chianti_to_nist_GF_ratio"] = row["gf"] / source["length_GF"]
            entry["chianti_to_nist_A_ratio"] = row["A_per_second"] / source["length_A_per_second"]
            gap = abs(levels[a]["observed_cm_inverse"] - levels[b]["observed_cm_inverse"])
            emitted_g = levels[physical_upper]["identity"][2] + 1
            physical_lambda = 1e8 / gap
            entry["physical_E1_A_from_nist_length_S_per_second"] = (
                2.02613e18 * source["length_S_atomic_units"] /
                (emitted_g * physical_lambda**3))
            if source.get("velocity_S_atomic_units") is not None:
                entry["conditional_physical_E1_A_from_nist_velocity_S_per_second"] = (
                    2.02613e18 * source["velocity_S_atomic_units"] /
                    (emitted_g * physical_lambda**3))
            entry["selected_observed_gap_cm_inverse"] = gap
            entry["nist_order_inverts_observed_order"] = (
                source["first_chianti_id"] != physical_lower)
            entry["nist_order_inverts_chianti_theoretical_order"] = (
                source["first_chianti_id"] != min((a, b),
                    key=lambda i: levels[i]["theoretical_cm_inverse"]))
        comparison.append(entry)
    matched = [row for row in comparison if row["nist_E1"] is not None]
    unmatched = [row for row in comparison if row["nist_E1"] is None]
    bound = [row for row in comparison if max(row["chianti"]["first_id"],
                                              row["chianti"]["second_id"]) <= 35]
    doubled = [row for row in matched if 1.9 < row["chianti_to_nist_GF_ratio"] < 2.1]
    if any(row["E1_selection_allowed"] for row in bound if row["nist_E1"] is None):
        raise ValueError("E1-allowed bound radiation has no matched NIST E1 record")
    if any(row["nist_E1"]["formatted_display_count"] not in (1, 2) for row in matched):
        raise ValueError("Unidentified formatted NIST E1 display multiplicity")
    if any((row["nist_E1"]["formatted_display_count"] == 2) !=
           (1.9 < row["chianti_to_nist_GF_ratio"] < 2.1) for row in matched):
        raise ValueError("CHIANTI E1 ratio and NIST display multiplicity disagree")
    chianti_pairs = {tuple(sorted((row["first_id"], row["second_id"])))
                     for row in radiation}
    nist_e1_absent_from_chianti = [
        record for record in nist if record["multipole"] == "E1" and
        tuple(sorted((record["first_chianti_id"],
                      record["second_chianti_id"]))) not in chianti_pairs]
    source_states = set()
    for record in nist:
        source_states.update((record["first_chianti_id"],
                              record["second_chianti_id"]))
    gauge = [record["length_S_atomic_units"] / record["velocity_S_atomic_units"]
             for record in nist if record["multipole"] == "E1" and
             record.get("velocity_S_atomic_units", 0) > 0]
    selected_gauge = [row["nist_E1"]["length_S_atomic_units"] /
                      row["nist_E1"]["velocity_S_atomic_units"]
                      for row in matched if row["nist_E1"].get("velocity_S_atomic_units", 0) > 0]
    bound_without_nist_E1 = [row for row in bound if row["nist_E1"] is None]
    state_table = [{"chianti_id": index,
                    "configuration": " ".join(orbital + (str(count) if count > 1 else "")
                                               for orbital, count in level["identity"][0]),
                    "term": level["identity"][1],
                    "J": level["identity"][2] / 2,
                    "g": level["identity"][2] + 1,
                    "parity": "odd" if parity(level["identity"]) else "even",
                    "observed_cm_inverse": level["observed_cm_inverse"],
                    "chianti_theoretical_cm_inverse": level["theoretical_cm_inverse"],
                    "appears_in_nist_original": index in source_states}
                   for index, level in sorted(levels.items())]
    return {"input_sha256": hashes, "chianti_levels": len(levels),
            "state_table": state_table,
            "states_appearing_in_nist_original": len(source_states),
            "nist_records": len(nist), "nist_unmapped_records": len(unmapped),
            "nist_unmapped_examples": unmapped[:10],
            "nist_multipole_counts": {multipole: sum(record["multipole"] == multipole
                                                for record in nist)
                                        for multipole in ("E1", "E2", "M1", "M2")},
            "nist_E1_records_absent_from_chianti": nist_e1_absent_from_chianti,
            "chianti_positive_gf_A_rows": len(comparison),
            "bound_positive_gf_A_rows": len(bound),
            "matched_nist_E1_rows": len(matched), "unmatched_nist_E1_rows": len(unmatched),
            "bound_matched_nist_E1_rows": sum(row["nist_E1"] is not None for row in bound),
            "bound_same_parity_rows_without_nist_E1": sum(
                row["nist_E1"] is None and not row["opposite_parity"] for row in bound),
            "bound_E1_forbidden_opposite_parity_rows_without_nist_E1": sum(
                row["nist_E1"] is None and row["opposite_parity"] for row in bound),
            "bound_non_E1_rows_all_have_physical_A_orientation": all(
                row["chianti_orientation_is_physical"] for row in bound_without_nist_E1),
            "bound_non_E1_rows_without_nist_match": [
                [row["chianti"]["first_id"], row["chianti"]["second_id"]]
                for row in bound_without_nist_E1 if not row["matched_multipoles"]],
            "nist_E1_order_inversions_vs_observed": [
                [row["chianti"]["first_id"], row["chianti"]["second_id"]]
                for row in matched if row["nist_order_inverts_observed_order"]],
            "nist_E1_order_inversions_vs_chianti_theoretical": [
                [row["chianti"]["first_id"], row["chianti"]["second_id"]]
                for row in matched if row["nist_order_inverts_chianti_theoretical_order"]],
            "nist_E1_length_to_velocity_S_ratio_range": [min(gauge), max(gauge)],
            "matched_E1_length_to_velocity_S_ratio_range": [min(selected_gauge),
                                                             max(selected_gauge)],
            "roughly_doubled_nist_GF_pairs": [
                [row["chianti"]["first_id"], row["chianti"]["second_id"]]
                for row in doubled],
            "roughly_doubled_nist_GF_count": len(doubled),
            "roughly_doubled_nist_GF_ratio_range": [
                min(row["chianti_to_nist_GF_ratio"] for row in doubled),
                max(row["chianti_to_nist_GF_ratio"] for row in doubled)],
            "unique_display_nist_GF_max_fractional_difference": max(
                abs(row["chianti_to_nist_GF_ratio"] - 1)
                for row in matched if row["nist_E1"]["formatted_display_count"] == 1),
            "doubled_pairs_with_two_formatted_display_rows": sum(
                row["nist_E1"]["formatted_display_count"] == 2 for row in doubled),
            "unmatched_pairs": [[r["chianti"]["first_id"], r["chianti"]["second_id"]]
                                for r in unmatched],
            "rows": comparison}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nist", type=Path, default=DEFAULT_NIST)
    parser.add_argument("--display", type=Path, default=DEFAULT_DISPLAY)
    parser.add_argument("--chianti-dir", type=Path, default=CHI_DIR)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("result.json"))
    args = parser.parse_args()
    result = compare(args.nist, args.display, args.chianti_dir)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print({key: result[key] for key in ("nist_records", "nist_unmapped_records",
                                      "chianti_positive_gf_A_rows", "matched_nist_E1_rows",
                                      "unmatched_nist_E1_rows")})


if __name__ == "__main__":
    main()
