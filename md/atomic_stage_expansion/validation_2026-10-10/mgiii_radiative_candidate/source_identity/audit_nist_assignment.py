"""Compare observed CHIANTI Mg III level energies with official NIST levels."""

import csv
import hashlib
import json
from pathlib import Path
import re


HERE = Path(__file__).resolve().parent
NIST = HERE / "mgiii_nist_levels_2026-10-10.tsv"
CHIANTI = Path.home() / "RT_Codes/CHIANTI/dbase/mg/mg_3/mg_3.elvlc"
URL = ("https://physics.nist.gov/cgi-bin/ASD/energy1.pl?"
       "spectrum=Mg%20III&units=0&format=3&output=0&page_size=500&"
       "multiplet_ordered=0&conf_out=on&term_out=on&level_out=on&j_out=on")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def orbit(configuration):
    match = re.search(r"(\d+)([spdfg])$", configuration)
    return match.group(0) if match else None


def main():
    nist = []
    with NIST.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                nist.append({"configuration": row["Configuration"],
                             "term": row["Term"],
                             "J": float(row["J"]),
                             "energy_cm1": float(row["Level (cm-1)"])})
            except (ValueError, TypeError):
                continue
    chianti = []
    for line in CHIANTI.read_text().splitlines():
        fields = line.split()
        if len(fields) < 9 or not fields[0].isdigit():
            continue
        index = int(fields[0])
        if index <= 27:
            continue
        observed = float(fields[-2])
        if observed <= 0:
            continue
        chianti.append({"index": index, "configuration": ".".join(fields[1:4]),
                        "term": fields[4] + fields[5], "J": float(fields[6]),
                        "observed_cm1": observed,
                        "theoretical_cm1": float(fields[-1])})
    if len(nist) != 111 or len(chianti) != 68:
        raise ValueError("Local level-source counts changed")
    matches = []
    for level in chianti:
        same_J = [row for row in nist if row["J"] == level["J"]]
        nearest = min(same_J, key=lambda row:
                      abs(row["energy_cm1"] - level["observed_cm1"]))
        at_energy = min(nist, key=lambda row:
                        abs(row["energy_cm1"] - level["observed_cm1"]))
        matches.append({**level,
                        "same_J_nearest": {**nearest,
                          "difference_cm1": abs(nearest["energy_cm1"] -
                                                level["observed_cm1"]),
                          "orbital_match": orbit(nearest["configuration"]) ==
                                           orbit(level["configuration"])},
                        "energy_nearest_any_J": {**at_energy,
                          "difference_cm1": abs(at_energy["energy_cm1"] -
                                                level["observed_cm1"])}})
    close = [row for row in matches if row["same_J_nearest"]["difference_cm1"] <= 1]
    wrong_orbit = [row for row in close if not row["same_J_nearest"]["orbital_match"]]
    mismatch = [row for row in matches if row not in close]
    result = {
        "source": {"retrieval_date": "2026-10-10",
                   "url": URL, "nist_sha256": digest(NIST),
                   "chianti_sha256": digest(CHIANTI)},
        "nist_level_rows": len(nist),
        "chianti_higher_observed_rows": len(chianti),
        "same_J_energy_within_1_cm1": len(close),
        "same_J_energy_and_orbital_match_within_1_cm1": len(close)-len(wrong_orbit),
        "close_energy_wrong_orbit": wrong_orbit,
        "not_matching_same_J_energy_within_1_cm1": mismatch,
        "level_127_143_note": (
            "CHIANTI 127 has no observed energy; selected 620693 cm^-1 is "
            "theoretical. CHIANTI 143 observed 620701.688 cm^-1 matches "
            "NIST 6d J=1 at 620702 cm^-1. Their 8.688 cm^-1 selected gap "
            "is not independently verified by NIST for level 127."),
        "level_167_observed_assignment_valid": False,
        "physical_scope": (
            "Energy and J matching corroborates labels only where configuration "
            "also matches. It does not verify wavefunction mixing or E1 strengths."),
    }
    (HERE / "result.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
