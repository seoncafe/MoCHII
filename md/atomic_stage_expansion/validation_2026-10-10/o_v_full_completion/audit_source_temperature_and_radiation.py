"""Measure the distributed O V source temperatures and photon-frequency metadata."""

import hashlib
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as C  # noqa: E402


CHIANTI = Path.home() / "RT_Codes/CHIANTI/dbase/o/o_5"
STOUT = Path.home() / "CLOUDY/c25.00/data/stout/o/o_5/o_5.coll"
PAPER = Path.home() / "MoCHII/references/Fernandez-Menchero_2014_AA_566_A104.pdf"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def first_finite_bt_temperature(row):
    abscissa = float(row["xs"][np.flatnonzero(row["xs"] > 0)[0]])
    scale = row["cups"]
    if row["ttype"] in (1, 4):
        scaled_temperature = scale ** (1 / (1 - abscissa)) - scale
    else:
        scaled_temperature = scale * abscissa / (1 - abscissa)
    return float(scaled_temperature * row["de"] / C.KB_OVER_RY)


def main():
    destination = HERE / "source_audit.json"
    if destination.exists():
        raise FileExistsError(destination)
    collisions = C.read_scups(CHIANTI / "o_5.scups")
    temperatures = np.asarray([first_finite_bt_temperature(row)
                               for row in collisions])
    stout_lines = STOUT.read_text().splitlines()
    stout_temperatures = [float(value) for value in stout_lines[1].split()[1:]]
    stout_first = {}
    for line in stout_lines[2:]:
        fields = line.split()
        if (len(fields) == 4 + len(stout_temperatures)
            and fields[:2] == ["CS", "ELECTRON"]):
            pair = (int(fields[2]), int(fields[3]))
            if pair[0] <= 40 and pair[1] <= 40:
                stout_first[pair] = float(fields[4])
    compared = []
    for row in collisions:
        pair = (row["ll"], row["ul"])
        if pair in stout_first:
            calculated = float(C.upsilon(row, np.array([5000.0]))[0])
            quoted = stout_first[pair]
            compared.append(abs(calculated / quoted - 1) if quoted else abs(calculated))
    levels, radiation, _ = C.bound_excitation_atom("o", 5)
    frequency = []
    weighted = []
    for lower, upper, wavelength, gf, rate in C.read_wgfa(CHIANTI / "o_5.wgfa"):
        if lower not in levels or upper not in levels or rate <= 0 or wavelength == 0:
            continue
        selected_gap = levels[upper][1] - levels[lower][1]
        frequency.append(abs(selected_gap * abs(wavelength) / 1e8 - 1))
        if gf > 0:
            inferred_gap = 1e8 * np.sqrt(rate * levels[upper][0] /
                                          (6.670251e15 * gf))
            weighted.append(abs(inferred_gap / selected_gap - 1))
    report = {
        "source_sha256": {str(path): digest(path) for path in
                          (CHIANTI / "o_5.scups", CHIANTI / "o_5.wgfa", STOUT, PAPER)},
        "bound_counts": [len(levels), len(radiation), len(collisions),
                         sum(len(row["xs"]) for row in collisions)],
        "inferred_first_finite_node_K": {
            "minimum": float(np.min(temperatures)),
            "median": float(np.median(temperatures)),
            "maximum": float(np.max(temperatures)),
            "rows": len(temperatures)},
        "stout_distributed_temperature_nodes_K": stout_temperatures,
        "matched_low_40_pairs": len(compared),
        "maximum_relative_low_40_stout_strength_difference_at_5000_K": max(compared),
        "raw_wavelength_vs_selected_energy": {
            "rows": len(frequency), "maximum_relative_gap_difference": max(frequency),
            "rows_exceeding_one_percent": int(sum(error > .01 for error in frequency))},
        "raw_A_over_gf_vs_selected_energy": {
            "rows": len(weighted), "maximum_relative_gap_difference": max(weighted),
            "rows_exceeding_one_percent": int(sum(error > .01 for error in weighted))},
        "paper_methods_temperature_sentence_K": [16000.0, 160000000.0],
        "paper_grid_matches_distributed_O_V_grid": False,
        "physical_model_accepted": False,
        "source_accuracy_verified": False,
    }
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
