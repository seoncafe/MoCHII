"""Check Si I source identity and exact validated native atom records."""

from collections import Counter, defaultdict
import hashlib
from pathlib import Path
import tempfile

import numpy as np

import fit_cooling as F
from chianti_cooling import bound_excitation_atom, read_excitation_atom
from silicon_source_records import (CHI_ELVLC_SHA256, CHI_SCUPS_SHA256,
    CHI_WGFA_SHA256, NIST_SHA256, NIST_SOURCE, chianti_levels, nist_records)
from write_excitation_atom import write_excitation_atom


ROOT = Path(__file__).resolve().parents[2]
CANDIDATE = (ROOT / "md/atomic_stage_expansion/validation_2026-10-10/"
             "si_i_electric_dipole_source_audit/candidate_v2/nlevel_si_1.txt")
CHI_STEM = Path("/home/kiseon/RT_Codes/CHIANTI/dbase/si/si_1/si_1")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def numerical_records(path):
    return b"\n".join(line for line in Path(path).read_bytes().splitlines()
                      if line.strip() and not line.lstrip().startswith(b"#")) + b"\n"


def main():
    inputs = {Path(f"{CHI_STEM}.elvlc"): CHI_ELVLC_SHA256,
              Path(f"{CHI_STEM}.wgfa"): CHI_WGFA_SHA256,
              Path(f"{CHI_STEM}.scups"): CHI_SCUPS_SHA256,
              NIST_SOURCE: NIST_SHA256}
    assert all(digest(path) == expected for path, expected in inputs.items())
    _, source_radiation, source_collisions = read_excitation_atom("si", 1)
    levels, radiation, collisions = bound_excitation_atom("si", 1)
    assert (len(levels), len(radiation), len(collisions),
            sum(len(row["xs"]) for row in collisions)) == (35, 250, 595, 8925)
    assert set(levels) == set(range(1, 36))
    assert all(row["bt_tmin_K"] == 1003.02 for row in collisions)
    assert all(rate > 0 and levels[lower][1] < levels[upper][1]
               for lower, upper, _, _, rate in radiation)
    assert len({(lower, upper) for lower, upper, *_ in radiation}) == 250
    source_pairs = {(row["ll"], row["ul"]): row for row in source_collisions
                    if row["ll"] in levels and row["ul"] in levels}
    assert len(source_pairs) == len(collisions)
    for row in collisions:
        original = source_pairs[row["ll"], row["ul"]]
        assert all(row[key] == original[key] for key in ("ttype", "cups", "de"))
        assert np.array_equal(row["xs"], original["xs"])
        assert np.array_equal(row["ys"], original["ys"])

    original = {(min(l, u), max(l, u)): (l, u, gf, rate)
                for l, u, _, gf, rate in source_radiation if l in levels and u in levels}
    assert len(original) == 229
    corrected = {(l, u): rate for l, u, _, _, rate in radiation}
    assert original[14, 15][:2] == (15, 14)
    assert (14, 15) in corrected
    assert abs(corrected[14, 15] - 5.631292575070227) < 1e-13
    assert len(set(corrected) - set(original)) == 21

    _, identities = chianti_levels(Path(f"{CHI_STEM}.elvlc"))
    records, unmapped = nist_records(NIST_SOURCE, identities)
    assert not unmapped
    assert Counter(row["multipole"] for row in records) == {
        "E1": 196, "E2": 13, "M1": 10, "M2": 16}
    grouped = defaultdict(list)
    for row in records:
        pair = tuple(sorted((row["first_chianti_id"], row["second_chianti_id"])))
        grouped[pair].append(row)
    doubled = [pair for pair, (_, _, gf, _) in original.items()
               for row in grouped[pair] if row["multipole"] == "E1"
               and 1.9 < gf / row["length_GF"] < 2.1]
    assert len(doubled) == 15
    assert sum(any(row["multipole"] == "E1" for row in group)
               and any(row["multipole"] == "M2" for row in group)
               and pair in original for pair, group in grouped.items()) == 11

    with tempfile.TemporaryDirectory() as directory:
        metadata = write_excitation_atom("si", 1, output_dir=directory)
        generated = Path(metadata["path"])
        assert metadata["passed"] is True
        assert tuple(metadata[key] for key in ("nlev", "nrad", "nups", "nodes")) == (
            35, 250, 595, 8925)
        assert metadata["bt_tmin_K"] == 1003.02
        assert numerical_records(generated) == numerical_records(CANDIDATE)
        assert b"BT_TMIN_K 1003.02\n" in numerical_records(generated)
        coefficient_path = Path(F.write_ion_file(
            "SiI", "si", 1, np.array([1e-20]), np.array([1000.0]),
            0.001, 0.001, directory))
        header = coefficient_path.read_text()
        assert "NIST 2002 MCHF" in header
        assert "1003.02 K" in header
        assert "not a published physical interval" in header
        assert "on 201 source temperatures" in header
        assert F.IONS["SiI"] == ("si", 1)
    assert all(digest(path) == expected for path, expected in inputs.items())
    print("Si I source identity, radiation, original collision nodes, and numeric atom pass")


if __name__ == "__main__":
    main()
