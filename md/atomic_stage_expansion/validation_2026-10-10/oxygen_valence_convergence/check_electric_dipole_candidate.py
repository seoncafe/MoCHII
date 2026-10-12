"""Focused physical identities for the isolated oxygen E1 transformation."""

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools" / "fitting"))
from chianti_cooling import read_excitation_atom  # noqa: E402
from electric_dipole_probabilities_candidate import (  # noqa: E402
    electric_dipole_probability_from_strength as probability,
    observed_energy_electric_dipole_probabilities as transform,
)

assert probability(2.0, 3.0, 3.0, 100.0, 200.0) == 16.0
assert probability(2.0, 3.0, 6.0, 100.0, 100.0) == 1.0
assert probability(2.0, 3.0, 6.0, 100.0, 200.0) == 8.0
for invalid in ((0.0, 3.0, 3.0, 100.0, 100.0),
                (2.0, -1.0, 3.0, 100.0, 100.0),
                (2.0, 3.0, 3.0, 0.0, 100.0)):
    try:
        probability(*invalid)
    except ValueError:
        pass
    else:
        raise AssertionError("Invalid physical factor was accepted")

for stage, expected_inverted in ((2, 31), (3, 48)):
    raw, radiation, collisions = read_excitation_atom("o", stage)
    stem = Path.home() / f"RT_Codes/CHIANTI/dbase/o/o_{stage}/o_{stage}"
    physical = transform(stage, raw, radiation, collisions, stem)
    assert len(physical) == len(radiation)
    assert all(rate > 0 and math.isfinite(rate) for _, _, _, _, rate in physical)
    if stage == 2:
        ground_m1 = [row for row in physical if row[0:2] == (3, 2)]
        assert len(ground_m1) == 1 and ground_m1[0][-1] == 1.32e-7
    else:
        rates = {(low, high): rate for low, high, _, _, rate in physical}
        assert rates[2, 4] == 6.791e-3 and rates[3, 4] == 2.046e-2
        old = {(low, high): rate for low, high, _, _, rate in radiation}
        for pair in ((1, 2), (1, 4), (2, 3), (2, 5), (3, 5), (4, 5)):
            assert rates[pair] == old[pair]
    try:
        transform(stage, raw, radiation[:-1], collisions, stem)
    except ValueError:
        pass
    else:
        raise AssertionError("Transition-count identity change was accepted")
    print(f"O {stage}: {len(physical)} physical radiative rows; {expected_inverted} source inversions")
print("Electric-dipole frequency, weight, direction, and source guards passed")
