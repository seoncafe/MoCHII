"""Check oxygen E1 line strength, source identity, and coronal cooling."""

import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "fitting"))
import chianti_cooling as C  # noqa: E402
from electric_dipole_probabilities import (  # noqa: E402
    electric_dipole_probability_from_strength as probability,
    observed_energy_electric_dipole_probabilities as transform,
)
from atomic_excitation_energies import selected_energies  # noqa: E402
from write_excitation_atom import _energy_pair  # noqa: E402

two_levels = {1: (1.0, 0.0), 2: (3.0, 100.0)}
assert _energy_pair(two_levels, 2, 1, collision=True) == (1, 2)
try:
    _energy_pair(two_levels, 2, 1)
except ValueError:
    pass
else:
    raise AssertionError("Unidentified radiative A was flipped to the physical direction")

assert probability(2.0, 3.0, 3.0, 100.0, 200.0) == 16.0
assert probability(2.0, 3.0, 6.0, 100.0, 100.0) == 1.0
assert probability(2.0, 3.0, 6.0, 100.0, 200.0) == 8.0
for invalid in ((0.0, 3.0, 3.0, 100.0, 100.0),
                (2.0, -1.0, 3.0, 100.0, 100.0),
                (2.0, 3.0, 3.0, 0.0, 100.0),
                (1e300, 3.0, 3.0, 1.0, 1e100)):
    try:
        probability(*invalid)
    except ValueError:
        pass
    else:
        raise AssertionError("Invalid E1 probability was accepted")

temperature = np.geomspace(1000.0, 100000.0, 401)
original_reader = C.read_excitation_atom
for stage, nlev, nrad, ncoll in ((2, 196, 4166, 19110), (3, 177, 3164, 15576)):
    stem = Path(C.ion_dir("o", stage)) / f"o_{stage}"
    levels = C.read_elvlc(str(stem) + ".elvlc")
    raw_radiation = C.read_wgfa(str(stem) + ".wgfa")
    collisions = C.read_scups(str(stem) + ".scups")
    assert (len(levels), len(raw_radiation), len(collisions)) == (nlev, nrad, ncoll)
    direct = transform(stage, levels, raw_radiation, collisions, stem)
    loaded_levels, loaded_radiation, loaded_collisions = original_reader("o", stage)
    assert loaded_levels == levels and len(loaded_collisions) == len(collisions)
    assert loaded_radiation == direct
    physical = selected_energies("o", stage, levels, direct, collisions)
    assert all(physical[upper][1] > physical[lower][1]
               for lower, upper, _, _, rate in direct if rate > 0)
    assert C.radiative_cascade_probabilities(physical, collisions, direct)[0] == [1]
    assert C.radiative_cascade_probabilities(physical, collisions, raw_radiation)[0] == [1]
    if stage == 2:
        fine = [row for row in direct if row[:2] == (3, 2)]
        assert len(fine) == 1 and fine[0][-1] == 1.32e-7
    else:
        rates = {(low, high): rate for low, high, _, _, rate in direct}
        source = {(low, high): rate for low, high, _, _, rate in raw_radiation}
        assert rates[2, 4] == 6.791e-3 and rates[3, 4] == 2.046e-2
        for pair in ((1, 2), (1, 4), (2, 3), (2, 5), (3, 5), (4, 5)):
            assert rates[pair] == source[pair]
    try:
        transform(stage, levels, raw_radiation[:-1], collisions, stem)
    except ValueError:
        pass
    else:
        raise AssertionError("Changed radiative source was accepted")
    try:
        transform(stage, levels, raw_radiation, collisions[:-1], stem)
    except ValueError:
        pass
    else:
        raise AssertionError("Changed collision source was accepted")

    def raw_reader(element, ion):
        if element.lower() == "o" and ion == stage:
            return levels, raw_radiation, collisions
        return original_reader(element, ion)

    try:
        C.read_excitation_atom = raw_reader
        before = C.cooling_low_density_limit("o", stage, temperature)
    finally:
        C.read_excitation_atom = original_reader
    after = C.cooling_low_density_limit("o", stage, temperature)
    error = float(np.max(abs(after / before - 1)))
    assert np.all(np.isfinite(after)) and np.all(after > 0) and error < 1e-14
    print(f"O {stage}: {nlev} levels, {nrad} radiation rows, {ncoll} collisions; "
          f"coronal relative difference {error:.3e}")
print("Oxygen E1 source, direction, low optical A, and coronal checks passed")
