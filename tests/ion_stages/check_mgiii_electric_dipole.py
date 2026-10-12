"""Check corrected bound Mg III energies, radiation, and source nodes."""

from hashlib import sha256
from pathlib import Path
import math
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "fitting"))
import chianti_cooling as C  # noqa: E402
import electric_dipole_probabilities as E1  # noqa: E402
import write_excitation_atom as W  # noqa: E402
from atomic_excitation_energies import selected_energies  # noqa: E402

SOURCE = Path(C.ion_dir("mg", 3)) / "mg_3"

assert E1.electric_dipole_probability_from_strength(2, 3, 3, 100, 200) == 16
assert E1.electric_dipole_probability_from_strength(2, 3, 6, 100, 100) == 1

raw_levels, raw_radiation, raw_collisions = C.read_excitation_atom("mg", 3)
assert (len(raw_levels), len(raw_radiation), len(raw_collisions)) == (209, 4753, 21736)
selected = selected_energies("mg", 3, raw_levels, raw_radiation, raw_collisions)
assert raw_levels[167] == (1.0, 628148.375, 634430.0)
assert selected[167] == (1.0, 634430.0)
assert all(selected[index][1] == (observed if observed > 0 or index == 1 else theory)
           for index, (g, observed, theory) in raw_levels.items() if index != 167)
try:
    changed = dict(raw_levels)
    changed[167] = (1.0, 628149.0, 634430.0)
    selected_energies("mg", 3, changed, raw_radiation, raw_collisions)
except ValueError:
    pass
else:
    raise AssertionError("Changed Mg III level 167 identity was accepted")
levels, radiation, collisions = C.bound_excitation_atom("mg", 3)
assert (len(levels), len(radiation), len(collisions)) == (167, 3364, 13861)
assert set(levels) == set(range(1, 168))
assert all(row["bt_tmin_K"] == 1800.0 for row in collisions)
assert sum(len(row["xs"]) for row in collisions) == 207915
raw_pairs = {(row["ll"], row["ul"]): row for row in raw_collisions}
for row in collisions:
    source = raw_pairs[row["ll"], row["ul"]]
    assert (row["xs"] == source["xs"]).all()
    assert (row["ys"] == source["ys"]).all()
assert E1.bound_mgiii_electric_dipole_probabilities(
    raw_levels, levels, [row for row in raw_radiation
                         if row[0] in levels and row[1] in levels],
    collisions, SOURCE) == radiation

raw_rates = {(low, high): rate for low, high, _, _, rate in raw_radiation}
physical_rates = {(low, high): rate for low, high, _, _, rate in radiation}
assert sum(max(low, high) <= 27 for low, high in physical_rates) == 117
reversed_count = 0
e1_count = unidentified_count = 0
for lower, upper, wavelength, gf, rate in raw_radiation:
    if lower not in levels or upper not in levels:
        continue
    physical_lower, physical_upper = sorted((lower, upper), key=lambda i: levels[i][1])
    reversed_count += physical_upper == lower
    physical_rate = physical_rates[physical_lower, physical_upper]
    if max(lower, upper) <= 27:
        assert physical_rate == rate
    elif gf > 0:
        theory_gap = raw_levels[upper][2] - raw_levels[lower][2]
        physical_gap = levels[physical_upper][1] - levels[physical_lower][1]
        expected = rate * (physical_gap / theory_gap) ** 3 * (
            raw_levels[upper][0] / raw_levels[physical_upper][0])
        assert math.isclose(physical_rate, expected, rel_tol=1e-14)
        e1_count += 1
    else:
        assert physical_rate == rate
        unidentified_count += 1
assert (reversed_count, e1_count, unidentified_count) == (14, 3135, 112)
assert all(levels[upper][1] > levels[lower][1]
           for lower, upper, _, _, rate in radiation if rate > 0)

source_hashes = E1.MGIII_SOURCE_SHA256
try:
    E1.MGIII_SOURCE_SHA256 = ("0" * 64,) * 3
    try:
        E1.bound_mgiii_electric_dipole_probabilities(
            raw_levels, levels, radiation, collisions, SOURCE)
    except ValueError as error:
        assert "source identity" in str(error)
    else:
        raise AssertionError("Changed Mg III source was accepted")
finally:
    E1.MGIII_SOURCE_SHA256 = source_hashes

with tempfile.TemporaryDirectory() as temporary:
    output = W.write_excitation_atom("mg", 3, output_dir=temporary)
    assert (output["nlev"], output["nrad"], output["nups"], output["nodes"]) == (
        167, 3364, 13861, 207915)
    generated = Path(output["path"]).read_bytes().split(b"NLEVEL_VERSION 2", 1)[1]
    if len(sys.argv) == 2:
        validated = Path(sys.argv[1]).read_bytes().split(b"NLEVEL_VERSION 2", 1)[1]
        assert generated == validated
        assert sha256(generated).hexdigest() == sha256(validated).hexdigest()
    else:
        assert len(sys.argv) == 1

print("Mg III source ID 167, 167 bound levels, 14 E1 reversals, first 27 A "
      "values, original collision nodes, and generated atom passed")
