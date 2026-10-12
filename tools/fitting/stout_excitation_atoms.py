"""Fe VII excitation atom from one consistent Stout ICFT source set.

The 189 bound levels, radiative rates, and effective collision strengths in
Cloudy c25.00 Stout derive from the Fe6+ calculation attributed in the files
to Witthoeft and Badnell (2008, A&A 481, 543).
The final Table 3 matches observed energies at levels 1--9 and calculated
energies at levels 10--189; printed Tables 4/5 rates match this source.
Source strengths cover 9800--9.8e7 K. At 1000--9800 K, the adopted approximation holds only
Upsilon at its lowest tabulated value; Maxwellian rates use actual T.
Low-temperature resonance accuracy is unverified. The source's 189 levels
are below the measured Fe VII ionization limit of 1007928 cm^-1.
The excitation model takes Fe VII from CHIANTI 11.0.2 (fe_7); this Stout atom
is kept for comparison (chianti_excitation_model_coverage.py).
"""

from hashlib import sha256
from itertools import combinations
from pathlib import Path
import math

import numpy as np

from chianti_cooling import HC_ERG_CM, K_B_ERG, RY_ERG


STEM = Path.home() / "CLOUDY/c25.00/data/stout/fe/fe_7/fe_7"
SOURCE_SHA256 = {
    "nrg": "42a060bca5c3dab121ca5c5cc4102235fc28976a1b1635962cafc3799f8c4415",
    "tp": "d3c426196550ccfe986b55c765d5407b27e1b1be6f5340678504a7dd0a417f74",
    "coll": "2fa169e62921cafda0490659a95d8a2d11183a3e14ad4635939879f321e0a72d",
}
IONIZATION_CM1 = 1007928.0
LEVEL_COUNT = 189
PAIR_COUNT = LEVEL_COUNT * (LEVEL_COUNT - 1) // 2
GRID_COUNT = 13


def _source_lines(stem, extension):
    path = Path(str(stem) + "." + extension)
    source = path.read_bytes()
    if sha256(source).hexdigest() != SOURCE_SHA256[extension]:
        raise ValueError(f"Fe VII Stout {extension} source changed: {path}")
    lines = source.decode("ascii").splitlines()
    if not lines or lines[0].strip() != "17 09 05":
        raise ValueError(f"Fe VII Stout {extension} format changed: {path}")
    return lines[1:]


def fevii_excitation_atom(stem=STEM):
    """Return checked levels, A rows, and collision rows for Fe VII.

    The return shape matches ``chianti_cooling.read_excitation_atom``. Type 2
    and ``xs`` define an invertible storage coordinate for native fits; the
    physical source is the temperature/Upsilon table in each collision row.
    ``chianti_cooling.upsilon`` applies the adopted constant cold continuation
    to that table while rate coefficients retain their actual temperature.
    """
    levels = {}
    for line in _source_lines(stem, "nrg"):
        if line.startswith("*"):
            break
        fields = line.split(maxsplit=3)
        if len(fields) < 4:
            raise ValueError("Invalid Fe VII Stout energy record")
        index, energy, weight = int(fields[0]), float(fields[1]), float(fields[2])
        if index in levels or not math.isfinite(energy) or not math.isfinite(weight):
            raise ValueError("Duplicate or nonfinite Fe VII Stout energy")
        if weight <= 0 or weight != int(weight):
            raise ValueError("Invalid Fe VII statistical weight")
        levels[index] = (weight, energy, energy)
    if set(levels) != set(range(1, LEVEL_COUNT + 1)) or levels[1][1] != 0:
        raise ValueError("Incomplete Fe VII Stout energy levels")
    if any(levels[i][1] >= levels[i + 1][1] for i in range(1, LEVEL_COUNT)):
        raise ValueError("Fe VII Stout levels are not in energy order")
    if levels[LEVEL_COUNT][1] >= IONIZATION_CM1:
        raise ValueError("Fe VII Stout atom includes a continuum level")

    expected_pairs = set(combinations(range(1, LEVEL_COUNT + 1), 2))
    radiative_pairs = set()
    radiation = []
    for line in _source_lines(stem, "tp"):
        if line.startswith("*"):
            break
        fields = line.split()
        if len(fields) != 4 or fields[0] != "A":
            raise ValueError("Invalid Fe VII Stout radiative record")
        lower, upper, rate = int(fields[1]), int(fields[2]), float(fields[3])
        pair = (lower, upper)
        if pair not in expected_pairs or pair in radiative_pairs or not math.isfinite(rate) or rate < 0:
            raise ValueError("Invalid Fe VII Stout radiative pair or rate")
        radiative_pairs.add(pair)
        if rate > 0:
            gap = levels[upper][1] - levels[lower][1]
            radiation.append((lower, upper, 1.0e8 / gap, 0.0, rate))
    if radiative_pairs != expected_pairs or len(radiation) != 17531:
        raise ValueError("Incomplete Fe VII Stout radiative rates")

    collision_pairs = set()
    transitions = []
    grid = None
    for line in _source_lines(stem, "coll"):
        if line.startswith("*"):
            break
        fields = line.split()
        if fields[0] == "TEMP":
            if grid is not None:
                raise ValueError("Duplicate Fe VII collision temperature grid")
            grid = np.asarray([float(value) for value in fields[1:]], dtype=float)
            if (grid.size != GRID_COUNT or np.any(~np.isfinite(grid))
                or np.any(np.diff(grid) <= 0) or grid[0] != 9800.0
                or grid[-1] != 9.8e7):
                raise ValueError("Invalid Fe VII collision temperature grid")
            continue
        if fields[:2] != ["CS", "ELECTRON"] or len(fields) != GRID_COUNT + 4 or grid is None:
            raise ValueError("Invalid Fe VII Stout collision record")
        lower, upper = int(fields[2]), int(fields[3])
        pair = (lower, upper)
        strength = np.asarray([float(value) for value in fields[4:]], dtype=float)
        if (pair not in expected_pairs or pair in collision_pairs
            or np.any(~np.isfinite(strength)) or np.any(strength <= 0)):
            raise ValueError("Invalid Fe VII Stout collision pair or strength")
        collision_pairs.add(pair)
        gap = levels[upper][1] - levels[lower][1]
        de = gap * HC_ERG_CM / RY_ERG
        et = K_B_ERG * grid / (de * RY_ERG)
        transitions.append({"ll": lower, "ul": upper, "de": de,
                            "gf": 0.0, "lim": -1.0, "ttype": 2,
                            "cups": 1.0, "xs": et / (1.0 + et),
                            "ys": strength, "temperature_grid": grid,
                            "upsilon_grid": strength})
    if collision_pairs != expected_pairs:
        raise ValueError("Incomplete Fe VII Stout collision strengths")
    return levels, radiation, transitions
