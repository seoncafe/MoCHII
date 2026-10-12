"""Checked N IV energies for the CHIANTI 136-state rate model.

Cloudy c25.00 Stout uses observed energies for states 1--29 and BSR energies
for the rest. This splice is an approximation to the original scattering
calculation: low-temperature resonance positions remain uncertain, and the
published comparison does not establish convergence above state 87. The
1,000--100,000 K fit interval is a numerical interval, not an accuracy claim.
The original CHIANTI transition indices and Einstein A values remain fixed.
"""

from hashlib import sha256
from pathlib import Path


STOUT_NRG = Path.home() / "CLOUDY/c25.00/data/stout/n/n_4/n_4.nrg"
STOUT_SHA256 = "e4ab82978db50833727e9b3cd4d4460834d63840a10e8ecfd829e2163d9f94cf"
N_STOUT = 238
N_BOUND = 136
N_RADIATIVE = 3272
N_COLLISION = 9180


def niv_energies(chianti_levels, radiation, collisions, source=STOUT_NRG):
    """Return all 238 (statistical weight, cm^-1) rows after source checks.

    Require the complete 136-state transition set before selecting a subset.
    This prevents an edited or mismatched source from producing output.
    """
    source = Path(source)
    try:
        raw = source.read_bytes()
    except OSError as error:
        raise ValueError(f"N IV Stout energy source unavailable: {source}") from error
    if sha256(raw).hexdigest() != STOUT_SHA256:
        raise ValueError(f"N IV Stout energy source hash changed: {source}")
    energies = {}
    for line in raw.decode("ascii").splitlines():
        if "\t" not in line:
            continue
        fields = line.split(maxsplit=3)
        if len(fields) < 3 or not fields[0].isdigit():
            continue
        index = int(fields[0])
        if index in energies:
            raise ValueError(f"N IV Stout duplicate level {index}")
        energies[index] = (float(fields[2]), float(fields[1]))
    expected = set(range(1, N_STOUT + 1))
    if set(energies) != expected or set(chianti_levels) != expected:
        raise ValueError("N IV Stout and CHIANTI must have matching levels 1--238")
    for index in expected:
        if energies[index][0] != chianti_levels[index][0]:
            raise ValueError(f"N IV statistical weight differs at level {index}")
    if len(radiation) != N_RADIATIVE or len(collisions) != N_COLLISION:
        raise ValueError("N IV radiative or collision source count changed")
    for lower, upper, _, _, rate in radiation:
        if not (1 <= lower < upper <= N_BOUND and rate > 0 and
                energies[upper][1] > energies[lower][1]):
            raise ValueError(f"N IV radiative pair {lower}->{upper} is inconsistent")
    for transition in collisions:
        lower, upper = transition["ll"], transition["ul"]
        if not (1 <= lower < upper <= N_BOUND and
                energies[upper][1] > energies[lower][1]):
            raise ValueError(f"N IV collision pair {lower}->{upper} is inconsistent")
    return energies


