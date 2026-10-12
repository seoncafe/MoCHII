"""Proposed oxygen radiative transformation for post-validation integration.

This is an isolated draft, not imported by production source readers. Original
CHIANTI files remain unchanged. A source-version guard prevents silent reuse
after any level identity or transition probability changes.
"""

from hashlib import sha256
import math
from pathlib import Path
import re

from atomic_excitation_energies import selected_energies

EINSTEIN = 6.6703e15
SOURCE_SHA256 = {
    2: ("ae169f011f2564828640380b85367ad330c2cff77b06c65249b659a02306685d",
        "570ea65bfc4c646fde8f1a7cb008eb9f206088c2eb765fb43c99b3d1d6e23620",
        "0614d5b3ebe5e655e82b89a38fc5cfdfb6d294faddc02b592f6bb37be6a78c01"),
    3: ("2bd4da43c855325c8d4c5c6640f539945106219e51231095596797801ec7e297",
        "19b9363a630a092eb24963f72ee9b201b8354ed3b0535ec0c4da0c3be922d2e3",
        "c45303a6d6a4dd931aece7c1023ac80d19692c35d78de534b4c145ddfe7e9333"),
}
SZ_OPTICAL = {(2, 4): 6.791e-3, (3, 4): 2.046e-2}


def electric_dipole_probability_from_strength(source_a, source_upper_g,
                                              physical_upper_g, source_frequency,
                                              physical_frequency):
    """Preserve E1 line strength under a physical photon-frequency change."""
    values = (source_a, source_upper_g, physical_upper_g,
              source_frequency, physical_frequency)
    if not all(math.isfinite(value) and value > 0 for value in values):
        raise ValueError("E1 probability requires finite positive source values")
    return (source_a * (physical_frequency / source_frequency) ** 3
            * source_upper_g / physical_upper_g)


def observed_energy_electric_dipole_probabilities(stage, levels, radiation,
                                                  collisions, source_stem):
    """Return physical ordered O II/O III radiative tuples for selected energies.

    E1 line strengths are held fixed while the selected physical frequency
    sets A and gf. Other A values remain source estimates, except the two
    identified O III optical ground rates from Storey and Zeippen (2000).
    Validated numerical use is 1,000--100,000 K; density tests are required
    through 1e10 cm^-3. This is not a precision claim for atomic A values.
    """
    if stage not in SOURCE_SHA256:
        raise ValueError("Only O II and O III source identities are supported")
    paths = [Path(f"{source_stem}.{ext}") for ext in ("elvlc", "wgfa", "scups")]
    actual_hash = tuple(sha256(path.read_bytes()).hexdigest() for path in paths)
    if actual_hash != SOURCE_SHA256[stage]:
        raise ValueError(f"O {stage} CHIANTI source identity changed")
    expected = (196, 4166) if stage == 2 else (177, 3164)
    if (len(levels), len(radiation)) != expected:
        raise ValueError(f"O {stage} CHIANTI transition set changed")
    selected = selected_energies("o", stage, levels, radiation, collisions)
    configurations = {}
    for line in paths[0].read_text().splitlines():
        if line.strip() == "-1":
            break
        configurations[int(line[:7])] = line[7:37]
    if set(configurations) != set(levels):
        raise ValueError(f"O {stage} level configurations are incomplete")

    def parity(index):
        occupations = re.findall(r"(\d+)([spdfgh])(\d*)", configurations[index])
        if not occupations:
            raise ValueError(f"O {stage} configuration cannot be classified: {index}")
        return sum((int(number) if number else 1) * "spdfgh".index(orbital)
                   for _, orbital, number in occupations) % 2

    physical = []
    inverted_e1 = 0
    inverted_non_e1 = 0
    for source_lower, source_upper, wavelength, gf, source_a in radiation:
        if source_a <= 0:
            physical.append((source_lower, source_upper, wavelength, gf, source_a))
            continue
        low, high = sorted((source_lower, source_upper), key=lambda i: selected[i][1])
        gap = selected[high][1] - selected[low][1]
        if gap <= 0:
            raise ValueError(f"O {stage} zero or negative physical photon energy")
        reversed_order = high == source_lower
        opposite = parity(source_lower) != parity(source_upper)
        if opposite:
            if gf <= 0:
                raise ValueError(f"O {stage} unidentified opposite-parity transition")
            j_lower = (levels[source_lower][0] - 1) / 2
            j_upper = (levels[source_upper][0] - 1) / 2
            if abs(j_lower - j_upper) > 1 or j_lower == j_upper == 0:
                raise ValueError(f"O {stage} E1 angular momentum mismatch")
            inferred = 1e8 * (source_a * levels[source_upper][0] / (EINSTEIN * gf)) ** .5
            frequencies = [1e8 / abs(wavelength)]
            theoretical_gap = levels[source_upper][2] - levels[source_lower][2]
            if theoretical_gap > 0:
                frequencies.append(theoretical_gap)
            source_frequency = min(frequencies, key=lambda value: abs(value / inferred - 1))
            if abs(source_frequency / inferred - 1) > .01:
                raise ValueError(f"O {stage} source E1 frequency mismatch")
            ratio = gap / source_frequency
            physical_a = electric_dipole_probability_from_strength(
                source_a, levels[source_upper][0], levels[high][0],
                source_frequency, gap)
            physical_gf = gf * ratio
            inverted_e1 += int(reversed_order)
        else:
            if reversed_order and (stage, source_lower, source_upper) != (2, 2, 3):
                raise ValueError(f"O {stage} unidentified same-parity inversion")
            physical_a, physical_gf = source_a, gf
            inverted_non_e1 += int(reversed_order)
        model_wavelength = 1e8 / gap
        if levels[low][1] <= 0 and low != 1 or levels[high][1] <= 0:
            model_wavelength = -model_wavelength
        if stage == 3 and (low, high) in SZ_OPTICAL:
            physical_a = SZ_OPTICAL[low, high]
        physical.append((low, high, model_wavelength, physical_gf, physical_a))
    if (inverted_e1, inverted_non_e1) != ((30, 1) if stage == 2 else (48, 0)):
        raise ValueError(f"O {stage} source inversion identities changed")
    return physical
