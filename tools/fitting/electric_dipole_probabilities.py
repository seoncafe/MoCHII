"""Physical E1 probabilities at selected observed/theoretical energies.

Original CHIANTI files remain unchanged. Source-version and transition-identity
guards prevent silent reuse after the atomic level assignments change.
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
MGIII_SOURCE_SHA256 = (
    "07e61cf56c782806986ca9a19c2197f2cfa85adaf568be4b0ab711850084851e",
    "bb00af4bf3369c4facb78b068e29df4d7dfb9c0a042ae2f21ca377bed44ab841",
    "9aba89f5d1f6ca158a88547420377d988df1363e98d0cbd0c948fb444a302402",
)


def electric_dipole_probability_from_strength(source_a, source_upper_g,
                                              physical_upper_g, source_frequency,
                                              physical_frequency):
    """Preserve E1 line strength under a physical photon-frequency change."""
    values = (source_a, source_upper_g, physical_upper_g,
              source_frequency, physical_frequency)
    if not all(math.isfinite(value) and value > 0 for value in values):
        raise ValueError("E1 probability requires finite positive source values")
    result = (source_a * (physical_frequency / source_frequency) ** 3
              * source_upper_g / physical_upper_g)
    if not math.isfinite(result) or result <= 0:
        raise ValueError("E1 probability is not finite and positive")
    return result


def observed_energy_electric_dipole_probabilities(stage, levels, radiation,
                                                  collisions, source_stem):
    """Return physical ordered O II/O III radiative tuples for selected energies.

    E1 line strengths are held fixed while the selected physical frequency
    sets A and gf. Other same-parity A values remain approximate positive
    source rates at the physical photon gap; their M1/E2 components are not
    all identified. The two O III ground optical rates are absolute totals
    from Storey and Zeippen (2000). Numerical balance was validated for
    1,000--100,000 K and n_e=1e-6--1e10 cm^-3. That test range is not an
    atomic-accuracy claim within or beyond those conditions. Returned gf for
    same-parity transitions remains source classification metadata, not a
    recalculated physical absorption oscillator strength. This includes the
    two O III optical lines whose absolute A totals replace the source A.
    """
    if stage not in SOURCE_SHA256:
        raise ValueError("Only O II and O III source identities are supported")
    paths = [Path(f"{source_stem}.{ext}") for ext in ("elvlc", "wgfa", "scups")]
    actual_hash = tuple(sha256(path.read_bytes()).hexdigest() for path in paths)
    if actual_hash != SOURCE_SHA256[stage]:
        raise ValueError(f"O {stage} CHIANTI source identity changed")
    expected = (196, 4166, 19110) if stage == 2 else (177, 3164, 15576)
    if (len(levels), len(radiation), len(collisions)) != expected:
        raise ValueError(f"O {stage} CHIANTI transition set changed")
    selected = selected_energies("o", stage, levels, radiation, collisions)
    configurations = {}
    terms = {}
    for line in paths[0].read_text().splitlines():
        if line.strip() == "-1":
            break
        index = int(line[:7])
        configurations[index] = line[7:37]
        terms[index] = (int(line[42:47]), line[47:52].strip())
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
            if reversed_order:
                if (configurations[2].strip() != "2s2 2p3"
                    or configurations[3].strip() != "2s2 2p3"
                    or terms[2] != (2, "D") or terms[3] != (2, "D")
                    or levels[2][0] != 4 or levels[3][0] != 6
                    or abs(selected[2][1] - 26830.6) > .001
                    or abs(selected[3][1] - 26810.8) > .001
                    or source_a != 1.32e-7):
                    raise ValueError("O II fine-structure M1 source identity changed")
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


def bound_mgiii_electric_dipole_probabilities(raw_levels, bound_levels,
                                               radiation, collisions, source_stem):
    """Orient bound Mg III photons and retain source E1 line strengths.

    Only 167 bound states are supported; the remaining raw states above the
    ionization limit require a separate branching model. The first 27 states
    retain Jonsson (2014) A values. For higher E1 transitions, selected photon
    frequency changes A as frequency cubed at fixed dipole strength. Higher
    gf=0 rates retain their source A as approximate effective rates because
    their multipole type is not identified. This approximation was tested for
    1,000--100,000 K and n_e=1e-6--1e10 cm^-3; it does not establish atomic
    accuracy there. BT collision data use a constant Upsilon below 1,800 K,
    and level-label and unidentified-multipole uncertainties remain. Returned
    gf for non-E1 rows is source metadata, not a physical absorption strength.
    """
    paths = [Path(f"{source_stem}.{ext}") for ext in ("elvlc", "wgfa", "scups")]
    actual_hashes = tuple(sha256(path.read_bytes()).hexdigest() for path in paths)
    if actual_hashes != MGIII_SOURCE_SHA256:
        raise ValueError("Mg III CHIANTI source identity changed")
    if (set(bound_levels) != set(range(1, 168)) or len(radiation) != 3364
            or len(collisions) != 13861):
        raise ValueError("Mg III bound transition set changed")
    configurations = {}
    for line in paths[0].read_text().splitlines():
        if line.strip() == "-1":
            break
        index = int(line[:7])
        configurations[index] = line[7:37]
    if not set(bound_levels) <= set(configurations):
        raise ValueError("Mg III bound configurations are incomplete")

    def parity(index):
        configuration = re.sub(r"\([^)]*\)", "", configurations[index])
        occupations = re.findall(r"(\d+)([spdfgh])(\d*)", configuration)
        if not occupations:
            raise ValueError(f"Mg III configuration cannot be classified: {index}")
        return sum((int(number) if number else 1) * "spdfgh".index(orbital)
                   for _, orbital, number in occupations) % 2

    physical = []
    first27_count = e1_count = unidentified_count = reversed_count = 0
    for source_lower, source_upper, wavelength, gf, source_a in radiation:
        if (source_lower not in bound_levels or source_upper not in bound_levels
                or not math.isfinite(source_a) or source_a <= 0):
            raise ValueError("Mg III bound radiative row changed")
        low, high = sorted((source_lower, source_upper),
                           key=lambda index: bound_levels[index][1])
        gap = bound_levels[high][1] - bound_levels[low][1]
        if not math.isfinite(gap) or gap <= 0:
            raise ValueError("Mg III nonpositive physical photon energy")
        reversed_order = high == source_lower
        if max(source_lower, source_upper) <= 27:
            if reversed_order:
                raise ValueError("Mg III Jonsson ground transition reversed")
            physical_a, physical_gf = source_a, gf
            first27_count += 1
        else:
            opposite = parity(source_lower) != parity(source_upper)
            j_lower = (raw_levels[source_lower][0] - 1) / 2
            j_upper = (raw_levels[source_upper][0] - 1) / 2
            e1 = (opposite and abs(j_upper - j_lower) <= 1
                  and j_upper + j_lower > 0 and gf > 0)
            if e1:
                theory_gap = raw_levels[source_upper][2] - raw_levels[source_lower][2]
                if theory_gap <= 0:
                    raise ValueError("Mg III nonpositive E1 source frequency")
                inferred = 1e8 * (source_a * raw_levels[source_upper][0]
                                  / (EINSTEIN * gf)) ** .5
                if not math.isfinite(inferred) or abs(inferred / theory_gap - 1) > .01:
                    raise ValueError("Mg III E1 source frequency mismatch")
                physical_a = electric_dipole_probability_from_strength(
                    source_a, raw_levels[source_upper][0], raw_levels[high][0],
                    theory_gap, gap)
                physical_gf = gf * gap / theory_gap
                e1_count += 1
            else:
                if gf > 0 or reversed_order:
                    raise ValueError("Mg III unidentified radiative transition")
                physical_a, physical_gf = source_a, gf
                unidentified_count += 1
        reversed_count += int(reversed_order)
        physical.append((low, high, 1e8 / gap, physical_gf, physical_a))
    if (first27_count, e1_count, unidentified_count, reversed_count) != (117, 3135, 112, 14):
        raise ValueError("Mg III radiative classifications changed")
    return physical


def configuration_parity(configuration):
    """Parity (0 even, 1 odd) of a CHIANTI configuration label.

    Parenthesized parent terms are removed before the orbital occupations are
    summed. A label without orbitals cannot be classified and is refused.
    """
    configuration = re.sub(r"\([^)]*\)", "", configuration)
    occupations = re.findall(r"(\d+)([spdfgh])(\d*)", configuration)
    if not occupations:
        raise ValueError(f"Configuration cannot be classified: {configuration!r}")
    return sum((int(number) if number else 1) * "spdfgh".index(orbital)
               for _, orbital, number in occupations) % 2


def reversed_electric_dipole_probabilities(raw_levels, selected, radiation,
                                           configurations):
    """Orient radiative rows whose selected energies reverse the source order.

    CHIANTI applies each .wgfa A value as a decay from its second index to its
    first. Where the selected energies (observed where measured, otherwise
    theoretical; CHIANTI's own rule) put the first index higher, that decay
    would create energy. For an E1 row (opposite parity, |Delta J| <= 1, not
    J = 0 - 0, gf > 0) the line strength S is a property of the pair of
    states, so the physical decay is from the higher selected state with
    A = A_source (nu_selected / nu_source)^3 g_source_upper / g_selected_upper,
    and gf scales as nu_selected / nu_source. nu_source is the frequency at
    which the source A and gf describe the same S,
    A = 6.6703e15 gf / (g_u lambda^2); it is the frequency of the source's
    own structure calculation, which may differ from both CHIANTI energy
    columns (CHIANTI Fe VII takes A from Tayal & Zatsarinny 2014 but its
    theoretical energies mostly from Li et al. 2018). Reversed rows that are
    not E1 have no identified multipole and are refused. Rows in the source
    order are returned unchanged. The transform assumes the state
    identifications are correct; it does not verify them.

    Returns (rows, count, unmatched) with rows as (lower, upper, wavelength,
    gf, A) in the selected physical order, count the number of transformed
    rows, and unmatched the transformed rows whose nu_source matches neither
    |wavelength| nor the theoretical gap within 2 percent, as
    (upper, lower, nu_source, |1/wavelength|, theoretical gap).
    """
    physical = []
    transformed = 0
    unmatched = []
    for source_lower, source_upper, wavelength, gf, source_a in radiation:
        if (source_a <= 0 or source_lower not in selected
                or source_upper not in selected
                or selected[source_upper][1] > selected[source_lower][1]):
            physical.append((source_lower, source_upper, wavelength, gf, source_a))
            continue
        gap = selected[source_lower][1] - selected[source_upper][1]
        if not math.isfinite(gap) or gap <= 0:
            raise ValueError(f"Zero selected photon energy {source_lower}-{source_upper}")
        j_lower = (raw_levels[source_lower][0] - 1) / 2
        j_upper = (raw_levels[source_upper][0] - 1) / 2
        e1 = (configuration_parity(configurations[source_lower])
              != configuration_parity(configurations[source_upper])
              and abs(j_lower - j_upper) <= 1 and j_lower + j_upper > 0 and gf > 0)
        if not e1:
            raise ValueError(f"Reversed non-E1 radiative row {source_lower}-{source_upper}; "
                             "identify its multipole before writing")
        source_frequency = 1e8 * (source_a * raw_levels[source_upper][0]
                                  / (EINSTEIN * gf)) ** .5
        listed = 1e8 / abs(wavelength) if wavelength != 0 else math.nan
        theory_gap = raw_levels[source_upper][2] - raw_levels[source_lower][2]
        if not (abs(listed / source_frequency - 1) <= .02
                or abs(theory_gap / source_frequency - 1) <= .02):
            unmatched.append((source_lower, source_upper, source_frequency,
                              listed, theory_gap))
        physical_a = electric_dipole_probability_from_strength(
            source_a, raw_levels[source_upper][0], raw_levels[source_lower][0],
            source_frequency, gap)
        physical.append((source_upper, source_lower, 1e8 / gap,
                         gf * gap / source_frequency, physical_a))
        transformed += 1
    return physical, transformed, unmatched
