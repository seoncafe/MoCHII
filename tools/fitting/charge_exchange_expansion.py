"""Audit charge-transfer recombination candidates at added stage boundaries.

The ten added rates have published mean heat contributions in KF99 Table 1.
They support selective staging under ground-state, low-density assumptions;
they do not identify product levels or constitute a complete plasma model.
Cloudy's last column is gas heat, not product excitation or KF96 rate data.
The first four rows per element correspond to reactant charges 1 through 4.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path


CLOUDY_DATA = Path.home() / "CLOUDY/c25.00/data/ctrecombdata.dat"

# Coefficients checked against Table 1 of the final KF96 journal article.
# Tuple order is (a [10^-9 cm^3 s^-1], b, c, d), with T4 = T / 10^4 K.
# KF96 Table 1 has N4+ (3.28, .52, -.52, -.19), whereas Cloudy's N4+
# row has (2.95, .55, -.39, -1.07); therefore N4+ is deliberately absent.
KF96_ADDED = {
    ("c", 4): (332.46, -0.11, -0.995, -0.00158),
    ("n", 3): (4.54, 0.57, -0.65, -0.89),
    ("o", 3): (3.98, 0.26, 0.56, -2.62),
    ("o", 4): (0.252, 0.63, 2.08, -4.16),
    ("ne", 3): (14.73, 0.0452, -0.84, -0.31),
    ("ne", 4): (6.47, 0.54, 3.59, -5.22),
    ("s", 4): (6.44, 0.13, 2.69, -5.69),
    ("ar", 4): (6.37, 0.85, 10.21, -6.22),
    ("mg", 3): (6.49, 0.53, 2.82, -7.63),
    ("mg", 4): (6.36, 0.55, 3.86, -5.19),
}

# KF96 final journal Table 1 intervals [K]. Cloudy is checked independently.
KF96_ADDED_INTERVALS = {
    ('c',4): (10.,1e5), ('n',3): (10.,1e5),
    ('o',3): (1e3,5e4), ('o',4): (1e3,3e4),
    ('ne',3): (5e3,5e4), ('ne',4): (1e3,3e4),
    ('s',4): (1e3,3e4), ('ar',4): (1e3,3e4), ('mg',3): (1e3,3e4),
    ('mg',4): (1e3,3e4),
}

# Mean kinetic heat [eV] for X^q+ + H(1s) -> X^(q-1)+ + H+.
# Kingdon & Ferland 1999, ApJ 516 L107, Table 1, DOI 10.1086/312008.
# Original publisher PDF SHA-256:
# 2c681fe32532baddee5d19a67e472cefc51f6a0135de01ab82eb9aa182e78def
# Ground-state reactants and radiative relaxation at low density are assumed.
# Individual product populations, diffuse cascade reabsorption, and high-density
# collisional relaxation are not determined by these constant mean energies.
# Rows may contain avoided-crossing estimates, rather than state-resolved data.
KF99_ADDED_HEAT = {
    ('c',4): 11.30, ('n',3): 6.40, ('o',3): 5.00, ('o',4): 8.47,
    ('ne',3): 5.82, ('ne',4): 8.60, ('s',4): 8.60,
    ('ar',4): 8.60, ('mg',3): 5.73, ('mg',4): 8.60,
}

# Further published total rates, separated from the nine Cloudy-matching rows.
# N4 uses the updated KF96 database distributed with Cloudy (Ferland et al.
# 2013, RMxAA 49, 137, Sec. 2.1: "an updated version of the Kingdon & Ferland
# (1996) database"; Cloudy c25.00 data/ctrecombdata.dat), 1e3-1e6 K, in place
# of the printed KF96 Table 1 fit (3.28, 0.52, -0.52, -0.19; 1e3-3e5 K, after
# Stancil et al. 1995), which it exceeds by up to 1.37x at 1e4 K (author's
# decision, 2026-10-11).
# Ca3/4 are Table 1 fits. Cl3/4 and Fe4 are Table 2 Landau-Zener estimates,
# evaluated at log10(T/K) = 3, 3.5, 4, 4.5 in KF96 Section 3.2. Adopt
# 1000--30000 K inside that sampled span; it is not a printed Table 2 interval.
# Fe4 has a radiative-transfer caveat: the nonradiative estimate may not dominate.
KF96_RECOMBINATION_TOTALS = {
    ("n",4): ((2.95,.55,-.39,-1.07),1e3,1e6),
    ("cl",3): ((1.88,.32,1.77,-5.70),1e3,3e4),
    ("cl",4): ((7.27,.29,1.04,-10.14),1e3,3e4),
    ("ca",3): ((.0317,2.12,12.06,-.40),1e3,3e4),
    ("ca",4): ((2.68,.69,-.68,-4.47),1e3,3e4),
    ("fe",4): ((14.60,.0357,-.92,-.37),1e3,3e4),
}
# KF99 Table 1 mean kinetic heat; the same ground-state/low-density assumptions
# documented above apply. These means do not determine reverse thermal rates.
KF99_RECOMBINATION_MEAN_HEAT = {
    ("n",4):11.00,("cl",3):8.00,("cl",4):9.00,
    ("ca",3):6.60,("ca",4):9.90,("fe",4):10.00,
}


# Existing MoCHII CX_HIGH form-6 rates. Intervals are read from the published
# KF96 Table 1, not from Cloudy. N2+ and O2+ use Barragan tables instead.
KF96_EXISTING = {
    ("c", 2): ((0.000167, 2.79, 304.72, -4.07), 5000., 50000.),
    ("c", 3): ((3.25, 0.21, 0.19, -3.29), 1000., 100000.),
    ("ne", 2): ((0.00001, 0., 0., 0.), 5000., 50000.),
    ("s", 2): ((0.00001, 0., 0., 0.), 1000., 30000.),
    ("s", 3): ((2.29, 0.0402, 1.59, -6.06), 1000., 30000.),
    ("ar", 2): ((0.00001, 0., 0., 0.), 1000., 30000.),
    ("ar", 3): ((4.57, 0.27, -0.18, -1.57), 1000., 30000.),
    ("mg", 2): ((0.0000858, 0.00249, 0.0293, -4.33), 1000., 30000.),
    ("fe", 2): ((1.26, 0.0772, -0.41, -7.31), 1000., 100000.),
    ("fe", 3): ((3.42, 0.51, -2.06, -8.99), 1000., 100000.),
    # Si2+ + H: the updated KF96 database distributed with Cloudy (Ferland et
    # al. 2013, Sec. 2.1), 500-1e5 K, in place of the printed KF96 fit
    # (1.23, 0.24, 3.17, 0.00418; 10-1e6 K), the same rule as N4+ (2026-10-11).
    ("si", 2): ((6.77, 0.0736, -0.43, -0.11), 500., 1e5),
    ("si", 3): ((0.49, -0.0874, -0.36, -0.79), 1000., 30000.),
    ("si", 4): ((7.58, 0.37, 1.06, -4.09), 1000., 50000.),
}

ATOMIC_NUMBER = {"c": 6, "n": 7, "o": 8, "ne": 10, "mg": 12,
                 "si": 14, "s": 16, "cl": 17, "ar": 18, "ca": 20,
                 "fe": 26}


@dataclass(frozen=True)
class RecombinationCandidate:
    element: str
    charge: int
    coefficients: tuple[float, float, float, float]
    min_temperature: float
    max_temperature: float
    cloudy_energy_defect_ev: float | None  # Positive heats gas for recombination.

    def rate(self, temperature: float) -> float:
        """KF96 fit in cm^3 s^-1, clamped to its stated fit interval."""
        if not math.isfinite(temperature) or temperature <= 0:
            raise ValueError("temperature must be positive and finite")
        t4 = min(max(temperature, self.min_temperature),
                 self.max_temperature) / 1e4
        a, b, c, d = self.coefficients
        return a * 1e-9 * t4**b * (1 + c * math.exp(d * t4))


def cloudy_recombination_rows(path: Path = CLOUDY_DATA):
    """Read Cloudy's four charge-transfer rows per element without inference.

    The C++ reader at atmdat_char_tran.cpp:929-941 reads rows by zero-based
    element and charge; HCTRecom(ion) uses charge = ion + 1. The file has 30
    element blocks, starting at H, after its magic number.
    """
    words = path.read_text().split()
    if int(words[0]) != 201903042:
        raise ValueError("unexpected Cloudy charge-transfer file version")
    values = [float(x) for x in words[1:]]
    if len(values) % 28:
        raise ValueError("incomplete Cloudy charge-transfer block")
    blocks = len(values) // 28
    if blocks < max(ATOMIC_NUMBER.values()):
        raise ValueError("Cloudy charge-transfer file ends before Ca")
    rows = {}
    for el, z in ATOMIC_NUMBER.items():
        for charge in range(1, 5):
            start = ((z - 1) * 4 + charge - 1) * 7
            row = values[start:start + 7]
            rows[el, charge] = RecombinationCandidate(
                el, charge, tuple(row[:4]), row[4], row[5], row[6])
    return rows


def published_added_rates(path: Path = CLOUDY_DATA):
    """Return only added rates matching the published KF96 coefficients.

    This checks rates and their published intervals; KF99_ADDED_HEAT supplies
    separately reviewed mean heat, without asserting product-state branching.
    """
    rows = cloudy_recombination_rows(path)
    verified = {}
    for key, expected in KF96_ADDED.items():
        row = rows[key]
        if any(not math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
               for a, b in zip(row.coefficients, expected)):
            raise ValueError(f"Cloudy/KF96 coefficient mismatch: {key}")
        if row.min_temperature <= 0 or row.max_temperature < row.min_temperature:
            raise ValueError(f"invalid charge-transfer temperature range: {key}")
        if (row.min_temperature,row.max_temperature) != KF96_ADDED_INTERVALS[key]:
            raise ValueError(f"Cloudy/KF96 temperature interval mismatch: {key}")
        if row.cloudy_energy_defect_ev != KF99_ADDED_HEAT[key]:
            raise ValueError(f"Cloudy/KF99 mean heat mismatch: {key}")
        verified[key] = row
    return verified


def kf96_added_rates():
    """Published added rate fits without a Cloudy runtime dependency.

    Constant KF99 heat is a separately documented mean-energy approximation.
    The reverse from excited products cannot be inferred from a total rate.
    """
    return {key: RecombinationCandidate(key[0],key[1],coefficients,
                *KF96_ADDED_INTERVALS[key],KF99_ADDED_HEAT[key])
            for key,coefficients in KF96_ADDED.items()}


def kf96_recombination_total_rates():
    """Further KF96 total recombination fits with their selected domains.

    For Table 2 estimates the selected 1000--30000 K domain lies inside the
    sampled temperature span, and carries the published Landau-Zener caveats.
    """
    return {key: RecombinationCandidate(key[0],key[1],coefficients,low,high,
                KF99_RECOMBINATION_MEAN_HEAT[key])
            for key,(coefficients,low,high) in KF96_RECOMBINATION_TOTALS.items()}


def published_existing_rates(path: Path = CLOUDY_DATA):
    """Return existing KF96 rates with their published temperature intervals.

    Every coefficient, Si2+ included (both now the updated KF96 database,
    Ferland et al. 2013), is checked against the Cloudy source as an
    independent transcription check. Cloudy's
    heat column is metadata only and must receive separate physical review.
    """
    rows = cloudy_recombination_rows(path)
    published = {}
    for key, (coefficients, low, high) in KF96_EXISTING.items():
        row = rows[key]
        if any(
            not math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)
            for a, b in zip(row.coefficients, coefficients)
        ):
            raise ValueError(f"Cloudy/KF96 coefficient mismatch: {key}")
        published[key] = RecombinationCandidate(
            key[0], key[1], coefficients, low, high,
            row.cloudy_energy_defect_ev)
    return published


def kf96_existing_rates():
    """Published rate fits without a Cloudy runtime dependency.

    The energy-defect field is None because KF96 Table 1 prints no kinetic
    energy and this function must not imply one. The generator uses only the
    coefficient and temperature fields; its separate CX_HEAT table supplies
    the existing thermal values, whose product-state assumptions require
    the separate provenance checks described in the audit.
    """
    return {
        key: RecombinationCandidate(key[0], key[1], coefficients, low, high,
                                    None)
        for key, (coefficients, low, high) in KF96_EXISTING.items()
    }
