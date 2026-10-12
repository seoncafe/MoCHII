"""Published Zhang et al. (2026) ground-level Fe RR/DR and author-data audit.

The coefficients and Eqs. (4)--(6) were verified in the final publisher PDF,
ApJS 285, 50, Tables 6/7, DOI 10.3847/1538-4365/ae80c0. Ground-level,
zero-density S1 DR and RR apply over q^2*(20--1e7) K. Compact DR fits differ
from author ADF09 totals at some nodes by more than the stated 2%; retain
this limitation rather than certify that accuracy. No density suppression or
metastable-parent recombination is supplied by these ground-level fits.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import zipfile
from pathlib import Path


ARCHIVE = (Path(__file__).resolve().parents[3] / "references"
           / "Zhang_2026_Fe3plus_to_Fe6plus_supplement.zip")
ARCHIVE_SHA256 = "ed44a5f76b824ec462c6807e3566d4bb82ebb0e4064e5856e398882a5f9d0033"
SOURCE_DOI = "10.3847/1538-4365/ae80c0"
DATA_DOI = "10.5281/zenodo.20264991"

# (c_i [cm^3 s^-1 K^(3/2)], E_i [K]), initial ground level, S1.
# Final publisher Table 6, independently verified including signs and units.
DR_FITS = {
    3: ((9.847e-7,3.927e-6,4.314e-5,1.218e-4,1.179e-3,2.724e-2),
        (205.8,2448.,14880.,48220.,330000.,682300.)),
    4: ((2.012e-6, 8.398e-6, 4.158e-5, 8.589e-5, 1.216e-3, 3.703e-2),
        (188.5, 1624., 9845., 45530., 262200., 743800.)),
    5: ((7.126e-6, 5.411e-5, 1.810e-4, 2.245e-4, 4.365e-3, 3.994e-2,
         5.979e-5),
        (501., 3001., 16270., 62810., 379200., 784700., 14360000.)),
    6: ((2.325e-5, 9.656e-5, 4.079e-4, 7.939e-3, 1.336e-1, 4.337e-3),
        (1017., 9864., 47770., 279400., 776800., 1638000.)),
}
RR_FITS = {
    3:(8.204e-11,.4252,143.6,2.877e7,.4732,3.001e5),
    4:(4.242e-10,.5864,26.78,3.415e7,.2466,3.415e5),
    5:(2.053e-9,.6747,4.608,3.438e7,.1254,3.840e5),
    6:(2.598e-9,.7035,6.016,3.394e7,.0792,3.937e5),
}
FINAL_PDF = ARCHIVE.parent / "Zhang_2026_ApJS_285_50.pdf"
FINAL_PDF_SHA256 = "ba91543363d26b11593c19ade6c731484a7ed8b39da7679493399c35cfa40d50"

SEQUENCES = {3:("CR","3D5",1),4: ("V", "3D4", 4), 5: ("TI", "3D3", 4),
             6: ("SC", "3D2", 3)}
NUMBER = re.compile(r"\d+\.\d+E[+-]\d+")


def dr_fit(charge: int, temperature_K: float) -> float:
    """Evaluate published zero-density, ground-level S1 DR [cm^3 s^-1].

    Published domain: q^2*(20--1e7) K. Elsewhere this is uncertified analytic
    extrapolation. The stated 2% accuracy is exceeded at some ADF09 nodes.
    """
    if charge not in DR_FITS or not math.isfinite(temperature_K) or temperature_K <= 0:
        raise ValueError("Charge must be 3, 4, 5, or 6 and temperature positive")
    coefficients, energies = DR_FITS[charge]
    return math.fsum(c * math.exp(-e / temperature_K)
                     for c, e in zip(coefficients, energies)) / temperature_K**1.5


def rr_fit(charge: int, temperature_K: float) -> float:
    """Published zero-density, ground-level RR [cm^3/s], Eqs. (5)--(6).

    The published domain is q^2*(20--1e7) K. Evaluation elsewhere is an
    analytic extrapolation without established accuracy, not endpoint holding.
    """
    if charge not in RR_FITS or not math.isfinite(temperature_K) or temperature_K<=0:
        raise ValueError("Unsupported charge or nonpositive temperature")
    a,b,t0,t1,c,t2=RR_FITS[charge]
    bt=b+c*math.exp(-t2/temperature_K)
    x=math.sqrt(temperature_K/t0)
    return a/(x*(1+x)**(1-bt)*(1+math.sqrt(temperature_K/t1))**(1+bt))


def verify_final_recombination_source():
    """Ensure selected coefficients refer to the checked final journal file."""
    if hashlib.sha256(FINAL_PDF.read_bytes()).hexdigest()!=FINAL_PDF_SHA256:
        raise ValueError("Zhang 2026 final journal PDF hash does not match verified source")


def read_ground_dr_grid(archive: Path, charge: int) -> tuple[tuple[float, float], ...]:
    """Read ADF09 final ALFT total for initial parent index 1 (ground level).

    ADAS ``xxdata_09`` documents ``diel_tot`` in cm^3 s^-1 and ``tea`` in K.
    ``ALFT(1)`` is the first initial parent, verified here against the parent
    level indexing and its zero excitation energy. The final ALFT section is
    used so resolved contributions are not inadvertently counted twice.
    """
    if charge not in SEQUENCES:
        raise ValueError("Only Fe3+ through Fe6+ are supported")
    member = f"supplement/adf09_Fe{charge}+_S1"
    with zipfile.ZipFile(archive) as files:
        source = files.read(member).decode("ascii")
    sequence, ground_config, parent_count = SEQUENCES[charge]
    header = source.splitlines()[0]
    if (f"SEQ='{sequence}" not in header or "NUCCHG=26" not in header
            or "ITYPE=3" not in header):
        raise ValueError(f"Unexpected ADF09 identity: {member}")
    parent = re.search(r"PARENT LEVEL INDEXING\s+BWNP=.*?NPRNT=\s*(\d+)", source)
    ground = re.search(r"^\s*1\s+3S2\s+3P6\s+(3D\d)\s+.*?\)\s+0\.0\s*$",
                       source, re.MULTILINE)
    if not parent or int(parent.group(1)) != parent_count or not ground or ground.group(1) != ground_config:
        raise ValueError(f"Ground parent identity cannot be verified: {member}")
    marker = "    T(K)     ALFT( 1)"
    if source.count(marker) != 1:
        raise ValueError(f"Missing or ambiguous ADF09 total: {member}")
    total = source.split(marker, 1)[1].split("C---", 1)[0]
    rows: list[tuple[float, float]] = []
    for line in total.splitlines():
        fields = line.split()
        if fields and NUMBER.fullmatch(fields[0]):
            if len(fields) != parent_count + 1 or not all(NUMBER.fullmatch(s) for s in fields):
                raise ValueError(f"Malformed ADF09 total row: {line}")
            values = tuple(map(float, fields))
            if values[0] <= 0 or values[1] <= 0:
                raise ValueError(f"Nonpositive ground DR grid row: {line}")
            rows.append((values[0], values[1]))
    if len(rows) != 19 or any(rows[i][0] <= rows[i-1][0] for i in range(1, len(rows))):
        raise ValueError(f"Unexpected ADF09 temperature grid: {member}")
    return tuple(rows)


def compare(archive: Path) -> dict:
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest != ARCHIVE_SHA256:
        raise ValueError(f"Author archive SHA-256 mismatch: {digest}")
    result = {"paper_doi": SOURCE_DOI, "data_doi": DATA_DOI,
              "archive_sha256": digest, "fit_table_status": "verified final publisher Tables 6/7",
              "charge_states": {}}
    for charge in (3, 4, 5, 6):
        rows = read_ground_dr_grid(archive, charge)
        deviations = [(temperature, dr_fit(charge, temperature) / dr - 1.)
                      for temperature, dr in rows]
        valid = [item for item in deviations if item[0] >= 20 * charge**2]
        worst = max(valid, key=lambda item: abs(item[1]))
        result["charge_states"][str(charge)] = {
            "transition": f"Fe {charge+1} to Fe {charge}",
            "grid_points": len(rows), "grid_temperature_K": [rows[0][0], rows[-1][0]],
            "published_fit_minimum_K": 20 * charge**2,
            "maximum_valid_grid_relative_difference": worst[1],
            "maximum_valid_grid_difference_temperature_K": worst[0],
            "grid": [{"temperature_K": t, "adf09_dr_cm3_s": dr,
                      "fit_dr_cm3_s": dr_fit(charge, t),
                      "relative_difference": deviation}
                     for (t, dr), (_, deviation) in zip(rows, deviations)],
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    args = parser.parse_args()
    print(json.dumps(compare(args.archive), indent=2))


if __name__ == "__main__":
    main()
