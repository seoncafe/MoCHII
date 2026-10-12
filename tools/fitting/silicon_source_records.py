"""Identity-map the original NIST MCHF Si I transition records.

The NIST source file is preserved byte for byte under data/atomic/source/si_i.
This parser interprets its level labels and multipole line strengths; it does
not modify the original scientific data.
"""

import re
from pathlib import Path

NIST_SHA256 = "438150a064b12c202433f590122da7eaf819287e4656d0a6c5ba88586744d904"
CHI_ELVLC_SHA256 = "1538e98411f8c2f03414aebd36f901f7b88bc79109b1f0aea7254af070bd4f57"
CHI_WGFA_SHA256 = "be54f54d2e91695b766c0cb545e744d0d540165add1231775b7361f82ab21c42"
CHI_SCUPS_SHA256 = "102b70aa6020b9ff778ab4a504ca46dcc3fdbfcad6c2b5f1dc35c8dc8dc6eb5b"
NIST_SOURCE = Path(__file__).resolve().parents[2] / "data/atomic/source/si_i/Si_14.16.cff-mchf.lsj.db"
LEVEL_RE = re.compile(r"^\s*(\d+)\s+(-\d+\.\d+)\s+(\S+)\s*$")
ENERGY_RE = re.compile(r"^\s*([\d.]+) CM-1\s+([\d.]+) ANGS\(VAC\)")
RATE_RE = re.compile(r"^\s*(E1|E2|M1|M2)\s+S\s*=\s*([\d.eE+-]+)\s+GF\s*=\s*([\d.eE+-]+)\s+AKI\s*=\s*([\d.eE+-]+)")
VELOCITY_RE = re.compile(r"^\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s*$")
ORBITAL_RE = re.compile(r"(\d[spdfg])(?:(\d+)|\((\d+)\))?")
ORBITAL_ANGULAR_MOMENTUM = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}


def level_identity(configuration, term, twice_j):
    orbitals = ORBITAL_RE.findall(configuration)
    if not orbitals:
        raise ValueError(f"Unparsed orbital configuration: {configuration}")
    return (tuple((orbital, int(bare_count or parenthesized_count or 1))
                  for orbital, bare_count, parenthesized_count in orbitals),
            term.rstrip("*"), int(twice_j))

def nist_identity(label, twice_j):
    configuration, term = label.rsplit("_", 1)
    return level_identity(configuration, term, twice_j)

def chianti_levels(path):
    levels = {}
    identities = {}
    for line in path.read_text().splitlines():
        if line.lstrip().startswith("-1"):
            break
        fields = line.split()
        if not fields:
            continue
        index = int(fields[0])
        configuration = " ".join(fields[1:-5])
        term = fields[-5] + fields[-4]
        twice_j = round(2 * float(fields[-3]))
        identity = level_identity(configuration, term, twice_j)
        if identity in identities:
            raise ValueError(f"Ambiguous CHIANTI identity: {identity}")
        identities[identity] = index
        levels[index] = {"identity": identity, "observed_cm_inverse": float(fields[-2]),
                         "theoretical_cm_inverse": float(fields[-1])}
    return levels, identities

def nist_records(path, identities):
    lines = path.read_text().splitlines()
    records = []
    unmapped = []
    for line_number, line in enumerate(lines):
        rate = RATE_RE.match(line)
        if rate is None:
            continue
        if line_number < 3:
            raise ValueError("Incomplete NIST transition record")
        first = LEVEL_RE.match(lines[line_number - 3])
        second = LEVEL_RE.match(lines[line_number - 2])
        energy = ENERGY_RE.match(lines[line_number - 1])
        if first is None or second is None or energy is None:
            raise ValueError(f"Unexpected NIST record at line {line_number + 1}")
        first_identity = nist_identity(first[3], int(first[1]))
        second_identity = nist_identity(second[3], int(second[1]))
        first_id, second_id = identities.get(first_identity), identities.get(second_identity)
        record = {
            "nist_line": line_number + 1,
            "first_label": first[3], "second_label": second[3],
            "first_chianti_id": first_id, "second_chianti_id": second_id,
            "first_g": int(first[1]) + 1, "second_g": int(second[1]) + 1,
            "first_total_energy_hartree": float(first[2]),
            "second_total_energy_hartree": float(second[2]),
            "multipole": rate[1], "source_gap_cm_inverse": float(energy[1]),
            "source_wavelength_angstrom": float(energy[2]),
            "length_S_atomic_units": float(rate[2]),
            "length_GF": float(rate[3]), "length_A_per_second": float(rate[4]),
        }
        if line_number + 1 < len(lines):
            velocity = VELOCITY_RE.match(lines[line_number + 1])
            if velocity is not None:
                record.update({"velocity_S_atomic_units": float(velocity[1]),
                               "velocity_GF": float(velocity[2]),
                               "velocity_A_per_second": float(velocity[3])})
        records.append(record)
        if first_id is None or second_id is None:
            unmapped.append(record)
    return records, unmapped

def parity(identity):
    return sum(ORBITAL_ANGULAR_MOMENTUM[orbital[-1]] * count
               for orbital, count in identity[0]) % 2
