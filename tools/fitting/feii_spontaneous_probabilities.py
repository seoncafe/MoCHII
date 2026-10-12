"""Complete the spontaneous decays of the MOCASSIN v2 Fe II atom from CHIANTI.

The 'mocassin_v2' Fe II atom (142 levels, the Zhang & Pradhan level list of
MOCASSIN feii.dat) carries 1144 source radiative rows plus 40 forbidden A
values for levels 1-13. No radiative row leaves levels 14-35 (a4P, b4P, a4H,
b4F, a4G, b4D; 13474-31483 cm^-1): those even-parity levels of 3d7 and
3d6 4s have only forbidden (M1/E2) and intercombination decays, which the
source omits, so every excitation into them stays there until a collision
removes it.

CHIANTI 11.0.2 fe_2 is on the same 142-level list. This script maps each
level by index, statistical weight (2J+1 of .elvlc) and energy (observed
.elvlc value, otherwise theoretical) and refuses the whole update if any
level fails; it then appends every CHIANTI .wgfa row whose level pair is
absent from the atom. Rows already present are kept with their v2 values;
every existing line of the file is kept unchanged. A row whose photon
energy on the atom's own level energies is not positive is refused.

Run:  python3 tools/fitting/feii_spontaneous_probabilities.py
          [--atom data/atomic/excitation/mocassin_v2/nlevel_fe_2.txt]
          [--report REPORT.json]
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

from atomic_output import atomic_write_text
from chianti_cooling import DBASE, read_elvlc, read_wgfa

ROOT = Path(__file__).resolve().parents[2]
ATOM = ROOT / 'data/atomic/excitation/mocassin_v2/nlevel_fe_2.txt'
SOURCE = Path(DBASE) / 'fe/fe_2'
ENERGY_TOLERANCE_CM = 0.01   # cm^-1; the two lists agree to 0.005 cm^-1
MARKER = '# CHIANTI fe_2 spontaneous decays added:'
L_OF = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def configuration_parity(configuration):
    """Parity (0 even, 1 odd) of a CHIANTI configuration string, e.g. 3d6(5D).4p."""
    core = re.sub(r'\([^)]*\)', '', configuration)
    total = 0
    for _, letter, count in re.findall(r'(\d+)([spdfg])(\d*)', core):
        total += L_OF[letter] * (int(count) if count else 1)
    return total % 2


def chianti_levels():
    levels = read_elvlc(SOURCE / 'fe_2.elvlc')
    parity = {}
    with open(SOURCE / 'fe_2.elvlc') as stream:
        for line in stream:
            if line.strip() == '-1':
                break
            parity[int(line[:7])] = configuration_parity(line[7:37].strip())
    return levels, parity


def read_atom_text(path):
    lines = Path(path).read_text().splitlines()
    nlev_line = next(i for i, s in enumerate(lines) if s.startswith('NLEV '))
    nlev = int(lines[nlev_line].split()[1])
    levels = {}
    for s in lines[nlev_line + 1:nlev_line + 1 + nlev]:
        i, g, e = s.split()
        levels[int(i)] = (float(g), float(e))
    nrad_line = next(i for i, s in enumerate(lines) if s.startswith('NRAD '))
    nrad = int(lines[nrad_line].split()[1])
    rows = {}
    for s in lines[nrad_line + 1:nrad_line + 1 + nrad]:
        l, u, a = s.split()[:3]
        key = (int(l), int(u))
        if key in rows or (key[1], key[0]) in rows:
            raise ValueError(f'Duplicate radiative pair {key} in {path}')
        rows[key] = float(a)
    return lines, nlev_line, levels, nrad_line, nrad, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--atom', type=Path, default=ATOM)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    lines, nlev_line, levels, nrad_line, nrad, rows = read_atom_text(args.atom)
    if any(s.startswith(MARKER) for s in lines):
        raise ValueError('CHIANTI decays were already added to this atom')

    source_levels, parity = chianti_levels()
    if sorted(source_levels) != sorted(levels):
        raise ValueError('Level index sets differ: no identity mapping')
    worst = 0.0
    refused_levels = []
    for i, (g, e) in levels.items():
        sg, e_obs, e_th = source_levels[i]
        e_source = e_obs if (e_obs > 0 or i == 1) else e_th
        difference = abs(e - e_source)
        # The mapped level must also be the nearest CHIANTI level of equal g.
        nearest = min((k for k in source_levels if source_levels[k][0] == g),
                      key=lambda k: abs((source_levels[k][1] if source_levels[k][1] > 0 or k == 1
                                         else source_levels[k][2]) - e))
        if sg != g or difference > ENERGY_TOLERANCE_CM or nearest != i:
            refused_levels.append(dict(index=i, g=g, energy=e, chianti_g=sg,
                                       chianti_energy=e_source, nearest_same_g=nearest))
        worst = max(worst, difference)
    if refused_levels:
        raise ValueError(f'Level identities do not map: {refused_levels}')

    present = set(rows) | {(u, l) for l, u in rows}
    raw = {}
    with open(SOURCE / 'fe_2.wgfa') as stream:
        for line in stream:
            if line.strip() == '-1':
                break
            fields = line.split()
            raw[(int(fields[0]), int(fields[1]))] = fields[4]
    source_rows = read_wgfa(SOURCE / 'fe_2.wgfa')
    if len(source_rows) != len(raw):
        raise ValueError('Duplicate CHIANTI radiative pairs need a summing rule')
    added, refused, kept_differences = [], [], []
    for l, u, wavelength, gf, a in source_rows:
        if (l, u) in present:
            mine = rows.get((l, u), rows.get((u, l)))
            kept_differences.append(mine / a)
            continue
        photon = levels[u][1] - levels[l][1]
        if a <= 0 or photon <= 0 or l not in levels or u not in levels:
            refused.append(dict(lower=l, upper=u, A=a, photon_cm=photon))
            continue
        added.append((l, u, raw[(l, u)], a, parity[l] != parity[u]))
    if refused:
        raise ValueError(f'Refused CHIANTI rows: {refused}')
    missing_in_chianti = sorted(k for k in rows if k not in raw and (k[1], k[0]) not in raw)

    upper_sinks = {u for _, u in rows} | {u for _, u, *_ in added}
    undecayed = [i for i in range(2, len(levels) + 1) if i not in upper_sinks]
    trapped_before = [i for i in range(14, 36) if i not in {u for _, u in rows}]
    n_e1 = sum(1 for *_, odd in added if odd)
    header = [
        f'{MARKER} {len(added)} rows ({len(added) - n_e1} parity-conserving forbidden, '
        f'{n_e1} parity-changing), levels mapped by index, g and energy (max |dE| {worst:.4f} cm^-1)',
        '# Added rows: every CHIANTI 11.0.2 fe_2.wgfa pair absent from the v2 rows; A copied unchanged '
        '(the source A of the upper level); existing v2 rows and values unchanged',
        f'# source SHA-256 {SOURCE / "fe_2.wgfa"}: {sha256(SOURCE / "fe_2.wgfa")}',
        f'# source SHA-256 {SOURCE / "fe_2.elvlc"}: {sha256(SOURCE / "fe_2.elvlc")}',
        f'# generated by tools/fitting/feii_spontaneous_probabilities.py, SHA-256 {sha256(__file__)}',
    ]
    first_keyword = next(i for i, s in enumerate(lines) if not s.startswith('#'))
    new_rows = [f'{l} {u} {text}' for l, u, text, _, _ in added]
    output = (lines[:first_keyword] + header + lines[first_keyword:nrad_line]
              + [f'NRAD {nrad + len(added)}'] + lines[nrad_line + 1:nrad_line + 1 + nrad]
              + new_rows + lines[nrad_line + 1 + nrad:])
    atomic_write_text(args.atom, '\n'.join(output) + '\n')

    report = dict(
        atom=str(args.atom), atom_sha256=sha256(args.atom),
        chianti_wgfa_sha256=sha256(SOURCE / 'fe_2.wgfa'),
        chianti_elvlc_sha256=sha256(SOURCE / 'fe_2.elvlc'),
        nrad_before=nrad, nrad_after=nrad + len(added), added_rows=len(added),
        added_parity_changing=n_e1,
        added_with_upper_level_14_35=sum(1 for _, u, *_ in added if 14 <= u <= 35),
        levels_14_35_without_decay_before=trapped_before,
        levels_without_decay_after=undecayed,
        maximum_level_energy_difference_cm=worst,
        shared_pairs=len(kept_differences),
        shared_pairs_ratio_v2_to_chianti_min=min(kept_differences),
        shared_pairs_ratio_v2_to_chianti_max=max(kept_differences),
        v2_pairs_absent_from_chianti=len(missing_in_chianti))
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
