"""Compare identified low-level CHIANTI, Stout, and installed oxygen rates."""

import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'tools/fitting'))
from chianti_cooling import read_elvlc, read_wgfa

CHIANTI = Path('/home/kiseon/RT_Codes/CHIANTI/dbase/o')
STOUT = Path('/home/kiseon/CLOUDY/c25.00/data/stout/o')
HERE = Path(__file__).resolve().parent


def installed(stage):
    path = ROOT / f'data/atomic/nlevel_o_{stage}.txt'
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('NRAD '))
    count = int(lines[start].split()[1])
    return path, {(int(x[0]), int(x[1])): float(x[2])
                  for line in lines[start + 1:start + count + 1]
                  if (x := line.split())}


def stout_rows(path, top):
    rates = {}
    types = {}
    for line in path.read_text().splitlines()[1:]:
        fields = line.split()
        if len(fields) < 4 or fields[0] != 'A':
            continue
        i, j = map(int, fields[1:3])
        if i > top or j > top:
            continue
        pair = (i, j)
        rates[pair] = rates.get(pair, 0.0) + float(fields[3])
        kind = fields[4] if len(fields) > 4 and not fields[4].startswith('#') else 'unspecified'
        types.setdefault(pair, []).append((kind, float(fields[3])))
    return rates, types


def main():
    result = {}
    files = []
    for stage, top in ((2, 5), (3, 6)):
        stem = f'o_{stage}'
        cp = CHIANTI / stem
        sp = STOUT / stem
        levels = read_elvlc(cp / f'{stem}.elvlc')
        wgfa = read_wgfa(cp / f'{stem}.wgfa')
        file, rates = installed(stage)
        stout, types = stout_rows(sp / f'{stem}.tp', top)
        files.extend((cp / f'{stem}.elvlc', cp / f'{stem}.wgfa',
                      sp / f'{stem}.nrg', sp / f'{stem}.tp', file))
        # CHIANTI and Stout keep the same identities but assign different
        # low-level indices to O II 2D and 2P fine-structure states.
        map_index = {1: 1, 2: 3, 3: 2, 4: 5, 5: 4} if stage == 2 else {i: i for i in range(1, top + 1)}
        source = {}
        for i, j, wavelength, gf, a in wgfa:
            if i > top or j > top or a <= 0:
                continue
            pair = (i, j) if levels[i][1] < levels[j][1] else (j, i)
            source[pair] = source.get(pair, 0.0) + a
        rows = []
        for pair in sorted(set(rates) | set(source)):
            stout_pair = tuple(sorted((map_index[pair[0]], map_index[pair[1]])))
            rows.append({'installed_pair': pair, 'installed_A_s1': rates.get(pair),
                         'CHIANTI_A_s1': source.get(pair),
                         'Stout_pair': stout_pair, 'Stout_A_s1': stout.get(stout_pair),
                         'Stout_components': types.get(stout_pair, []),
                         'installed_over_Stout': rates[pair] / stout[stout_pair]
                         if pair in rates and stout_pair in stout else None})
        reversed_rows = [(i, j, a) for i, j, _, _, a in wgfa
                         if i <= top and j <= top and a > 0 and levels[j][1] < levels[i][1]]
        result[stem] = {'level_count': top, 'CHIANTI_to_Stout_index': map_index,
                        'selected_energy_reversed_WGFA_rows': reversed_rows,
                        'radiative_rows': rows}
    result['sha256'] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in files}
    (HERE / 'ground_rate_comparison.json').write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
