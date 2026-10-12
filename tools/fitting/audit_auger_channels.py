#!/usr/bin/env python3
"""Audit shell-resolved Auger destinations inside the 100 eV photon band.

Read the original phfit2 tables through the repository parser and the Cloudy
transcription of Kaastra & Mewe (1993) Table 2. This is a reachability audit,
not a replacement for shell-resolved transport rates or Auger electron spectra.
"""
import argparse
import json
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools/fitting'))
import make_element_data as source

TARGET = {'c': 5, 'n': 6, 'o': 5, 'ne': 5, 's': 7, 'ar': 7,
          'mg': 4, 'fe': 7, 'si': 5, 'cl': 7, 'ca': 6}
CLOUDY_YIELD = Path.home() / 'CLOUDY/c25.00/data/mewe_nelectron.dat'
CDS_TABLE2 = (Path.home() / 'RT_Codes/ExoAtmosphere/references/'
              'kaastra_mewe_1993/table2')
OUTPUT = ROOT / 'md/atomic_stage_expansion/validation_2026-10-10/auger_audit'
SHELLS = ('1s', '2s', '2p', '3s', '3p', '3d', '4s')
# KM93 resolves p orbitals into two subshells; PH1 and Cloudy combine each pair.
CDS_SHELLS = {1: (1,), 2: (2,), 3: (3, 4), 4: (5,), 5: (6, 7)}


def read_yields(path):
    rows = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        a = line.split()
        if len(a) != 14:
            raise ValueError(f'Invalid Cloudy yield row: {line}')
        z, stage, shell, maximum = map(int, a[:4])
        probabilities = [float(x) for x in a[4:]]
        if abs(sum(probabilities) - 1.0) > 1e-3:
            raise ValueError(f'Yield normalization failed: {(z, stage, shell)}')
        if any(p > 0 for p in probabilities[maximum:]):
            raise ValueError(f'Yield exceeds declared maximum: {(z, stage, shell)}')
        rows[(z, stage, shell)] = probabilities
    return rows


def read_cds_table2(path):
    rows = {}
    for line in path.read_text().splitlines():
        if len(line) != 76:
            raise ValueError(f'Invalid CDS Table 2 row: {line}')
        z, stage, shell = int(line[0:2]), int(line[2:5]), int(line[5:7])
        probabilities = [int(line[26+5*i:31+5*i])/10000 for i in range(10)]
        rows[(z, stage, shell)] = {'edge_eV': float(line[7:14]),
                                    'mean_auger_eV': float(line[14:21]),
                                    'epsilon': int(line[21:25])/1000,
                                    'probabilities': probabilities}
    return rows


def audit(upper):
    ph1, ninn, ntot = source.parse_ph1_phfit2()
    yields = read_yields(CLOUDY_YIELD)
    cds = read_cds_table2(CDS_TABLE2)
    current = {e: n for e, _, n in source.ELEMENTS}
    channels = []
    shell_rows = []
    missing_yields = []
    for element, z, _ in source.ELEMENTS:
        for stage in range(1, TARGET[element]):
            ne = z - stage + 1
            einn, _ = source.phfit2_shells(ph1, ninn, ntot, z, ne)
            outer = source.ETH[(z, ne)]
            nout = ntot[ne-1]
            if z == ne and z > 18:
                nout = 7
            if z == ne + 1 and z in (20, 21, 22, 25, 26):
                nout = 7
            occupied = [s for s in range(1, nout+1)
                        if (z, ne, s) in ph1 and ph1[(z, ne, s)][2] > 0]
            for shell in occupied:
                edge = ph1[(z, ne, shell)][0]
                if shell == occupied[-1]:
                    edge = outer
                activation = max(outer, edge, einn)
                if activation >= upper:
                    continue
                # The PH2 outer fit is used below E_inn; PH1 shell rows
                # are physically scored only once the full shell sum applies.
                p = yields.get((z, stage, shell))
                if p is None:
                    missing_yields.append([element, stage, shell])
                    continue
                is_inner = shell != occupied[-1]
                shell_rows.append({'element': element, 'stage': stage,
                                   'shell': SHELLS[shell-1], 'edge_eV': edge,
                                   'activation_eV': activation, 'inner': is_inner,
                                   'current_stage': stage < current[element]})
                for electrons, probability in enumerate(p, start=1):
                    if probability <= 0:
                        continue
                    destination = stage + electrons
                    minimum = sum(source.ETH[(z, z - j + 1)]
                                  for j in range(stage, destination)) if destination <= z+1 else float('inf')
                    row = {'element': element, 'source_stage': stage,
                           'shell': SHELLS[shell-1], 'shell_index': shell,
                           'shell_edge_eV': edge, 'activation_eV': activation,
                           'ejected_electrons': electrons,
                           'yield': probability, 'destination_stage': destination,
                           'represented_in_target': destination <= TARGET[element],
                           'minimum_sequential_IP_eV': minimum,
                           'energetically_possible_below_upper': max(activation, minimum) < upper,
                           'channel_onset_eV': max(activation, minimum),
                           'source_in_current_roster': stage < current[element],
                           'inner': is_inner}
                    if is_inner and electrons > 1 and shell in CDS_SHELLS:
                        keys = [(z, stage, s) for s in CDS_SHELLS[shell]]
                        present = [key for key in keys if key in cds]
                        row['cds_table2_shells'] = list(CDS_SHELLS[shell])
                        row['cds_table2_probabilities'] = [
                            cds[key]['probabilities'][electrons-1] for key in present]
                        row['cds_table2_edges_eV'] = [cds[key]['edge_eV'] for key in present]
                        row['cds_mean_auger_eV'] = [cds[key]['mean_auger_eV'] for key in present]
                        row['cds_energy_discrepancy'] = [cds[key]['epsilon'] for key in present]
                        row['cds_table2_mapping_complete'] = len(present) == len(keys)
                        if len(present) == len(keys) and any(
                                abs(cds[key]['probabilities'][electrons-1] - probability) > 2e-5
                                for key in keys):
                            raise ValueError(f'Cloudy and CDS yields differ for {(z, stage, shell, electrons)}')
                    channels.append(row)
    relevant = [c for c in channels if c['inner'] and c['ejected_electrons'] > 1]
    potential = [c for c in relevant if c['energetically_possible_below_upper']]
    return {'upper_energy_eV': upper,
            'criterion': 'PH1 active at max(stage IP, shell edge, PHINNER) < upper; positive yield; sequential final-stage IP sum < upper',
            'counts': {'ph1_shells_active': len(shell_rows),
                       'positive_inner_nonadjacent_channels_raw': len(relevant),
                       'energetically_possible_inner_nonadjacent_channels': len(potential),
                       'energetically_possible_unrepresented_products': sum(not c['represented_in_target'] for c in potential),
                       'by_element': dict(sorted(Counter(c['element'] for c in potential).items()))},
            'missing_yields': missing_yields, 'channels': channels}


def plot_counts(data, output):
    cache_dir = tempfile.TemporaryDirectory(prefix='auger_mpl_')
    os.environ.setdefault('MPLCONFIGDIR', cache_dir.name)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels = list(TARGET)
    counts = data['counts']['by_element']
    fig, ax = plt.subplots(figsize=(9, 3.6))
    ax.bar(labels, [counts.get(x, 0) for x in labels], color='#325a83')
    ax.set_xlabel('Element')
    ax.set_ylabel('Reachable inner-shell channels')
    ax.set_title('Positive nonadjacent yields below 100 eV')
    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)
    cache_dir.cleanup()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upper-eV', type=float, default=100.0)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    args = parser.parse_args()
    result = audit(args.upper_eV)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / 'auger_channels.json').write_text(json.dumps(result, indent=2) + '\n')
    plot_counts(result, args.output_dir / 'auger_channels.pdf')
    print(json.dumps(result['counts'], indent=2))
    if result['missing_yields']:
        raise SystemExit(f'Missing {len(result["missing_yields"])} yield rows')


if __name__ == '__main__':
    main()
