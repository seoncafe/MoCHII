#!/usr/bin/env python3
"""Evaluate KM93 Eq. (7) for the audited vacancies; do not install heat data.

The correction uses the paper's own Lotz vacancy energy and tabulated epsilon.
It is an approximate unconditional mean, not an electron spectrum or a
correction to the PH1 threshold. Table precision is not physical uncertainty.
"""
import argparse
import hashlib
import json
from pathlib import Path

from audit_auger_channels import CDS_TABLE2, ROOT, read_cds_table2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise RuntimeError('Preserve existing audit products')
    audit_path = ROOT / 'md/atomic_stage_expansion/validation_2026-10-10/auger_audit_cds_checked/auger_channels.json'
    final_article = Path.home() / 'MoCHII/references/Kaastra_Mewe_1993_AAS_97_443.pdf'
    channels = json.loads(audit_path.read_text())['channels']
    table = read_cds_table2(CDS_TABLE2)
    atomic_numbers = {'mg': 12, 'fe': 26, 'ca': 20}
    rows = {}
    selected = []
    for channel in channels:
        if not (channel['inner'] and channel['ejected_electrons'] > 1
                and channel['energetically_possible_below_upper']):
            continue
        selected.append(channel)
        for shell in channel['cds_table2_shells']:
            key = (atomic_numbers[channel['element']], channel['source_stage'], shell)
            source = table[key]
            energy = source['edge_eV']
            mean = source['mean_auger_eV']
            epsilon = source['epsilon']
            corrected = mean - epsilon * energy
            # Conservative storage-rounding bound: EA and I have one decimal,
            # epsilon has three. This is NOT an atomic-model error estimate.
            storage_bound = 0.05 + abs(epsilon)*0.05 + abs(energy)*0.0005 + 0.05*0.0005
            rows[key] = {
                'element': channel['element'], 'source_stage': key[1],
                'KM_shell_index': shell, 'PH1_shell': channel['shell'],
                'PH1_edge_eV': channel['shell_edge_eV'],
                'KM_vacancy_energy_eV': energy, 'KM_mean_auger_eV': mean,
                'KM_epsilon': epsilon, 'KM_equation7_mean_auger_eV': corrected,
                'storage_rounding_bound_eV': storage_bound,
                'negative_beyond_storage_rounding': corrected + storage_bound < 0,
                'probability_sum_at_stored_precision': sum(source['probabilities']),
            }
    assert len(selected) == 9 and len(rows) == 14
    negative = [r for r in rows.values() if r['negative_beyond_storage_rounding']]
    assert len(negative) == 1
    assert (negative[0]['element'], negative[0]['source_stage'],
            negative[0]['KM_shell_index']) == ('ca', 2, 5)
    sources = (Path(__file__).resolve(), ROOT/'tools/fitting/audit_auger_channels.py',
               audit_path, CDS_TABLE2, final_article)
    result = {
        'scope': 'Source Eq. (7) arithmetic for nine historical reachable PH1 channels; no installed rates, heat, or ion fractions',
        'published_source': 'Kaastra and Mewe 1993, A&AS 97, 443-482, journal pages 445-447, Sections 2.2 and 3',
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        'physical_PH1_channels': len(selected), 'KM_subshell_rows': len(rows),
        'negative_corrected_means': len(negative),
        'production_modified': False, 'rows': list(rows.values()),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('Evaluated 14 KM subshell means for 9 PH1 channels; Ca II M1 corrected mean is negative.')


if __name__ == '__main__':
    main()
