"""Correct Mg III E1 rates in a preserved private atom; do not invoke the writer."""

import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'tools/fitting'))
import chianti_cooling as C
import write_excitation_atom as W

AUDIT = HERE.parent/'mgiii_radiative_consistency/result.json'
def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    audit = json.loads(AUDIT.read_text())
    original, _, _ = C.read_excitation_atom('mg', 3)
    levels, radiation, collisions = C.bound_excitation_atom('mg', 3)
    if (len(levels), len(radiation), len(collisions)) != (167, 3364, 13861):
        raise ValueError('Full bound Mg III source counts changed')
    original_hashes = {p: file_hash(p) for p in audit['source_sha256']}
    if original_hashes != audit['source_sha256']:
        raise ValueError('Mg III raw source hash changed after prior audit')
    source_rows = {(r['source_lower'], r['source_upper']): r for r in audit['rows']}
    source_ids = sorted(levels)
    compact = {source_id: index for index, source_id in enumerate(source_ids, 1)}
    # The production writer must reject a reversed source A. This preserved
    # diagnostic atom predates that guard and is corrected below in isolation.
    base_path = HERE/'original_A/nlevel_mg_3.txt'
    if not base_path.is_file():
        raise FileNotFoundError('Preserved private raw Mg III atom is required')
    text = base_path.read_text()
    lines = text.splitlines()
    source_id_line = next(line for line in lines
                          if line.startswith('# source level IDs:'))
    baseline_ids = [int(x) for x in source_id_line.split(':', 1)[1].split()]
    if baseline_ids != source_ids:
        raise ValueError('Preserved atom selected different source levels')
    total_nodes = int(next(line.split()[1] for line in lines
                           if line.startswith('NODES ')))
    nrad_header = next(i for i, line in enumerate(lines) if line.startswith('NRAD '))
    nrad = int(lines[nrad_header].split()[1])
    original_rates = {tuple(map(int, line.split()[:2])): float(line.split()[2])
                      for line in lines[nrad_header+1:nrad_header+1+nrad]}
    if len(original_rates) != nrad or nrad != len(radiation):
        raise ValueError('Unexpected original radiative row count')
    corrected = {}
    row_records = []
    for source_lower, source_upper, wavelength, gf, rate in radiation:
        if rate <= 0:
            continue
        details = source_rows[(source_lower, source_upper)]
        if source_lower not in compact or source_upper not in compact:
            raise ValueError('Unexpected inactive radiative row')
        first_energy = levels[source_lower][1]
        second_energy = levels[source_upper][1]
        lower, upper = ((source_lower, source_upper) if first_energy < second_energy
                        else (source_upper, source_lower))
        gap_cm1 = levels[upper][1] - levels[lower][1]
        if gap_cm1 <= 0:
            raise ValueError('Nonpositive emitted photon energy')
        pair = (compact[lower], compact[upper])
        original_rate = float(rate)
        if max(source_lower, source_upper) <= 27:
            classification = 'Jonsson_2014_original'
            adopted_rate = original_rate
        elif gf > 0 and details['dipole_selection']:
            # Source gf is an oscillator strength, proportional to frequency
            # times the dipole line strength. Preserve A/nu^3 from the source
            # theoretical transition and use the physical emitting weight.
            classification = 'Liang_Badnell_2010_E1_dipole_strength'
            theory_gap = details['theoretical_gap_cm-1']
            if theory_gap <= 0:
                raise ValueError('Nonpositive theoretical E1 transition energy')
            adopted_rate = (original_rate * (gap_cm1/theory_gap)**3 *
                            original[source_upper][0]/original[upper][0])
        else:
            # The source contains no multipole type or line strength for these
            # rows. Retaining A is an explicit effective-rate approximation.
            classification = 'higher_level_unidentified_multipole_original_A'
            adopted_rate = original_rate
        if not math.isfinite(adopted_rate) or adopted_rate <= 0:
            raise ValueError('Invalid dipole-strength-corrected spontaneous rate')
        if pair in corrected:
            raise ValueError('Duplicate Mg III physical radiative pair')
        corrected[pair] = adopted_rate
        row_records.append({
            'source_lower': source_lower, 'source_upper': source_upper,
            'physical_lower': lower, 'physical_upper': upper,
            'classification': classification,
            'selected_gap_cm1': gap_cm1,
            'theory_gap_cm1': details['theoretical_gap_cm-1'],
            'source_A_s-1': original_rate,
            'candidate_A_s-1': adopted_rate,
            'candidate_over_source_A': adopted_rate/original_rate,
            'source_gf': gf,
            'source_wavelength_A': wavelength,
        })
    if set(corrected) != set(original_rates):
        raise ValueError('Physical radiative pairs differ from source writer')
    if max(abs(original_rates[pair]/row['source_A_s-1']-1)
           for pair, row in zip(corrected, row_records)) > 1e-12:
        raise ValueError('Original writer radiative rate differs from source')
    candidate_lines = lines[:]
    for i in range(nrad):
        first, second, _ = candidate_lines[nrad_header+1+i].split()
        pair = int(first), int(second)
        candidate_lines[nrad_header+1+i] = (
            f'{first} {second} {corrected[pair]:.17g}')
    candidate_lines.insert(0, '# PRIVATE MG III RADIATIVE CANDIDATE; no publication or production acceptance.')
    candidate_lines.insert(1, '# Higher-level E1 rows preserve source dipole strength through A proportional to photon frequency cubed.')
    candidate_lines.insert(2, '# Higher-level gf=0 rates retain source A as an unidentified-multipole approximation.')
    candidate_lines.insert(3, '# First 27 Jonsson (2014) spontaneous rates retain their source values.')
    candidate_path = HERE/'nlevel_mg_3_line_strength_candidate.txt'
    candidate_path.write_text('\n'.join(candidate_lines)+'\n')
    if original_hashes != {p: file_hash(p) for p in audit['source_sha256']}:
        raise ValueError('Mg III raw source changed during candidate generation')
    kinds = sorted(set(r['classification'] for r in row_records))
    groups = {kind: [r for r in row_records if r['classification'] == kind]
              for kind in kinds}
    ratios = np.array([r['candidate_over_source_A'] for r in
                       groups['Liang_Badnell_2010_E1_dipole_strength']])
    ratio_quantiles = {str(q): float(np.quantile(ratios, q)) for q in
                       [0, .01, .1, .5, .9, .99, 1]}
    original_sum = {}
    corrected_sum = {}
    unidentified_sum = {}
    for row in row_records:
        upper = row['physical_upper']
        original_sum[upper] = original_sum.get(upper, 0.) + row['source_A_s-1']
        corrected_sum[upper] = corrected_sum.get(upper, 0.) + row['candidate_A_s-1']
        if row['classification'] == 'higher_level_unidentified_multipole_original_A':
            unidentified_sum[upper] = (unidentified_sum.get(upper, 0.) +
                                       row['candidate_A_s-1'])
    ratio_total = np.array([corrected_sum[i]/original_sum[i] for i in corrected_sum])
    multipole_fraction = {str(i): unidentified_sum.get(i, 0.)/corrected_sum[i]
                          for i in corrected_sum if unidentified_sum.get(i, 0.) > 0}
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    axes[0].scatter([r['selected_gap_cm1'] for r in groups['Liang_Badnell_2010_E1_dipole_strength']],
                    ratios, s=7, alpha=.5)
    axes[0].set_xscale('log'); axes[0].set_yscale('log')
    axes[0].set_xlabel(r'Selected photon gap [cm$^{-1}$]')
    axes[0].set_ylabel(r'Candidate $A$ / source $A$')
    axes[1].scatter([int(i) for i in multipole_fraction], list(multipole_fraction.values()),
                    s=9, alpha=.6)
    axes[1].set_xlabel('Physical emitting source level')
    axes[1].set_ylabel('Unidentified-multipole fraction of total A')
    axes[1].set_ylim(0, 1.01)
    fig.savefig(HERE/'radiative_correction_diagnostic.pdf')
    plt.close(fig)
    result = {
        'status': 'private radiative candidate; no finite-density validation or production acceptance',
        'source_hashes': original_hashes,
        'candidate_builder_hash': file_hash(Path(__file__)),
        'selection_code_hashes': {str(p): file_hash(p) for p in
                                  W._selection_code_files('mg', 3)},
        'prior_audit_hash': file_hash(AUDIT),
        'baseline_atom': str(base_path), 'baseline_hash': file_hash(base_path),
        'candidate_atom': str(candidate_path), 'candidate_hash': file_hash(candidate_path),
        'levels': len(levels), 'radiative_rows': len(row_records),
        'collision_rows': len(collisions), 'collision_nodes': total_nodes,
        'radiative_class_counts': {key: len(value) for key, value in groups.items()},
        'source_reversed_E1_rows': sum(r['selected_gap_cm1'] > 0 and
                                       levels[r['source_upper']][1] < levels[r['source_lower']][1]
                                       for r in groups['Liang_Badnell_2010_E1_dipole_strength']),
        'E1_rate_ratio_quantiles': ratio_quantiles,
        'E1_rates_outside_factor_two': int(np.sum((ratios < .5) | (ratios > 2))),
        'emitting_level_total_A_ratio_quantiles': {
            str(q): float(np.quantile(ratio_total, q)) for q in [0, .01, .5, .99, 1]},
        'unidentified_multipole_emitting_levels': len(multipole_fraction),
        'unidentified_multipole_branch_fraction_max': max(multipole_fraction.values()),
        'unidentified_multipole_branch_fraction_above_1_percent':
            sum(value > .01 for value in multipole_fraction.values()),
        'largest_rate_corrections': sorted(
            [r for r in groups['Liang_Badnell_2010_E1_dipole_strength']],
            key=lambda r: abs(math.log(r['candidate_over_source_A'])), reverse=True)[:20],
    }
    (HERE/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in
                      ['status','levels','radiative_rows','collision_rows',
                       'collision_nodes','radiative_class_counts',
                       'source_reversed_E1_rows','E1_rate_ratio_quantiles',
                       'E1_rates_outside_factor_two',
                       'unidentified_multipole_branch_fraction_max']}, indent=2))


if __name__ == '__main__':
    main()
