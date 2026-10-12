"""Evaluate charge-only Auger branches and rigorous event-energy bounds."""

import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'tools/fitting'))
import make_element_data as physical_ionization  # noqa: E402


SOURCE_AUDIT = ROOT / 'md/atomic_stage_expansion/validation_2026-10-10/auger_audit_cds_checked/auger_channels.json'
CDS = Path('/home/kiseon/RT_Codes/ExoAtmosphere/references/kaastra_mewe_1993/table2')
OUTPUT = Path(__file__).with_name('charge_energy_bounds.json')
TARGET = {'c': 5, 'n': 6, 'o': 5, 'ne': 5, 's': 7, 'ar': 7,
          'mg': 4, 'fe': 7, 'si': 5, 'cl': 7, 'ca': 6}
ATOMIC_NUMBER = {'mg': 12, 'fe': 26, 'ca': 20}
CDS_SHELL = {'2s': [2], '2p': [3, 4], '3s': [5], '3p': [6, 7]}
PHYSICAL_SHELL = {'2s': 2, '2p': 3, '3s': 4, '3p': 5}


def cds_table2():
    rows = {}
    for line in CDS.read_text().splitlines():
        assert len(line) == 76
        key = (int(line[:2]), int(line[2:5]), int(line[5:7]))
        p = [int(line[26 + 5 * i:31 + 5 * i]) / 10000 for i in range(10)]
        rows[key] = {'vacancy_eV': float(line[7:14]),
                     'mean_Auger_eV': float(line[14:21]), 'electron_count': p}
    return rows


def stage_indices():
    assert sum(TARGET.values()) == 64
    return {(element, stage): 1 + sum(list(TARGET.values())[:i]) + stage - 1
            for i, (element, count) in enumerate(TARGET.items())
            for stage in range(1, count + 1)}


def budget_at_energy(energy, ph1_edge, binding_two, km_min, km_max,
                     incompatible_primary, probabilities_equal):
    if energy < ph1_edge:
        status = 'below_PH1_absorption'
    elif energy < binding_two:
        status = 'final_charge_two_forbidden'
    elif incompatible_primary:
        status = 'PH1_primary_exceeds_final_charge_energy_budget'
    elif energy < km_max:
        status = 'KM_vacancy_unavailable_or_unresolved'
    elif not probabilities_equal:
        status = 'KM_subshell_charge_probabilities_differ'
    else:
        status = 'KM_charge_marginal_lookup_only'
    available = energy - binding_two
    primary = energy - ph1_edge
    return {
        'photon_eV': energy,
        'status': status,
        'total_emitted_electron_kinetic_energy_interval_eV':
            [0.0, available] if available >= 0 else None,
        'PH1_primary_kinetic_energy_eV': primary if primary >= 0 else None,
        'conditional_total_electron_interval_if_PH1_primary_eV':
            [primary, available] if primary >= 0 and available >= primary else None,
        'KM_vacancy_minimum_eV': km_min,
        'KM_vacancy_maximum_eV': km_max,
        'PH1_primary_leaves_nonnegative_Auger_photon_excitation_budget':
            primary <= available if primary >= 0 and available >= 0 else None,
    }


def main():
    source = json.loads(SOURCE_AUDIT.read_text())
    cds = cds_table2()
    global_index = stage_indices()
    rows = []
    at_100_generators = []
    for channel in source['channels']:
        if channel.get('ejected_electrons') != 2 or not channel.get('inner'):
            continue
        element = channel['element']
        stage = channel['source_stage']
        shell = channel['shell']
        z = ATOMIC_NUMBER[element]
        physical_shell = PHYSICAL_SHELL[shell]
        assert physical_shell == channel['shell_index']
        assert channel['destination_stage'] == stage + 2
        assert stage + 2 <= TARGET[element]
        km_keys = [(z, stage, s) for s in CDS_SHELL[shell]]
        km_rows = [cds[key] for key in km_keys]
        probabilities = [r['electron_count'] for r in km_rows]
        probabilities_equal = all(p == probabilities[0] for p in probabilities[1:])
        p1, p2 = probabilities[0][:2]
        assert probabilities_equal and sum(probabilities[0][2:]) == 0
        assert abs(p1 + p2 - 1) <= 1e-12
        assert abs(channel['cds_table2_probabilities'][0] - p2) <= 1e-12
        ph1_edge = channel['shell_edge_eV']
        binding_one = physical_ionization.ETH[(z, z - stage + 1)]
        binding_two = channel['minimum_sequential_IP_eV']
        next_binding = physical_ionization.ETH[(z, z - stage)]
        assert abs(binding_one + next_binding - binding_two) < 1e-10
        km_min = min(r['vacancy_eV'] for r in km_rows)
        km_max = max(r['vacancy_eV'] for r in km_rows)
        incompatible_primary = ph1_edge < binding_two
        nodes = sorted(set([ph1_edge, binding_two, km_min, km_max, 100.0]))
        budget = [budget_at_energy(e, ph1_edge, binding_two, km_min, km_max,
                                   incompatible_primary, probabilities_equal)
                  for e in nodes if e <= 100.0]
        supported = budget[-1]['status'] == 'KM_charge_marginal_lookup_only'
        entry = {
            'element': element, 'source_stage': stage,
            'physical_PH1_shell_index': physical_shell, 'shell_label': shell,
            'KM_shell_indices': CDS_SHELL[shell],
            'global_source_index_1based': global_index[(element, stage)],
            'global_plus_one_index_1based': global_index[(element, stage + 1)],
            'global_plus_two_index_1based': global_index[(element, stage + 2)],
            'PH1_edge_eV': ph1_edge,
            'KM_vacancy_edges_eV': [r['vacancy_eV'] for r in km_rows],
            'ground_state_binding_plus_one_eV': binding_one,
            'ground_state_binding_plus_two_eV': binding_two,
            'KM_charge_one_probability': p1,
            'KM_charge_two_probability': p2,
            'KM_subshell_probabilities_identical': probabilities_equal,
            'PH1_primary_residual_after_final_ground_binding_eV':
                ph1_edge - binding_two,
            'invalid_if_KM_mean_Auger_assigned_to_plus_two_eV':
                [binding_two + r['mean_Auger_eV'] - ph1_edge for r in km_rows],
            'energy_budget_at_boundaries': budget,
            'charge_marginal_lookup_possible_at_100eV': supported,
            'charge_marginal_not_a_thermal_model': True,
        }
        if supported:
            entry['unit_shell_absorption_generator_entries_1based'] = [
                [entry['global_source_index_1based'], entry['global_source_index_1based'], -1.0],
                [entry['global_plus_one_index_1based'], entry['global_source_index_1based'], p1],
                [entry['global_plus_two_index_1based'], entry['global_source_index_1based'], p2],
            ]
            assert abs(sum(x[2] for x in entry['unit_shell_absorption_generator_entries_1based'])) < 1e-12
            at_100_generators.append(entry['unit_shell_absorption_generator_entries_1based'])
        rows.append(entry)
    assert len(rows) == 9 and len(at_100_generators) == 6
    result = {
        'scope': 'read-only charge-marginal and event-energy bounds, no heat or ion fractions',
        'number_of_stages': 64,
        'number_of_adjacent_transitions': sum(TARGET.values()) - len(TARGET),
        'stage_order': list(TARGET),
        'source_sha256': {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [SOURCE_AUDIT, CDS, ROOT / 'tools/fitting/make_element_data.py']},
        'channels': rows,
        'unit_shell_generators_available_at_100eV': len(at_100_generators),
        'unit_shell_generators_unresolved_at_100eV': len(rows) - len(at_100_generators),
        'probability_use': 'CDS Table 2 charge marginals only at or above both PH1 and all mapped KM vacancy edges, with a nonnegative PH1 primary-to-final-ground budget; this is not a validated energy-dependent yield',
        'thermal_energy_computed': False,
    }
    assert result['number_of_adjacent_transitions'] == 53
    OUTPUT.write_text(json.dumps(result, indent=2) + '\n')
    print(OUTPUT)


if __name__ == '__main__':
    main()
