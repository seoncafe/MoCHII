"""Diagnose one omitted charge state using original RR/DR temperature nodes.

The analytic upper bound applies only to an adjacent stationary electron-only
chain with one added state. It excludes photoionization, charge exchange, Auger
transitions, excited-state ionization, and all further omitted charge states.
It is not an acceptance criterion for the production 64-stage truncation.
"""

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'tools/fitting'))
import make_element_data as sources
import validate_terminal_ionization as ionization


def main():
    output = HERE / 'native_node_bounds.json'
    if output.exists():
        raise FileExistsError(f'Preserve the existing diagnostic: {output}')
    provenance_path = HERE.parent / 'terminal_reverse_provenance/rate_provenance.py'
    spec = importlib.util.spec_from_file_location('terminal_rate_provenance', provenance_path)
    provenance = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(provenance)
    originals = {kind: provenance.source_comparison(path)
                 for kind, path in provenance.MASTER.items()}
    rr, dr, ci = sources.parse_badnell_rr(), sources.parse_badnell_dr(), sources.parse_cf()
    channels = []
    for element, (z, terminal) in ionization.ROSTER.items():
        if element == 'fe':
            channels.append({'element': element, 'status': 'not_evaluated',
                             'reason': 'Fe VIII uses a separate experimental DR source; '
                                       'the Badnell original-node audit does not apply.'})
            continue
        key = (z, z - terminal, 1)
        rr_nodes, dr_nodes = (np.asarray(originals[kind][key]) for kind in ('RR', 'DR'))
        temperatures = np.intersect1d(rr_nodes[:, 0], dr_nodes[:, 0])
        temperatures = temperatures[(temperatures >= 1000) & (temperatures <= 100000)]
        assert len(temperatures) > 0
        rr_nodes = rr_nodes[np.isin(rr_nodes[:, 0], temperatures)]
        dr_nodes = dr_nodes[np.isin(dr_nodes[:, 0], temperatures)]
        assert np.array_equal(rr_nodes[:, 0], temperatures)
        assert np.array_equal(dr_nodes[:, 0], temperatures)
        assert np.all(np.isfinite(rr_nodes)) and np.all(np.isfinite(dr_nodes))
        assert np.all(rr_nodes[:, 1] > 0) and np.all(dr_nodes[:, 1] >= 0)
        original = rr_nodes[:, 1] + dr_nodes[:, 1]
        printed_fit = rr_nodes[:, 2] + dr_nodes[:, 2]
        row, rr_origin, dr_origin = ionization.source_terminal(element, z, terminal, rr, dr, ci)
        assert rr_origin == dr_origin == 'Badnell M=1'
        evaluated_fit = np.asarray([ionization.recombination(row, t) for t in temperatures])
        assert np.all(evaluated_fit > 0) and np.all(printed_fit > 0)
        log_ci = ionization.ionization_log(row, temperatures)
        log_ratio = log_ci - np.log(original)
        # r/(1+r) is stable even when an extremely cold CI coefficient underflows.
        bound = np.exp(-np.logaddexp(0, -log_ratio))
        fit_bound = np.exp(-np.logaddexp(0, -(log_ci - np.log(evaluated_fit))))
        assert np.all(np.isfinite(bound)) and np.all((bound >= 0) & (bound <= 1))
        residual = np.abs(evaluated_fit / original - 1)
        printed_residual = np.abs(printed_fit / original - 1)
        peak = int(np.argmax(bound))
        channels.append({
            'element': element, 'terminal_stage': terminal, 'diagnostic_next_stage': terminal + 1,
            'source_key_Z_N_M': list(key), 'CI_key_Z_N': [z, z - terminal + 1],
            'terminal_CI_model': 'Voronov source fit without a calibrated terminal HYBRID factor',
            'temperature_K': temperatures.tolist(),
            'RR_original_cm3_s': rr_nodes[:, 1].tolist(),
            'DR_original_cm3_s': dr_nodes[:, 1].tolist(),
            'reverse_original_cm3_s': original.tolist(),
            'reverse_evaluated_fit_cm3_s': evaluated_fit.tolist(),
            'log_CI_over_original_reverse': log_ratio.tolist(),
            'one_added_state_fraction_upper_bound_original': bound.tolist(),
            'one_added_state_fraction_upper_bound_fit': fit_bound.tolist(),
            'maximum_original_bound': float(bound[peak]),
            'maximum_original_bound_temperature_K': float(temperatures[peak]),
            'maximum_total_reverse_fit_residual': float(np.max(residual)),
            'maximum_printed_total_reverse_fit_residual': float(np.max(printed_residual)),
            'residual_interpretation': 'Fit discrepancy against rounded original nodes; '
                                       'not an atomic-physics uncertainty estimate.',
        })
    paths = [Path(__file__), provenance_path, Path(sources.__file__), Path(ionization.__file__),
             Path(sources.BADNELL_RR), Path(sources.BADNELL_DR),
             Path(sources.VERNER_DIR) / 'cfit.f', *provenance.MASTER.values()]
    result = {
        'production_stage_count': 64, 'production_adjacent_transition_count': 53,
        'production_modified': False,
        'model': 'Adjacent stationary electron-only chain with exactly one added state; '
                 'if r=C_terminal/alpha_reverse, its fraction is <=r/(1+r).',
        'derivation': 'Normalize the retained chain to sum S with terminal weight w<=S. '
                      'The added state has weight r*w and fraction r*w/(S+r*w).',
        'scope_limits': ['No bound on all further charge states.', 'No charge exchange or photons.',
                         'No thermal-error bound.', 'No physical uncertainty bound from fit residuals.',
                         'Native nodes only; no interpolation of original RR/DR columns.',
                         'No production temperature limit or tolerance is selected.'],
        'channels': channels,
        'source_sha256': {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
    }
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({r['element']: r.get('maximum_original_bound') for r in channels}))


if __name__ == '__main__':
    main()
