"""Measure current Ar VI source cooling against forty selected levels."""

import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT/'tools/fitting'))
import chianti_cooling as C
from validate_excitation_density_profiles import full_excitation_balance
from validate_source_excitation_atom import DENSITY, source_subset
from write_excitation_atom import _selection_code_files, _source_files


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--element', choices=['ar'], required=True)
    parser.add_argument('--stage', type=int, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    if (args.element, args.stage) not in {('ar', 6)}:
        parser.error('Only Ar VI is in this source audit')
    resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    temperature = np.geomspace(1e3, 1e5, 801)
    levels, radiation, collisions = C.bound_excitation_atom(args.element, args.stage)
    ids = sorted(levels)[:40]
    selected, selected_radiation, selected_collisions = source_subset(
        args.element, args.stage, ids)
    pchip, zero_node, negative_node = [], [], []
    for transition in selected_collisions:
        pair = C.excitation_levels(selected, transition)
        if isinstance(C.scaled_collision_strength_interpolator(transition),
                      PchipInterpolator):
            pchip.append(pair)
        values = np.asarray(transition['ys'])
        if np.any(values == 0):
            zero_node.append(pair)
        if np.any(values < 0):
            negative_node.append(pair)
    full = full_excitation_balance(levels, radiation, collisions,
                                   temperature, DENSITY)
    chosen = full_excitation_balance(selected, selected_radiation,
                                     selected_collisions, temperature, DENSITY)
    saved = {'temperature_K': temperature}
    cutoff, zero_counts, minimum_positive, loss_counts = [], [], [], []
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True,
                                constrained_layout=True)
    differences, difference_axes = plt.subplots(2, 2, figsize=(9, 7),
                                               sharex=True, constrained_layout=True)
    for index, density in enumerate(DENSITY):
        reference, full_population = full[density]
        selected_cooling, selected_population = chosen[density]
        if (not np.all(np.isfinite(reference)) or np.any(reference < 0)
            or not np.all(np.isfinite(selected_cooling)) or np.any(selected_cooling < 0)):
            raise ValueError('Invalid source cooling')
        normal = reference > 64*np.finfo(float).tiny
        if not np.any(normal):
            raise ValueError('No normal source cooling in the numerical use interval')
        error = np.abs(selected_cooling[normal]/reference[normal]-1)
        cutoff.append(float(np.max(error)))
        zero_counts.append(dict(full=int(np.sum(reference == 0)),
                                selected=int(np.sum(selected_cooling == 0)),
                                normal=int(np.sum(normal))))
        loss_counts.append(int(np.sum((reference > 0) & (selected_cooling == 0))))
        minimum_positive.append(float(np.min(reference[reference > 0])))
        saved[f'cooling_density_{index}'] = reference
        saved[f'cooling_40_density_{index}'] = selected_cooling
        saved[f'populations_full_density_{index}'] = full_population
        saved[f'populations_40_density_{index}'] = selected_population
        label = r'$n_e = ' + f'{density:g}' + r'$ cm$^{-3}$'
        axis = axes.flat[index]
        axis.loglog(temperature, reference, label='Full bound source')
        axis.loglog(temperature, selected_cooling, '--', label='40 selected states')
        axis.set_title(label)
        axis.set_xlabel(r'$T$ [K]')
        axis.set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
        axis.legend(fontsize=8)
        axis = difference_axes.flat[index]
        axis.semilogx(temperature[normal],
                      100*(selected_cooling[normal]/reference[normal]-1))
        axis.set_title(label)
        axis.set_xlabel(r'$T$ [K]')
        axis.set_ylabel(r'Difference [\%]')
    figure.savefig(work/'cooling_profiles.pdf')
    differences.savefig(work/'cutoff_differences.pdf')
    plt.close(figure)
    plt.close(differences)
    np.savez_compressed(work/'full_source_profiles.npz', **saved)
    dependencies = list(_source_files(args.element, args.stage))
    dependencies += list(_selection_code_files(args.element, args.stage))
    dependencies += [Path(__file__), ROOT/'tools/fitting/chianti_cooling.py',
                     ROOT/'tools/fitting/validate_excitation_density_profiles.py']
    result = dict(passed=True, element=args.element, stage=args.stage,
                  full_levels=len(levels), selected_levels=40,
                  full_collision_rows=len(collisions),
                  selected_collision_rows=len(selected_collisions),
                  affected_pchip_pairs=pchip, selected_zero_node_pairs=zero_node,
                  selected_negative_node_pairs=negative_node,
                  temperature_points=len(temperature),
                  temperature_K=[float(temperature[0]), float(temperature[-1])],
                  densities_cm3=DENSITY, maximum_relative_cutoff_errors=cutoff,
                  source_zero_counts=zero_counts,
                  positive_full_selected_zero_counts=loss_counts,
                  minimum_positive_full_cooling=minimum_positive,
                  selected_40_meets_cutoff=max(cutoff) <= 0.003,
                  sources={str(path): digest(path) for path in dependencies},
                  profiles_sha256=digest(work/'full_source_profiles.npz'))
    (work/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
