"""Compare selected bound excitation atoms over 100--1,000,000 K.

The diagnostic measures source-level cutoff error. It does not validate the
scattering calculation, omitted processes, or the native serialized atom.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

import chianti_cooling as C
from atomic_output import atomic_write_text
from validate_excitation_density_profiles import full_excitation_balance


IONS = {
    'ArV': ('ar', 5, 110),
    'MgIV': ('mg', 4, 110),
    'FeVII': ('fe', 7, 170),
    'NeV': ('ne', 5, 25),
}
DENSITIES = (1e-6, 100.0, 1e6, 1e10)
TEMPERATURE = np.geomspace(100.0, 1e6, 401)
TARGET = 0.003
FLOOR = 64.0 * np.finfo(float).tiny


def source_paths(element: str, stage: int) -> list[Path]:
    stem = Path(C.ion_dir(element, stage)) / f'{element}_{stage}'
    return [Path(str(stem) + suffix) for suffix in ('.elvlc', '.wgfa', '.scups')]


def solve_blocks(levels, radiation, collisions, block_size):
    cooling = np.empty((len(DENSITIES), len(TEMPERATURE)))
    populations = np.empty((len(DENSITIES), len(TEMPERATURE), len(levels)))
    for start in range(0, len(TEMPERATURE), block_size):
        stop = min(start + block_size, len(TEMPERATURE))
        result = full_excitation_balance(
            levels, radiation, collisions, TEMPERATURE[start:stop], DENSITIES)
        for row, density in enumerate(DENSITIES):
            cooling[row, start:stop] = result[density][0]
            populations[row, start:stop] = result[density][1]
    if not np.all(np.isfinite(cooling)) or np.any(cooling < 0):
        raise ValueError('Nonfinite or negative source cooling')
    if not np.all(np.isfinite(populations)) or np.any(populations < -64*np.finfo(float).eps):
        raise ValueError('Nonfinite or negative source populations')
    if not np.allclose(populations.sum(axis=2), 1.0, rtol=0, atol=1e-9):
        raise ValueError('Source populations are not normalized')
    return cooling, populations


def strength_audit(collisions):
    zero_count = 0
    minimum_positive = np.inf
    maximum = 0.0
    for transition in collisions:
        strength = C.upsilon(transition, TEMPERATURE)
        if strength.shape != TEMPERATURE.shape or not np.all(np.isfinite(strength)):
            raise ValueError('Invalid source effective-strength array')
        if np.any(strength < 0):
            raise ValueError('Negative source effective strength')
        zero_count += int(np.count_nonzero(strength == 0))
        positive = strength[strength > 0]
        if positive.size:
            minimum_positive = min(minimum_positive, float(positive.min()))
            maximum = max(maximum, float(positive.max()))
    return {
        'sampled_zero_effective_strengths': zero_count,
        'sampled_minimum_positive_effective_strength': (
            minimum_positive if np.isfinite(minimum_positive) else None),
        'sampled_maximum_effective_strength': maximum,
    }


def make_plots(out, ion, full, selected, valid, ratio_floored,
               density_ratio_valid, density_ratio_full, density_ratio_selected,
               full_count, selected_count):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, constrained_layout=True)
    for row, density in enumerate(DENSITIES):
        ax = axes.flat[row]
        ax.loglog(TEMPERATURE, np.where(full[row] > 0, full[row], np.nan),
                  label=f'Full: {full_count} states')
        ax.loglog(TEMPERATURE, np.where(selected[row] > 0, selected[row], np.nan),
                  '--', label=f'Selected: {selected_count} states')
        ax.set_title(r'$n_e = ' + f'{density:g}' + r'$ cm$^{-3}$')
        ax.set_xlabel(r'$T$ [K]')
        ax.set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
        ax.legend(fontsize=8)
    fig.suptitle(ion + ' bound source cooling')
    fig.savefig(out / 'cooling_profiles.pdf')
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, constrained_layout=True)
    for row, density in enumerate(DENSITIES):
        ax = axes.flat[row]
        ax.semilogx(TEMPERATURE, 100.0 * (ratio_floored[row] - 1.0),
                    color='0.75', label='Ratio with full floor')
        ax.semilogx(TEMPERATURE[valid[row]],
                    100.0 * (selected[row, valid[row]] / full[row, valid[row]] - 1.0),
                    color='C0', label='Resolved full cooling')
        ax.axhline(100 * TARGET, color='C3', linestyle=':')
        ax.axhline(-100 * TARGET, color='C3', linestyle=':')
        ax.set_title(r'$n_e = ' + f'{density:g}' + r'$ cm$^{-3}$')
        ax.set_xlabel(r'$T$ [K]')
        ax.set_ylabel(r'Cutoff difference [\%]')
        ax.legend(fontsize=7)
    fig.suptitle(ion + ' source cutoff differences')
    fig.savefig(out / 'cutoff_differences.pdf')
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, constrained_layout=True)
    for row, density in enumerate(DENSITIES):
        ax = axes.flat[row]
        mask = density_ratio_valid[row]
        ax.semilogx(TEMPERATURE[mask],
                    100.0 * (density_ratio_selected[row, mask] /
                             density_ratio_full[row, mask] - 1.0),
                    label='Selected/full density-ratio difference')
        ax.axhline(100 * TARGET, color='C3', linestyle=':')
        ax.axhline(-100 * TARGET, color='C3', linestyle=':')
        ax.set_title(r'$n_e = ' + f'{density:g}' + r'$ cm$^{-3}$')
        ax.set_xlabel(r'$T$ [K]')
        ax.set_ylabel(r'Density-ratio difference [\%]')
        ax.legend(fontsize=7)
    fig.suptitle(ion + ' density dependence relative to the lowest density')
    fig.savefig(out / 'density_ratio_differences.pdf')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ion', choices=IONS, required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--block-size', type=int, default=16)
    args = parser.parse_args()
    if args.block_size < 1 or args.block_size > len(TEMPERATURE):
        parser.error('--block-size must be between 1 and 401')
    out = args.work_dir.resolve() / args.ion
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'result.json').exists():
        raise FileExistsError(f'Refusing to overwrite completed diagnostic: {out}')
    begin = time.monotonic()
    element, stage, count = IONS[args.ion]
    levels, radiation, collisions = C.bound_excitation_atom(element, stage)
    ordered = sorted(levels)
    if len(ordered) < count:
        raise ValueError('Selected count exceeds bound source atom')
    keep = set(ordered[:count])
    selected_levels = {i: levels[i] for i in ordered if i in keep}
    selected_radiation = [r for r in radiation if r[0] in keep and r[1] in keep]
    selected_collisions = [tr for tr in collisions if tr['ll'] in keep and tr['ul'] in keep]
    for label, atom_levels, atom_radiation in (
        ('full', levels, radiation), ('selected', selected_levels, selected_radiation)
    ):
        if any(a > 0 and atom_levels[u][1] <= atom_levels[l][1]
               for l, u, _w, _g, a in atom_radiation):
            raise ValueError(f'{label} source has nonphysical radiative direction')
    source_strengths = strength_audit(collisions)
    full, full_pop = solve_blocks(levels, radiation, collisions, args.block_size)
    selected, selected_pop = solve_blocks(
        selected_levels, selected_radiation, selected_collisions, args.block_size)
    valid = full > FLOOR
    both_zero = (full == 0) & (selected == 0)
    full_zero_selected_positive = (full == 0) & (selected > 0)
    selected_zero_full_positive = (selected == 0) & (full > 0)
    ratio_floored = np.maximum(selected, FLOOR) / np.maximum(full, FLOOR)
    density_ratio_valid = (full > FLOOR) & (selected > FLOOR) & \
                          (full[0:1] > FLOOR) & (selected[0:1] > FLOOR)
    density_ratio_full = np.full(full.shape, np.nan)
    density_ratio_selected = np.full(selected.shape, np.nan)
    maxima = []
    density_ratio_maxima = []
    for row in range(len(DENSITIES)):
        error = np.abs(selected[row, valid[row]] / full[row, valid[row]] - 1.0)
        maxima.append(float(error.max()) if error.size else None)
        mask = density_ratio_valid[row]
        density_ratio_full[row, mask] = full[row, mask] / full[0, mask]
        density_ratio_selected[row, mask] = selected[row, mask] / selected[0, mask]
        difference = np.abs(density_ratio_selected[row, mask] /
                            density_ratio_full[row, mask] - 1.0)
        density_ratio_maxima.append(float(difference.max()) if difference.size else None)
    accepted = (all(value is not None and value <= TARGET for value in maxima)
                and all(value is not None and value <= TARGET for value in density_ratio_maxima)
                and bool(np.all(valid)) and bool(np.all(density_ratio_valid)))
    np.savez_compressed(
        out / 'full_domain_profiles.npz', temperature_K=TEMPERATURE,
        density_cm3=np.asarray(DENSITIES), cooling_full=full,
        cooling_selected=selected, populations_full=full_pop,
        populations_selected=selected_pop, full_above_floor=valid,
        both_zero=both_zero, full_zero_selected_positive=full_zero_selected_positive,
        selected_zero_full_positive=selected_zero_full_positive,
        ratio_to_full_floor=ratio_floored,
        density_ratio_valid=density_ratio_valid,
        density_ratio_full=density_ratio_full,
        density_ratio_selected=density_ratio_selected)
    make_plots(out, args.ion, full, selected, valid, ratio_floored,
               density_ratio_valid, density_ratio_full, density_ratio_selected,
               len(levels), count)
    artifacts = ('full_domain_profiles.npz', 'cooling_profiles.pdf',
                 'cutoff_differences.pdf', 'density_ratio_differences.pdf')
    artifacts_complete = all((out / name).is_file() and (out / name).stat().st_size > 0
                             for name in artifacts)
    if not artifacts_complete:
        raise ValueError('Incomplete full-domain diagnostic artifacts')
    hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in source_paths(element, stage)}
    elapsed = time.monotonic() - begin
    result = {
        'passed': True,
        'artifacts_complete': artifacts_complete,
        'model_accepted_at_100_to_1e6': bool(accepted),
        'relative_target': TARGET,
        'temperature_K': [100.0, 1e6],
        'temperature_samples': len(TEMPERATURE),
        'densities_cm3': list(DENSITIES),
        'full_bound_collision_connected_levels': len(levels),
        'selected_levels': count,
        'full_radiative_rows': len(radiation),
        'selected_radiative_rows': len(selected_radiation),
        'full_collision_rows': len(collisions),
        'selected_collision_rows': len(selected_collisions),
        'maximum_resolved_relative_differences': maxima,
        'maximum_resolved_density_ratio_relative_differences': density_ratio_maxima,
        'density_ratio_resolved_counts': density_ratio_valid.sum(axis=1).tolist(),
        'density_ratio_unresolved_counts': (~density_ratio_valid).sum(axis=1).tolist(),
        'full_above_floor_counts': valid.sum(axis=1).tolist(),
        'both_zero_counts': both_zero.sum(axis=1).tolist(),
        'full_zero_selected_positive_counts': full_zero_selected_positive.sum(axis=1).tolist(),
        'selected_zero_full_positive_counts': selected_zero_full_positive.sum(axis=1).tolist(),
        'full_below_or_at_floor_counts': (~valid).sum(axis=1).tolist(),
        'full_floor_erg_cm3_s': FLOOR,
        'elapsed_s': elapsed,
        'source_sha256': hashes,
        'source_strength_audit': source_strengths,
        'comparison_rule': 'Absolute cutoff error uses full cooling above the numerical floor. Density-ratio error compares Lambda(ne)/Lambda(lowest density) for both atoms only where all four cooling values exceed the floor. Full-domain acceptance requires every sample in both comparisons to be resolved.',
        'physical_scope': 'Bound collision-connected levels, spontaneous emission, electron collisions and Maxwellian detailed balance; no other processes or scattering-accuracy claim.',
    }
    atomic_write_text(out / 'result.json', json.dumps(result, indent=2) + '\n')
    table = [
        f'# {args.ion} excitation temperature-range diagnostic', '',
        '| Input | Full source and selected atom |', '|---|---|',
        f'| Temperature | 401 logarithmic samples from 100 to 1,000,000 K |',
        f'| Electron density | {", ".join(f"{x:g}" for x in DENSITIES)} cm^-3 |',
        f'| Levels | Full {len(levels)}; selected first {count} source indices |',
        f'| Radiative transitions | Full {len(radiation)}; selected {len(selected_radiation)} |',
        f'| Electron collision transitions | Full {len(collisions)}; selected {len(selected_collisions)} |',
        '| Potential and boundaries | Identical source atom; no gravitational potential or spatial domain in this statistical equilibrium calculation |',
        '| Spectral input | None; no radiation field |',
        '| Species and reactions | Same ion; spontaneous emission and electron excitation/de-excitation only |',
        '',
        f'Diagnostic completed: {result["passed"]}. Full-domain source cutoff accepted at 0.3%: {result["model_accepted_at_100_to_1e6"]}.',
        f'Numerical floor: {FLOOR:.17e} erg cm^3 s^-1. Counts below the floor are recorded in `result.json`; they provide no relative-error evidence.',
        'The full cooling, absolute cutoff, and density-ratio curves are in the three PDF files. All values and numerical masks are in `full_domain_profiles.npz`.',
        '',
    ]
    atomic_write_text(out / 'report.md', '\n'.join(table))
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
