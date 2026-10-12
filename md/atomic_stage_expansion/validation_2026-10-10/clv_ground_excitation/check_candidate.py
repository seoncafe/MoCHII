"""Audit indexed Cl V journal values and test a provisional two-level atom.

The journal PDF remains inaccessible. These transcribed Table 5 numbers are
diagnostic inputs, not a verified publication import or a production atom.
"""

import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RAW = Path('/home/kiseon/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/data/clv.dat')
DRIVER = ROOT / ('md/atomic_stage_expansion/validation_2026-10-10/'
                 'full_bound_native/optimized_build/check_source_excitation_atom.x')
PAPER_T = np.array([1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000,
                    6000, 7000, 8000, 9000, 10000, 11000, 12000,
                    13000, 14000, 15000, 16000, 18000, 20000,
                    24000, 28000, 32000, 36000, 40000], dtype=float)
PAPER_UPS = np.array([2.59, 2.63, 2.72, 2.83, 2.94, 3.05, 3.15,
                      3.33, 3.49, 3.62, 3.73, 3.82, 3.90, 3.98,
                      4.04, 4.10, 4.16, 4.21, 4.26, 4.35, 4.43,
                      4.57, 4.68, 4.78, 4.85, 4.92])
HC = 1.9864458571489287e-16
KB = 1.380649e-16
A21 = 2.98e-2
GAP_CM1 = 1490.5
DENSITY = np.array([1e-6, 100, 1e6, 1e10])


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def raw_nodes():
    lines = RAW.read_text().splitlines()
    initial = next(i for i, line in enumerate(lines) if line.strip() == '2 15 0')
    temperatures = np.array([float(x) for x in lines[initial+3:initial+18]])
    rows = lines[initial+19:initial+34]
    strengths = np.array([float(row.split()[2]) for row in rows])
    return temperatures, strengths


def candidate_atom():
    lines = [
        '# PROVISIONAL: indexed final-journal Table 5; original publisher PDF inaccessible.',
        '# Saraph and Storey 1999, A&AS 134, 369-375, DOI 10.1051/aas:1999437.',
        '# Tabulated collision data end at 40000 K; scattering range stated sufficient to 30000 K.',
        '# Energy and A from Cloudy c25 Stout/NIST cl_5.nrg and cl_5.tp.',
        'NLEVEL_VERSION 2', 'NLEV 2', '1 2 0',
        f'2 4 {GAP_CM1:.17g}', 'NRAD 1', f'1 2 {A21:.17g}',
        f'NODES {len(PAPER_T)}', 'NUPS 1',
        f'1 2 TEMPERATURE_TABLE 0 0 0 LOG_PCHIP_T {len(PAPER_T)}',
    ]
    lines += [f'{t:.17g} {u:.17g}' for t, u in zip(PAPER_T, PAPER_UPS)]
    output = HERE/'clv_published_table_candidate.txt'
    output.write_text('\n'.join(lines)+'\n')
    return output


def grid_file(path, mode, temperatures):
    path.write_text(f'{mode} {len(temperatures)}\n'+
                    ' '.join(f'{t:.17g}' for t in temperatures)+'\n')


def main():
    raw_t, raw_u = raw_nodes()
    matching = np.isin(raw_t, PAPER_T)
    published = dict(zip(PAPER_T, PAPER_UPS))
    comparisons = [{'temperature_K': float(t), 'mocassin_upsilon': float(u),
                    'published_printed_upsilon': float(published[t]),
                    'difference': float(u-published[t])}
                   for t, u in zip(raw_t[matching], raw_u[matching])]
    if max(abs(x['difference']) for x in comparisons) > .0050000001:
        raise ValueError('MOCASSIN nodes disagree with printed journal values')
    atom = candidate_atom()
    grid = np.geomspace(1000, 40000, 801)
    strength_t = np.unique(np.r_[grid[::4], PAPER_T])
    grid_file(HERE/'strength_grid.txt', 1, strength_t)
    grid_file(HERE/'profile_grid.txt', 2, grid)
    strength_output = HERE/'native_strengths.txt'
    profile_output = HERE/'native_profiles.txt'
    if not strength_output.exists():
        subprocess.run([str(DRIVER), str(atom), str(HERE/'strength_grid.txt'),
                        str(strength_output)], check=True, timeout=60)
    if not profile_output.exists():
        subprocess.run([str(DRIVER), str(atom), str(HERE/'profile_grid.txt'),
                        str(profile_output)], check=True, timeout=60)
    interp = PchipInterpolator(np.log(PAPER_T), np.log(PAPER_UPS), extrapolate=False)
    overlap_t = np.geomspace(2000, 40000, 801)
    mocassin_interp = CubicSpline(np.log10(raw_t), raw_u, bc_type='natural')
    mocassin_interp_difference = float(np.max(np.abs(
        mocassin_interp(np.log10(overlap_t)) /
        np.exp(interp(np.log(overlap_t))) - 1)))
    native_u = np.loadtxt(strength_output)
    expected_u = np.exp(interp(np.log(strength_t)))
    strength_err = float(np.max(np.abs(native_u[:, 1]/expected_u-1)))
    if not np.allclose(native_u[:, 0], strength_t, rtol=1e-14, atol=0):
        raise ValueError('Native Upsilon temperatures differ')
    rows = np.loadtxt(profile_output)
    if rows.shape != (4*len(grid), 5):
        raise ValueError('Incomplete two-level native profile')
    gap_erg = GAP_CM1*HC
    ups = np.exp(interp(np.log(grid)))
    q12 = 8.629e-6*ups/(2*np.sqrt(grid))*np.exp(-gap_erg/(KB*grid))
    q21 = 8.629e-6*ups/(4*np.sqrt(grid))
    errors = []
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    for j, ne in enumerate(DENSITY):
        measured = rows[j::4]
        expected_pop2 = ne*q12/(A21+ne*(q12+q21))
        expected_cooling = expected_pop2*A21*gap_erg/ne
        if not np.allclose(measured[:, 0], grid, rtol=1e-14, atol=0):
            raise ValueError('Native cooling temperatures differ')
        if not np.all(measured[:, 1] == ne):
            raise ValueError('Native density grid differs')
        err = float(np.max(np.abs(measured[:, 2]/expected_cooling-1)))
        errors.append(err)
        if err > 1e-10 or np.max(abs(measured[:, 4]-expected_pop2)) > 1e-10:
            raise ValueError(f'Native two-level solution differs at ne={ne:g}')
        ax.loglog(grid, expected_cooling, label=r'$n_e='+f'{ne:g}'+r'$ cm$^{-3}$')
        ax.loglog(grid, measured[:, 2], '--', color=ax.lines[-1].get_color())
    ax.set_xlabel(r'$T$ [K]')
    ax.set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
    ax.legend()
    fig.savefig(HERE/'candidate_cooling_profiles.pdf')
    plt.close(fig)
    result = {
        'status': 'provisional numerical candidate; publisher PDF still inaccessible',
        'published_table_K': [1000, 40000],
        'scattering_energy_supported_temperature_K': 30000,
        'mocassin_temperature_K': [float(raw_t[0]), float(raw_t[-1])],
        'mocassin_original_nodes': len(raw_t),
        'published_printed_nodes': len(PAPER_T),
        'overlapping_raw_nodes': len(comparisons),
        'comparisons': comparisons,
        'maximum_absolute_raw_vs_printed_difference': max(abs(x['difference']) for x in comparisons),
        'maximum_relative_mocassin_vs_table_interpolation_difference_2000_to_40000_K': mocassin_interp_difference,
        'maximum_relative_native_strength_error': strength_err,
        'maximum_relative_native_cooling_error': max(errors),
        'native_grid_points': len(grid),
        'density_cm3': DENSITY.tolist(),
        'candidate_atom': str(atom),
        'source_hashes': {str(p): sha256(p) for p in [RAW, DRIVER, Path(__file__)]},
        'candidate_hash': sha256(atom),
        'publisher_pdf': 'https://aas.aanda.org/articles/aas/pdf/1999/02/ds1586.pdf',
        'published_text_index': 'https://www.researchgate.net/publication/253614097_Atomic_data_from_the_IRON_Project_XXX_Collision_data_for_the_Po_1over_2_-_Po_3over_2_fine-structure_lines_of_P_iii_S_iv_and_Cl_v',
    }
    (HERE/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['status','published_printed_nodes',
          'overlapping_raw_nodes','maximum_absolute_raw_vs_printed_difference',
          'maximum_relative_mocassin_vs_table_interpolation_difference_2000_to_40000_K',
          'maximum_relative_native_strength_error','maximum_relative_native_cooling_error']},
          indent=2))


if __name__ == '__main__':
    main()
