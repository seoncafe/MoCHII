"""Read-only checks of local Cl I and Cl V collision-source metadata.

This diagnostic does not adopt either source for MoCHII. In particular, the
Cl I Stout electron label conflicts with the cited paper's proton substitution.
"""

import hashlib
import json
import math
from pathlib import Path


CLOUDY = Path('/home/kiseon/CLOUDY/c25.00/data/stout/cl/cl_1')
MOCASSIN = Path('/home/kiseon/MOCASSIN/mocassin/data/clv.dat')
OUT = Path(__file__).with_name('result.json')
HC_OVER_K_CM = 1.438776877
COLL_PREF = 8.629e-6


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def substantive_lines(path):
    return [line.strip() for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]


def main():
    coll = substantive_lines(CLOUDY / 'cl_1.coll')
    temp = [float(v) for v in coll[1].split()[1:]]
    electron = [float(v) for v in coll[3].split()[4:]]
    assert coll[3].split()[:4] == ['RATE', 'ELECTRON', '1', '2']
    assert len(temp) == len(electron) == 9
    assert all(v == 4.7e-8 for v in electron)
    assert temp == sorted(temp)

    moc = substantive_lines(MOCASSIN)
    assert moc[0] == '5'
    assert moc[6].split() == ['2', '15', '0']
    tclv = [float(v) for v in moc[9:24]]
    assert len(tclv) == 15 and tclv == sorted(tclv)
    assert moc[24] == '0'
    assert moc[25] == '1 2 2.718'
    upsclv = [float(v.split()[2]) for v in moc[25:40]]
    assert len(upsclv) == 15
    assert moc[40] == '0 0 0'
    assert moc[41] == '1 2 2.98E-2'
    assert moc[-1] == '2 4 1490.8'

    delta_cl1 = 882.351 * HC_OVER_K_CM
    q_down = electron[0]
    examples = {}
    for t in (30., 1000., 3000.):
        examples[str(int(t))] = {
            'q_up_if_proton_colliders_cm3_s': q_down * (2. / 4.) * math.exp(-delta_cl1 / t),
            'pseudo_electron_upsilon_if_mislabeled': q_down * 2. * math.sqrt(t) / COLL_PREF,
        }
    result = {
        'cl_i': {
            'collision_file_sha256': sha256(CLOUDY / 'cl_1.coll'),
            'stout_label': 'RATE ELECTRON 1 2',
            'published_identity': 'proton deexcitation coefficient substituted in gamma_e column',
            'reported_q_down_cm3_s': q_down,
            'reported_temperature_grid_K': temp,
            'g_lower_upper': [4., 2.],
            'energy_gap_cm_1': 882.351,
            'energy_gap_over_k_K': delta_cl1,
            'A_down_s_1': 0.0123,
            'examples_only_not_electron_adoption': examples,
        },
        'cl_v': {
            'mocassin_file_sha256': sha256(MOCASSIN),
            'temperature_grid_K': tclv,
            'upsilon_grid': upsclv,
            'g_lower_upper': [2., 4.],
            'energy_gap_cm_1': 1490.8,
            'energy_gap_over_k_K': 1490.8 * HC_OVER_K_CM,
            'A_down_s_1': 0.0298,
            'q_down_at_10000K_cm3_s': COLL_PREF * upsclv[tclv.index(10000.)] / (4. * 100.),
        },
    }
    OUT.write_text(json.dumps(result, indent=2) + '\n')
    print(OUT)


if __name__ == '__main__':
    main()
