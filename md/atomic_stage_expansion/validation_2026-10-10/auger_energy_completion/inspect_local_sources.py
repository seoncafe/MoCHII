"""Record read-only Auger source coverage and threshold contradictions."""

import hashlib
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
OUTPUT = Path(__file__).with_name('source_inventory.json')
AUDIT = ROOT / 'md/atomic_stage_expansion/validation_2026-10-10/auger_audit_cds_checked/auger_channels.json'
CLOUDY = Path('/home/kiseon/CLOUDY/c25.00')
KM = Path('/home/kiseon/RT_Codes/ExoAtmosphere/references/kaastra_mewe_1993')
JAC = Path('/home/kiseon/RT_Codes/JAC.jl-master')
SOURCES = [
    AUDIT,
    CLOUDY / 'data/auger_spectrum.dat',
    CLOUDY / 'data/mewe_nelectron.dat',
    CLOUDY / 'data/mewe_fluor.dat',
    CLOUDY / 'source/grains.cpp',
    CLOUDY / 'source/ion_solver.cpp',
    CLOUDY / 'source/cont_gammas.cpp',
    KM / 'ReadMe', KM / 'table2', KM / 'table3',
    JAC / 'src/module-Cascade-inc-photoionization.jl',
    JAC / 'src/module-Cascade-inc-stepwise-decay.jl',
    JAC / 'src/module-Cascade-inc-simulations.jl',
    Path('/nfs/mocafe/kiseon/MoCHII/references/Kaastra_Mewe_1993_AAS_97_443.pdf'),
]


def grain_auger_coverage(path):
    numeric = []
    for line in path.read_text().splitlines():
        stripped = line.split('#', 1)[0].strip()
        if stripped:
            numeric.append(stripped)
    assert numeric[0] == '20060126'
    at = 1
    elements = []
    while numeric[at] != '-1':
        index = int(numeric[at]); at += 1
        n_shells = int(numeric[at]); at += 1
        shells = []
        for _ in range(n_shells):
            shell = int(numeric[at]); at += 1
            edge = float(numeric[at]); at += 1
            n_energies = int(numeric[at]); at += 1
            entries = []
            for _ in range(n_energies):
                energy, average_count = map(float, numeric[at].split())
                entries.append([energy, average_count]); at += 1
            shells.append({'shell_index_0based': shell, 'edge_eV': edge,
                           'mean_spectrum_points': len(entries)})
        elements.append({'element_index_0based': index, 'atomic_number': index + 1,
                         'shells': shells})
    assert at == len(numeric) - 1
    return elements


def main():
    channels = json.loads(AUDIT.read_text())['channels']
    two_electron = [c for c in channels if c.get('ejected_electrons') == 2 and c.get('inner')]
    assert len(two_electron) == 9
    rows = []
    for c in two_electron:
        assert c['cds_table2_mapping_complete']
        rows.append({
            'element': c['element'], 'stage': c['source_stage'], 'shell': c['shell'],
            'ph1_edge_eV': c['shell_edge_eV'],
            'sequential_final_threshold_eV': c['minimum_sequential_IP_eV'],
            'diagnostic_reachability_onset_eV': c['channel_onset_eV'],
            'km_initial_vacancy_edges_eV': c['cds_table2_edges_eV'],
            'km_charge_two_probabilities': c['cds_table2_probabilities'],
            'km_mean_Auger_energies_eV': c['cds_mean_auger_eV'],
        })
    result = {
        'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCES},
        'reachable_charge_two_channels_diagnostic_only': rows,
        'grain_auger_spectrum_coverage': grain_auger_coverage(CLOUDY / 'data/auger_spectrum.dat'),
        'grain_spectrum_invoked_in': str(CLOUDY / 'source/grains.cpp'),
        'cloudy_gas_charge_yields_invoked_in': str(CLOUDY / 'source/ion_solver.cpp'),
        'joint_charge_electron_photon_energies_found': False,
        'final_article_pages_visually_checked': [446, 447],
        'julia_executable_available': shutil.which('julia') is not None,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + '\n')
    print(OUTPUT)


if __name__ == '__main__':
    main()
