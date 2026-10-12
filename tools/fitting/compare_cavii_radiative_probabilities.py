"""Ca VII radiative probabilities: CHIANTI 11.0.2 ca_7 against Landi & Bhatia (2003).

Landi & Bhatia (2003, ApJ 589, 1075) print A values of a 5-configuration and
a 25-configuration calculation (Tables 3-5); the 55 printed pairs and both
columns were transcribed from the final article in
md/atomic_stage_expansion/validation_2026-10-11/cavii_published_review_checked/
published_source_review.json. This program compares every CHIANTI ca_7.wgfa
row present there with both columns, and checks for every row that the photon
energy on CHIANTI's selected energies (observed where measured, otherwise
theoretical) is positive and that the frequency at which A and gf describe the
same line strength, A = 6.6703e15 gf / (g_u lambda^2), agrees with the listed
wavelength. Section 2.1 of the article states that A and gf were corrected to
observed energies where both levels are observed.

The consequence of the column choice is measured with the installed atom
(data/atomic/excitation/chianti/nlevel_ca_7.txt): level populations at
T = 1e4, 2e4 K and n_e = 1e2, 1e4, 1e6 cm^-3 with the CHIANTI A values, and
with every Table 3 pair replaced by its 25-configuration value; the strongest
optical/IR lines and the total cooling are compared.

Run: python3 compare_cavii_radiative_probabilities.py --output DIR
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

import chianti_cooling as C
from chianti_cooling import ion_dir, read_elvlc, read_wgfa
from validate_chianti_excitation_model import file_upsilon, parse_atom

ROOT = Path(__file__).resolve().parents[2]
REVIEW = ROOT / ('md/atomic_stage_expansion/validation_2026-10-11/'
                 'cavii_published_review_checked/published_source_review.json')
EINSTEIN = 6.6703e15
ATOM = ROOT / 'data/atomic/excitation/chianti/nlevel_ca_7.txt'


def line_emissivities(atom, einstein, temperature, density):
    """Populations and line power per ion [erg s^-1] of a file atom."""
    g, energy = atom['g'], atom['E']
    n = len(g)
    rates = np.zeros((n, n))  # rates[j, i]: transitions i -> j per second
    for (lower, upper), a in einstein.items():
        rates[lower-1, upper-1] += a
    for row in atom['coll']:
        lower, upper = row['l']-1, row['u']-1
        if energy[upper] < energy[lower]:
            lower, upper = upper, lower
        common = C.COLL_PREF*file_upsilon(atom, row, temperature)[0]/np.sqrt(temperature)
        gap = (energy[upper]-energy[lower])*C.HC_ERG_CM
        rates[upper, lower] += density*common/g[lower]*np.exp(-gap/(C.K_B_ERG*temperature))
        rates[lower, upper] += density*common/g[upper]
    matrix = rates - np.diag(rates.sum(axis=0))
    matrix[0, :] = 1.0
    rhs = np.zeros(n); rhs[0] = 1.0
    population = np.linalg.solve(matrix, rhs)
    power = {(l, u): population[u-1]*a*(energy[u-1]-energy[l-1])*C.HC_ERG_CM
             for (l, u), a in einstein.items()}
    return population, power


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    output = parser.parse_args(argv).output
    stem = Path(ion_dir('ca', 7)) / 'ca_7'
    levels = read_elvlc(f'{stem}.elvlc')
    rows = read_wgfa(f'{stem}.wgfa')
    selected = {i: (g, obs if (obs > 0 or i == 1) else th) for i, (g, obs, th) in levels.items()}
    review = json.loads(REVIEW.read_text())
    published = {}
    for table, record in review['tables'].items():
        for row in record['rows']:
            published[(row['lower'], row['upper'])] = (
                table, row['published_five_configuration_A_s_inverse'],
                row['published_twenty_five_configuration_A_s_inverse'])
    radiation, comparison = [], []
    for lower, upper, wavelength, gf, a in rows:
        gap = selected[upper][1] - selected[lower][1]
        both_observed = all(levels[i][1] > 0 or i == 1 for i in (lower, upper))
        record = dict(lower=lower, upper=upper, A=a, gf=gf, wavelength_A=wavelength,
                      selected_gap_cm=gap, both_observed=both_observed)
        if gf > 0 and a > 0:
            source_frequency = 1e8 * (a * levels[upper][0] / (EINSTEIN * gf)) ** .5
            record['A_gf_frequency_cm'] = source_frequency
            record['A_gf_frequency_over_listed'] = source_frequency * abs(wavelength) / 1e8
            record['A_gf_frequency_over_selected_gap'] = source_frequency / gap if gap > 0 else None
        radiation.append(record)
        if (lower, upper) in published:
            table, five, twenty_five = published[lower, upper]
            comparison.append(dict(lower=lower, upper=upper, table=table, chianti_A=a,
                                   five_configuration_A=five, twenty_five_configuration_A=twenty_five,
                                   chianti_over_five=a / five, chianti_over_twenty_five=a / twenty_five))
    nonpositive = [r for r in radiation if r['selected_gap_cm'] <= 0]
    def matches(r, key):
        return abs(r[key] - 1) <= 2e-3
    tally = {}
    for r in comparison:
        column = ('25' if matches(r, 'chianti_over_twenty_five') else
                  '5' if matches(r, 'chianti_over_five') else 'neither')
        tally.setdefault(r['table'], {}).setdefault(column, 0)
        tally[r['table']][column] += 1
    e1_ratio = [r['A_gf_frequency_over_listed'] for r in radiation if 'A_gf_frequency_over_listed' in r]
    result = dict(
        source_sha256={str(Path(f'{stem}.{e}')): sha256(Path(f'{stem}.{e}').read_bytes()).hexdigest()
                       for e in ('elvlc', 'wgfa', 'scups')},
        published_review_sha256=sha256(REVIEW.read_bytes()).hexdigest(),
        levels=len(levels), radiative_rows=len(rows),
        nonpositive_selected_photon_energies=len(nonpositive),
        rows_compared_with_published=len(comparison),
        matching_column_by_table=tally,
        rows_with_gf=len(e1_ratio),
        A_gf_frequency_over_listed_range=[min(e1_ratio), max(e1_ratio)] if e1_ratio else None,
        comparison=comparison, radiation=radiation)
    atom = parse_atom(ATOM)
    chianti_a = {(l, u): a for l, u, a in atom['rad']}
    twenty_five = dict(chianti_a)
    for (lower, upper), (table, five, twenty) in published.items():
        if table == '3' and (lower, upper) in twenty_five:
            twenty_five[lower, upper] = twenty
    sensitivity = []
    for temperature in (1e4, 2e4):
        for density in (1e2, 1e4, 1e6):
            _, base = line_emissivities(atom, chianti_a, temperature, density)
            _, alternative = line_emissivities(atom, twenty_five, temperature, density)
            total = sum(base.values())
            strongest = sorted(base, key=base.get, reverse=True)[:8]
            sensitivity.append(dict(
                temperature_K=temperature, density_cm3=density,
                cooling_per_ion_25conf_over_chianti=sum(alternative.values())/total,
                lines=[dict(lower=l, upper=u,
                            wavelength_A=1e8/(atom['E'][u-1]-atom['E'][l-1]),
                            share_of_cooling=base[l, u]/total,
                            power_25conf_over_chianti=alternative[l, u]/base[l, u])
                       for l, u in strongest]))
    result['table3_25conf_substitution_sensitivity'] = sensitivity
    output.mkdir(parents=True, exist_ok=True)
    (output / 'cavii_radiative_probabilities.json').write_text(json.dumps(result, indent=1) + '\n')
    print(f"rows {len(rows)}, nonpositive gaps {len(nonpositive)}, compared {len(comparison)}")
    print('matching column by table:', tally)
    if e1_ratio:
        print(f'A/gf frequency over listed wavelength: {min(e1_ratio):.5f} .. {max(e1_ratio):.5f}')
    for record in sensitivity:
        print(f"T={record['temperature_K']:.0e} ne={record['density_cm3']:.0e} "
              f"cooling 25conf/CHIANTI={record['cooling_per_ion_25conf_over_chianti']:.4f}")
        for line in record['lines'][:5]:
            print(f"   {line['lower']}-{line['upper']} {line['wavelength_A']:.1f} A "
                  f"share {line['share_of_cooling']:.3f} ratio {line['power_25conf_over_chianti']:.4f}")
    for r in comparison:
        print(f"  T{r['table']} {r['lower']:2d}-{r['upper']:2d} A={r['chianti_A']:.4e}"
              f"  /5conf={r['chianti_over_five']:.4f}  /25conf={r['chianti_over_twenty_five']:.4f}")


if __name__ == '__main__':
    main()
