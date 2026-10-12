"""Bound Ca VII frequency rescaling under an explicit multipole assumption.

Diagnostic only: no atom is generated or accepted. The bounds assume that
the stored total A contains only E1, M1, E2, and M2 components with fixed
line strengths. They do not bound atomic-data error or configuration mixing.
"""
import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from audit_cavii_published_data import THEORY
from mocassin_v2_atom_inventory import MOCASSIN, read_atom


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Preserve the earlier source assessment')
    source = MOCASSIN/'data/cavii.dat'
    radiation = Path.home()/'CLOUDY/c13.05/data/chianti/ca/ca_7/ca_7.wgfa'
    atom = read_atom(source, include_tables=True)
    levels = {r['index']:r for r in atom['levels']}
    old = {}
    for line in radiation.read_text().splitlines():
        if line.strip() == '-1':
            break
        fields = line.split()
        lower, upper = map(int, fields[:2])
        wavelength, gf, rate = map(float, fields[2:5])
        old[lower,upper] = (wavelength, gf, rate)
    result = []
    counts = Counter()
    for transition in atom['spontaneous_probabilities']:
        rate = transition['A_s_inverse']
        if rate <= 0:
            continue
        lower, upper = transition['lower'], transition['upper']
        wavelength, gf, previous_rate = old[lower,upper]
        assert rate == previous_rate and wavelength != 0
        jl = (levels[lower]['statistical_weight']-1)/2
        ju = (levels[upper]['statistical_weight']-1)/2
        # Published Table 2: indices 1-5 are 3s2 3p2 (even); the remaining
        # states are 3s 3p3 or 3s2 3p 3d (odd). Index is a state identity.
        same_parity = (lower <= 5) == (upper <= 5)
        candidates = ('M1','E2') if same_parity else ('E1','M2')
        allowed = [m for m in candidates if abs(jl-ju) <= int(m[1]) <= jl+ju]
        assert allowed, (lower,upper)
        gap = THEORY[upper-1]-THEORY[lower-1]
        assert gap > 0
        ratio = gap*abs(wavelength)/1e8
        factors = [ratio**(2*int(m[1])+1) for m in allowed]
        counts['unique_component' if len(allowed)==1 else 'multiple_components'] += 1
        result.append(dict(lower=lower,upper=upper,source_A_s_inverse=rate,
            source_wavelength_angstrom=wavelength,source_gf=gf,
            theoretical_gap_cm_inverse=gap,frequency_ratio=ratio,
            allowed_components_under_assumption=allowed,
            conditional_minimum_rescaling=min(factors),
            conditional_maximum_rescaling=max(factors)))
    assert len(result)==310
    collision_reorientation = []
    for r in atom['collision_tables']:
        if any(r['upsilon']) and THEORY[r['upper']-1] <= THEORY[r['lower']-1]:
            collision_reorientation.append([r['lower'],r['upper']])
    paths=[Path(__file__),source,radiation,Path(__file__).with_name('audit_cavii_published_data.py')]
    summary = dict(production_modified=False,corrected_atom_accepted=False,
        multipole_assumption='Only E1, M1, E2, M2; fixed component line strengths; no configuration change',
        formula='A_new(component) = A_source(component) * (nu_new/nu_source)^(2L+1)',
        limitation='Conditional rescaling bounds are not bounds on atomic accuracy, collision-model error, or unknown configuration provenance.',
        radiation_pairs=len(result),component_counts=dict(counts),
        minimum_conditional_rescaling=min(r['conditional_minimum_rescaling'] for r in result),
        maximum_conditional_rescaling=max(r['conditional_maximum_rescaling'] for r in result),
        theoretical_nonpositive_radiative_gaps=0,
        collision_pairs_requiring_energy_orientation=collision_reorientation,
        transitions=result,source_sha256={str(p):sha(p) for p in paths})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:summary[k] for k in ('radiation_pairs','component_counts',
        'minimum_conditional_rescaling','maximum_conditional_rescaling',
        'collision_pairs_requiring_energy_orientation')},indent=2))


if __name__ == '__main__':
    main()
