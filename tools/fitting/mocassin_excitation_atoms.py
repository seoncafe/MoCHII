"""Convert coherent MOCASSIN v2 excitation tables into native source atoms.

This retains the supplied finite-level physical model, not a claim of atomic
accuracy. Positive Upsilon nodes use log(T)-log(Upsilon) PCHIP rather than
v2's cubic spline. Below the first temperature, the adopted constant Upsilon
continuation is used with the actual Maxwellian temperature. No continuation
above the final node is supplied. Source files are never changed.
"""
import argparse
import hashlib
import json
from pathlib import Path
from mocassin_v2_atom_inventory import inventory, read_atom, MOCASSIN


def native_atom(record, root=MOCASSIN):
    path = root/record['relative_file']
    atom = read_atom(path, include_tables=True)
    if atom['nonpositive_photon_energy_pairs']:
        raise ValueError('Nonpositive radiative energy: identify physical levels and A before import')
    temperature = atom['temperatures_K']
    if len(temperature) < 2:
        raise ValueError('A single collision temperature does not determine a temperature curve')
    if any(t <= 0 for t in temperature) or any(b <= a for a,b in zip(temperature, temperature[1:])):
        raise ValueError('Invalid temperature grid')
    levels = {row['index']: row for row in atom['levels']}
    if min(row['energy_cm_inverse'] for row in levels.values()) != 0:
        raise ValueError('Missing zero-energy ground state')
    collisions = []
    for row in atom['collision_tables']:
        if not any(row['upsilon']):
            continue
        if any(v <= 0 for v in row['upsilon']):
            raise ValueError('Mixed zero/positive table needs a different interpolation rule')
        lo,up = sorted((row['lower'],row['upper']),
                       key=lambda i:(levels[i]['energy_cm_inverse'], i))
        collisions.append((lo,up,row['upsilon']))
    collisions.sort(key=lambda row:(row[0],row[1]))
    radiation = [row for row in atom['spontaneous_probabilities'] if row['A_s_inverse'] > 0]
    if len({(row['lower'],row['upper']) for row in radiation}) != len(radiation):
        raise ValueError('Duplicate radiative pair')
    def num(v):
        return format(v, '.17g')
    text = [f"# MOCASSIN v2 supplied excitation model: {record['element']}_{record['stage']}",
            f'# source SHA-256 {path}: {hashlib.sha256(path.read_bytes()).hexdigest()}',
            '# Atomic accuracy and original publication correspondence remain unverified.',
            '# Preserve source levels, positive spontaneous rates, and effective collision strengths.',
            '# Interpolation: log(T)-log(Upsilon) PCHIP; differs from the v2 cubic spline.',
            '# Below the first node hold Upsilon constant; Maxwellian rates use actual T.',
            '# Cold continuation is an approximation; no continuation above the final node.',
            f'# Supplied temperature nodes: {num(temperature[0])}..{num(temperature[-1])} K.',
            'NLEVEL_VERSION 2', f"NLEV {atom['nlev']}"]
    if atom['packed_radiation_reader_defect']:
        text.insert(1, '# Correct the packed Fe II reader orientation: source lower,upper,A gives downward emission upper -> lower.')
    for i,row in sorted(levels.items()):
        text.append(f"{i} {num(row['statistical_weight'])} {num(row['energy_cm_inverse'])}")
    text.append(f'NRAD {len(radiation)}')
    for row in radiation:
        text.append(f"{row['lower']} {row['upper']} {num(row['A_s_inverse'])}")
    text.extend((f'NODES {len(collisions)*len(temperature)}', f'NUPS {len(collisions)}'))
    for lo,up,values in collisions:
        text.append(f'{lo} {up} TEMPERATURE_TABLE 0 0 0 LOG_PCHIP_T {len(temperature)}')
        text.extend(f'{num(t)} {num(v)}' for t,v in zip(temperature, values))
    # This is complete for the supplied finite-level model only.
    text.extend(('# Density ratio refers only to the complete supplied finite-level model.',
                 'DENSITY_RATIO 1'))
    metadata = {k:record[k] for k in ('element','stage','source_sha256','relative_file')}
    metadata.update(nlev=atom['nlev'], nrad=len(radiation), nups=len(collisions),
                    temperature_K=temperature, packed_reader_corrected=atom['packed_radiation_reader_defect'])
    return '\n'.join(text)+'\n', metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError('Preserve earlier candidate products')
    output.mkdir(parents=True)
    result = {'scope':'Candidate source models only; not production installation or atomic accuracy verification',
              'atoms':[], 'requires_physical_review':[]}
    for row in inventory()['ions']:
        if not row['line_file_available']:
            continue
        try:
            text, metadata = native_atom(row)
        except ValueError as exc:
            result['requires_physical_review'].append({'element':row['element'],'stage':row['stage'],'reason':str(exc)})
            continue
        path = output/f"nlevel_{row['element']}_{row['stage']}.txt"
        path.write_text(text)
        metadata.update(path=str(path), native_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        result['atoms'].append(metadata)
    (output/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f"Generated {len(result['atoms'])} candidate atoms; physical review required for {len(result['requires_physical_review'])} atoms")


if __name__ == '__main__':
    main()
