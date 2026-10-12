"""Consistent physical excitation energies without changing transition identities.

O V: observed levels 1--40, then the original 2014 ICFT target energies
mapped from Stout to retained CHIANTI levels by weight, energy and term.
Cl VI: observed levels 1--44, then original CHIANTI ICFT target energies.
These splices preserve low-level observations and avoid mixed high-level
order reversals. They do not shift scattering resonances or certify accuracy
at low temperature. Higher-level convergence must be assessed separately.
"""
from hashlib import sha256
from pathlib import Path
import re

from niv_level_energies import niv_energies

STOUT_OV=Path.home()/'CLOUDY/c25.00/data/stout/o/o_5/o_5.nrg'
STOUT_OV_SHA256='78689cda442d424a95c5e7e1967a32d88ee02cc839b1214a380a2321f7b6eaee'
ENERGY_HASHES={('o',5):'041886ff6e3e44a5d920d81250cc381f601d068019c08f6bd551cdbe2cb69247',
               ('cl',6):'7c219e83004418daeacbed19deaab7cfebf16523a6375c5e940196bf8d608348'}
MGIII_LEVEL_SHA256='07e61cf56c782806986ca9a19c2197f2cfa85adaf568be4b0ab711850084851e'


def orbital_configuration(label):
    """Normalize orbital occupations, excluding closed 1s and parent terms."""
    label=re.sub(r'\([^)]*\)','',label)
    return tuple((int(n),l.upper(),int(occ or '1'))
                 for n,l,occ in re.findall(r'(\d+)([spdfgSPDFG])(\d*)',label)
                 if not (n=='1' and l.upper()=='S'))


def icft_energies(element,stage,levels,radiation,collisions):
    """Checked observed/ICFT energy splice for retained O V or Cl VI levels."""
    key=(element,stage)
    path=Path.home()/f'RT_Codes/CHIANTI/dbase/{element}/{element}_{stage}/{element}_{stage}.elvlc'
    raw=path.read_bytes()
    if sha256(raw).hexdigest()!=ENERGY_HASHES[key]:
        raise ValueError(f'ICFT level source changed: {path}')
    count,nrad,ncoll,observed_count=(166,2995,13695,40) if element=='o' else (171,3760,14535,44)
    if set(levels)!=set(range(1,count+1)) or len(radiation)!=nrad or len(collisions)!=ncoll:
        raise ValueError('ICFT retained level/transition set changed')
    if element=='cl':
        result={i:(g,obs if i<=observed_count else theory) for i,(g,obs,theory) in levels.items()}
    else:
        source=STOUT_OV.read_bytes()
        if sha256(source).hexdigest()!=STOUT_OV_SHA256:
            raise ValueError('O V high-precision Stout source changed')
        records={}
        for line in source.decode('ascii').splitlines()[1:]:
            if line.startswith('*'):break
            fields=line.split(maxsplit=3)
            if len(fields)>=4 and fields[0].isdigit():
                records[int(fields[0])]=(float(fields[2]),float(fields[1]),fields[3].strip('"'))
        if set(records)!=set(range(1,239)):raise ValueError('Incomplete O V source')
        result={};used=set()
        for line in raw.decode('ascii').splitlines():
            if line.strip()=='-1':break
            i=int(line[:7]);g,obs,theory=levels[i]
            candidates=([i] if i<=observed_count else
                        [j for j,(weight,energy,label) in records.items()
                         if weight==g and abs(energy-theory)<=.51])
            if len(candidates)!=1:raise ValueError(f'O V ambiguous theoretical level {i}')
            j=candidates[0];weight,energy,label=records[j]
            config,term=label.rsplit(' ',1)
            identity=str(int(line[42:47]))+line[47:52].strip()+str(int(float(line[52:57])))
            if (j in used or weight!=g or orbital_configuration(line[7:37])!=orbital_configuration(config)
                or re.sub(r'[^A-Z0-9]','',term.upper())!=identity
                or (i<=observed_count and abs(energy-obs)>.001)):
                raise ValueError(f'O V source identity mismatch {i}->{j}')
            used.add(j);result[i]=(weight,energy)
    if result[1][1]!=0 or any(e<=0 for i,(g,e) in result.items() if i!=1):
        raise ValueError('Invalid physical ICFT energies')
    for lower,upper,_,_,rate in radiation:
        if rate>0 and result[upper][1]<=result[lower][1]:
            raise ValueError(f'ICFT radiative energy direction {lower}->{upper}')
    for tr in collisions:
        if result[tr['ul']][1]<=result[tr['ll']][1]:
            raise ValueError(f'ICFT collision energy direction {tr["ll"]}->{tr["ul"]}')
    return result


def mgiii_excitation_energies(levels,radiation,collisions):
    """Reject a misassigned Mg III J=0 observation at source level 167.

    CHIANTI assigns 628148.375 cm^-1 to 2s2 2p5 7p 1S, J=0, although
    the NIST ASD Mg III level table associates 628149 cm^-1 with a 7d J=1
    state. Use this same J=0 source state's theoretical 634430 cm^-1 energy.
    The theoretical value is an identified-state estimate, not a measured
    energy or a claim of atomic accuracy. All other selected energies and
    original collision nodes are preserved.
    Source comparison: https://physics.nist.gov/PhysRefData/ASD/levels_form.html
    See md/atomic_stage_expansion/mgiii_radiative_consistency.md.
    """
    path=Path.home()/'RT_Codes/CHIANTI/dbase/mg/mg_3/mg_3.elvlc'
    source=path.read_bytes()
    if sha256(source).hexdigest()!=MGIII_LEVEL_SHA256:
        raise ValueError('Mg III level source changed')
    if (set(levels)!=set(range(1,210)) or len(radiation)!=4753
            or len(collisions)!=21736 or levels[167]!=(1.0,628148.375,634430.0)):
        raise ValueError('Mg III full source identities changed')
    line=next((row for row in source.decode('ascii').splitlines()
               if row[:7].strip()=='167'),None)
    if (line is None or line[7:37].strip()!='2s2 2p5 7p'
            or int(line[42:47])!=1 or line[47:52].strip()!='S'
            or float(line[52:57])!=0.0):
        raise ValueError('Mg III level 167 term or angular momentum changed')
    result={i:(g,obs if (obs>0 or i==1) else theory)
            for i,(g,obs,theory) in levels.items()}
    result[167]=(levels[167][0],levels[167][2])
    return result


def selected_energies(element,stage,levels,radiation,collisions):
    """Physical energies shared by excitation, detailed balance and emission."""
    element=element.lower()
    if (element,stage)==('mg',3):
        return mgiii_excitation_energies(levels,radiation,collisions)
    if element=='n' and stage==4:return niv_energies(levels,radiation,collisions)
    if element=='s' and stage==6:
        from valence_excitation_energies import svi_fine_structure_energies
        return svi_fine_structure_energies(levels)
    if (element,stage) in [('o',4),('ne',5),('s',7),('ar',5),('mg',4)]:
        from valence_excitation_energies import valence_energies
        return valence_energies(element,stage,levels,radiation,collisions)
    if (element,stage) in ENERGY_HASHES:return icft_energies(element,stage,levels,radiation,collisions)
    return {i:(g,obs if (obs>0 or i==1) else theory) for i,(g,obs,theory) in levels.items()}
