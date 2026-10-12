"""Observed valence complexes and original scattering-target higher energies.

Preserve CHIANTI identities and A directions. Observed/theory splices do not
shift scattering resonances and have unquantified low-temperature accuracy.
Unresolved observed terms use theoretical splitting around the statistical-
weight term centroid, never an arbitrary positive energy floor.
"""
from hashlib import sha256
from pathlib import Path
import numpy as np

SOURCE_HASHES = {('o', '4'): {'elvlc': '73281e91f605977186a25798bed25c31f74e9a76a727e2954178538c7b899eb3', 'wgfa': '1ab5bb9f28302adff302058279c13292fe363d028bf0fbe1182da01d44a7013d', 'scups': '5bd169d70125a3a52553f0d2b60346a5796eb192bc3b0ae3b67eabc0747406d4'}, ('ne', '5'): {'elvlc': '91a020c35411e0508d82b9e835f5021e02305877bf7b8ab495996aa40668b1db', 'wgfa': '131de666bc09f714025791449a676b2c4019011803d9cfcc0a833bc61a8f4822', 'scups': '387bdb0d55dec9cce35edf07a9d99faa31b3e4814f24495ed5b4398b14415584'}, ('s', '7'): {'elvlc': 'd80fc3a98e1f5e39d6deb22b477ffb353fee5d6978554e86849f0c4fc8ef8b26', 'wgfa': 'de9e5a4d2475c2f657f09f9b432b857a0a8ad20aaa36bde58cb6fa2f9e6e5d20', 'scups': 'ea3c71f8e8fd9dab1391f8487991dbdde098a7271989c15d64c0d75015ecb25a', 'stout_nrg': 'e2b34144eed47b7498622424b869cd25950f2b5e061d3ed824843320515d3707'}, ('ar', '5'): {'elvlc': 'b9fa8cf6970696cf69abba11274fbc452f61e3eb82f99226c05e1ca139e576b0', 'wgfa': '25cab8b9436852066de244aac909ceccfde0a8a7c230c18a6b37b9a8c5821323', 'scups': '84554aacde47f1822bbeea3f793a77894f4a265deadb531044a59ab7f02c71be'}, ('mg', '4'): {'elvlc': '60b2fdc1d841f283efd7b02a2451d46960018b9ab745c06ca8f1ba50647224df', 'wgfa': 'c81a65bb80902e52163f91277de264569b33f24d6bbe87bb4309baa3aaf93fb1', 'scups': '4605c20dccc3e5cea9a4a4fa5b824ceb0f45db8725f95594305b91aceed6adae', 'stout_nrg': 'd72b9f0df6142289f0ef1c82a53a0072d3301d1acd67c7ef162119231b8c4066'}}
VALENCE_N = {'o':2,'ne':2,'s':3,'ar':3,'mg':2}


def valence_energies(element,stage,raw,radiation,collisions):
    """Select identified observed low complexes, coherent theory above."""
    element=element.lower();key=(element,str(stage))
    stem=Path.home()/f'RT_Codes/CHIANTI/dbase/{element}/{element}_{stage}/{element}_{stage}'
    hashes=SOURCE_HASHES[key]
    for suffix in ['elvlc','wgfa','scups']:
        if sha256(Path(str(stem)+'.'+suffix).read_bytes()).hexdigest()!=hashes[suffix]:
            raise ValueError('Valence excitation source changed')
    # Read original labels to select whole principal-shell configurations.
    import re
    result={}
    for line in Path(str(stem)+'.elvlc').read_text().splitlines():
        if line.strip()=='-1':break
        i=int(line[:7]);g,obs,theory=raw[i]
        configuration=line[7:37];shells=[int(n) for n in re.findall(r'(\d+)[spdfg]',configuration)]
        observed=(shells and max(shells)<=VALENCE_N[element]) or (element=='mg' and '4f' in configuration)
        result[i]=(g,obs if observed and (obs>0 or i==1) else theory)
    if element=='s':
        path=Path.home()/'CLOUDY/c25.00/data/stout/s/s_7/s_7.nrg'
        if sha256(path.read_bytes()).hexdigest()!=hashes['stout_nrg']:raise ValueError('S VII observed source changed')
        for line in path.read_text().splitlines()[1:]:
            if line.startswith('*'):break
            fields=line.split(maxsplit=3);i=int(fields[0])
            if i in (23,24):
                if raw[i][0]!=float(fields[2]):raise ValueError('S VII state identity changed')
                # Source identity: 2s2.2p5.3d 3D*1 and 1D*2 respectively.
                label='3D*1' if i==23 else '1D*2'
                if '2s2.2p5.3d '+label not in fields[3]:raise ValueError('S VII term identity mismatch')
                result[i]=(raw[i][0],float(fields[1]))
    if element=='mg':
        ids=[133,134];centroid=sum(raw[i][0]*raw[i][2] for i in ids)/sum(raw[i][0] for i in ids)
        if raw[133][1]!=raw[134][1]:raise ValueError('Mg IV unresolved term changed')
        for i in ids:result[i]=(raw[i][0],raw[i][1]+raw[i][2]-centroid)
    # Only bound radiation belongs to the ordinary excitation model.
    from ionization_potentials import ELEMENTS,read_ionization_potentials,HC_EV_CM
    z=dict(ELEMENTS)[element];edge=read_ionization_potentials()[z,z-stage+1]/HC_EV_CM
    bound={i for i,(g,e) in result.items() if e<edge}
    for l,u,w,g,a in radiation:
        if l in bound and u in bound and a>0 and result[u][1]<=result[l][1]:
            raise ValueError('Valence energy vector reverses source radiation')
    for tr in collisions:
        l,u=tr['ll'],tr['ul']
        if l in bound and u in bound and result[u][1]==result[l][1]:
            raise ValueError('Valence source has an unresolved equal collision gap')
    return result


def svi_fine_structure_energies(raw):
    """Observed unresolved centers with original theoretical term splittings."""
    result={i:(g,obs if obs>0 or i==1 else theory) for i,(g,obs,theory) in raw.items()}
    for ids in [(11,12),(18,19),(20,21),(27,28),(29,30),(31,32)]:
        a,b=ids
        if raw[a][1]!=raw[b][1]:raise ValueError('S VI unresolved observation changed')
        centroid=sum(raw[i][0]*raw[i][2] for i in ids)/sum(raw[i][0] for i in ids)
        for i in ids:result[i]=(raw[i][0],raw[i][1]+raw[i][2]-centroid)
    return result
