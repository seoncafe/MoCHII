"""Reader for Badnell adf48 (RR) / adf09 (DR) final-state-resolved partial
recombination files (IC), initial parent PRTI=1 only."""
import gzip, re
import numpy as np

NUM = re.compile(r'[-+]?\d+\.\d+(?:[Ee][-+]?\d+|[-+]\d+)')
LEV = re.compile(r'^\s*(\d+)\s+(\d+)\s+(.+?)\s+\(\s*(\d+)\)\s*(\d+)\(\s*([\d.]+)\)\s+(-?[\d.]+)')
PAR2 = re.compile(r'^\s*(\d+)\s+([0-9A-Z ]+?)\s+\(\s*([\d.]+)\)\s+(-?[\d.]+)')
PAR = re.compile(r'^\s*(\d+)\s+(.+?)\s+\(\s*(\d+)\)\s*(\d+)\(\s*([\d.]+)\)\s+(-?[\d.]+)')

def fnum(s):
    s = s.strip()
    if re.match(r'.*\d[-+]\d+$', s) and 'E' not in s.upper():
        s = re.sub(r'(\d)([-+]\d+)$', r'\1E\2', s)
    return float(s)

def read(path):
    lines = gzip.open(path, 'rt', errors='replace').read().splitlines()
    d = dict(parents={}, levels={}, nrep={}, nlrep={}, prtf={}, indx={})
    m = re.search(r'NUCCHG=\s*(\d+)', lines[0]); d['Z'] = int(m.group(1))
    i = 0; sect = None
    while i < len(lines):
        L = lines[i]
        if 'PARENT CONFG INDEXING' in L:
            d['parents'] = {}; sect = 'par'
        elif 'AVERAGED CONFG INDEXING' in L:
            sect = 'cfglev'
        elif 'PARENT LEVEL INDEXING' in L:
            d['BWNP'] = float(re.search(r'BWNP=\s*([\d.]+)', L).group(1)); sect = 'par'
        elif 'RESOLVED LEVEL INDEXING' in L:
            d['BWNR'] = float(re.search(r'BWNR=\s*([\d.]+)', L).group(1)); sect = 'lev'
        elif 'NL-SHELL INDEXING' in L: sect = 'nl'
        elif 'N-SHELL INDEXING' in L: sect = 'n'
        elif re.match(r'^\s*PRTI=', L): sect = 'prti'; break
        elif sect == 'par':
            m = PAR.match(L)
            m2 = PAR2.match(L)
            if m2 and not m:
                d['parents'][int(m2.group(1))] = dict(code=m2.group(2).strip(), S=0, L=0, J=float(m2.group(3)), E=float(m2.group(4)))
            if m and not L.strip().startswith('INDP'):
                d['parents'][int(m.group(1))] = dict(code=m.group(2).strip(), S=int(m.group(3)), L=int(m.group(4)),
                                                     J=float(m.group(5)), E=float(m.group(6)))
        elif sect == 'lev':
            m = LEV.match(L)
            if m:
                d['levels'][int(m.group(1))] = dict(indp=int(m.group(2)), code=m.group(3).strip(), S=int(m.group(4)),
                                                    L=int(m.group(5)), J=float(m.group(6)), E=float(m.group(7)))
        elif sect == 'nl':
            t = L.split()
            if len(t) >= 3 and t[0].isdigit() and t[1].isdigit() and t[2].isdigit():
                d['nlrep'][int(t[0])] = (int(t[1]), int(t[2]))
        elif sect == 'n':
            t = L.split()
            if len(t) >= 2 and t[0].isdigit() and t[1].isdigit():
                d['nrep'][int(t[0])] = int(t[1])
        i += 1
    # PRTI blocks
    prti = None; sub = None; prtf = None; te = []
    while i < len(lines):
        L = lines[i]
        m = re.match(r'^\s*PRTI=\s*(\d+)', L)
        if m:
            prti = int(m.group(1)); sub = 'indx'; i += 1; continue
        if prti is None or prti > 1:
            if prti is not None and prti > 1: break
            i += 1; continue
        if 'INDX TE=' in L:
            te = [fnum(x) for x in L.split('TE=')[1].split()] + [fnum(x) for x in lines[i+1].split()[2:]]
            d['TE'] = np.array(te); i += 2; continue
        m = re.search(r'PRTF=\s*(\d+)', L)
        if m:
            prtf = int(m.group(1)); d['prtf'].setdefault(prtf, {'nl': {}, 'n': {}}); i += 1; continue
        s = L.strip()
        if s.startswith('ILREP'): sub = 'ilrep'; i += 1; continue
        if s.startswith('INREP'): sub = 'inrep'; i += 1; continue
        if s.startswith('C---') or (s.startswith('----') and 'TE' not in s and len(s.split()) > 3):
            break
        t = L.split()
        if len(t) == 11 and t[0].isdigit() and i + 1 < len(lines):
            v = [fnum(x) for x in t[1:]] + [fnum(x) for x in lines[i+1].split()]
            if len(v) == 19:
                k = int(t[0])
                if sub == 'indx': d['indx'][k] = np.array(v)
                elif sub == 'inrep': d['prtf'][prtf]['n'][d['nrep'][k]] = np.array(v)
                elif sub == 'ilrep': d['prtf'][prtf]['nl'][d['nlrep'][k]] = np.array(v)
            i += 2; continue
        i += 1
    return d
