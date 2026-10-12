#!/usr/bin/env python3
"""Ground-term statistical weights of every metal stage MoCHII carries.

WHAT THIS PRODUCES
------------------
data/atomic/ground_term_weights.txt: for each element (C N O Ne S Ar Mg Fe Si
Cl Ca) and stage I-VII, the statistical weight g = (2S+1)(2L+1) of the ground
term and the energy of the stage's first excited term [eV].

WHO READS IT
------------
src/milne_recomb_spectrum_mod.f90 (metal ground continua of the diffuse
field).  The Milne relation gives the ground-state radiative recombination
X^(i+1) + e -> X^i(ground) + h nu from the ground photoionization cross
section with the factor g(X^i ground) / (2 g(X^(i+1) ground)), so both
weights are needed for every transition.  The excited-term energy is not
used by the physics; it bounds the part of the Milne integral where the
outer-shell cross section also contains photoionization into excited terms
of the ion (see the validity note in milne_recomb_spectrum_mod), and setup
prints that bound.

SOURCES
-------
Ground terms: NIST Atomic Spectra Database, ground states of the neutral
atoms and ions (the LS term of the ground level; Fe I a 5D, Fe II a 6D,
Fe III 5D, Fe IV 6S, Fe V 5D, Fe VI 4F, Fe VII 3F; every other stage follows
its isoelectronic sequence).  Written below explicitly and checked against
CHIANTI 11.0.2 .elvlc: the levels that share the configuration, 2S+1 and L
of level 1 are summed as 2J+1.  The script stops on any disagreement.
First excited term: CHIANTI 11.0.2 .elvlc, the lowest level whose
configuration or LS term differs from level 1 (observed energy, else the
theoretical one).  Where CHIANTI has no .elvlc for the stage (Cl I, Cl V,
Cl VII, Ar I, Ar II, Ca I, Ca III, Ca IV, Fe I) the energy is written as -1
(not available); C VII is bare and has no excited bound term below the
continuum that matters here (also -1).

Run:  python3 make_ground_term_weights.py [CHIANTI dbase dir]
      (default ~/RT_Codes/CHIANTI/dbase; writes data/atomic/ground_term_weights.txt)
"""
import os
import sys
import datetime

CM_TO_EV = 1.239841984e-4          # h c [eV cm]

# (element, Z) and ground terms of stages I..VII as (2S+1, L letter)
GROUND = {
    ('c', 6):  ['3P', '2P', '1S', '2S', '1S', '2S', '1S'],
    ('n', 7):  ['4S', '3P', '2P', '1S', '2S', '1S', '2S'],
    ('o', 8):  ['3P', '4S', '3P', '2P', '1S', '2S', '1S'],
    ('ne', 10): ['1S', '2P', '3P', '4S', '3P', '2P', '1S'],
    ('mg', 12): ['1S', '2S', '1S', '2P', '3P', '4S', '3P'],
    ('si', 14): ['3P', '2P', '1S', '2S', '1S', '2P', '3P'],
    ('s', 16): ['3P', '4S', '3P', '2P', '1S', '2S', '1S'],
    ('cl', 17): ['2P', '3P', '4S', '3P', '2P', '1S', '2S'],
    ('ar', 18): ['1S', '2P', '3P', '4S', '3P', '2P', '1S'],
    ('ca', 20): ['1S', '2S', '1S', '2P', '3P', '4S', '3P'],
    ('fe', 26): ['5D', '6D', '5D', '6S', '5D', '4F', '3F'],
}
LVAL = {'S': 0, 'P': 1, 'D': 2, 'F': 3, 'G': 4}


def term_weight(term):
    return int(term[0])*(2*LVAL[term[1]] + 1)


def read_elvlc(path):
    """CHIANTI 8+ .elvlc, fixed format (i7,a30,a5,i5,a5,f5.1,2f15.3)."""
    rows = []
    with open(path) as f:
        for line in f:
            if line[:7].strip() == '-1':
                break
            idx = int(line[0:7])
            conf = line[7:37].strip()
            mult = int(line[42:47])
            lchar = line[47:52].strip()
            j = float(line[52:57])
            eobs = float(line[57:72])
            eth = float(line[72:87])
            rows.append((idx, conf, mult, lchar, j, eobs, eth))
    return rows


def ground_from_elvlc(rows):
    _, conf0, m0, l0, _, _, _ = rows[0]
    g = 0
    dE = None
    for (_, conf, m, l, j, eobs, eth) in rows:
        same = conf == conf0 and m == m0 and l == l0
        if same:
            g += int(round(2*j + 1))
        else:
            e = eobs if eobs > 0.0 else eth
            if dE is None or e < dE:
                dE = e
    return g, f'{m0}{l0}', (dE*CM_TO_EV if dE is not None else -1.0)


def main():
    dbase = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.expanduser('~/RT_Codes/CHIANTI/dbase')
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(here, '..', '..', 'data', 'atomic', 'ground_term_weights.txt')
    lines = []
    nchk = 0
    for (el, Z), terms in GROUND.items():
        for s, term in enumerate(terms, start=1):
            g = term_weight(term)
            dE = -1.0
            path = os.path.join(dbase, el, f'{el}_{s}', f'{el}_{s}.elvlc')
            if os.path.exists(path):
                gc, tc, dE = ground_from_elvlc(read_elvlc(path))
                if gc != g or tc != term:
                    sys.exit(f'{el} {s}: NIST ground term {term} (g={g}) '
                             f'disagrees with CHIANTI {tc} (g={gc})')
                nchk += 1
            lines.append(f'GROUND {el:2s} {Z:3d} {s:2d} {g:3d} {term} {dE:12.6f}')
    with open(out, 'w') as f:
        f.write('# MoCHII ground-term statistical weights of the metal stages\n')
        f.write('# GROUND element Z stage g (2S+1)L  E_1[eV]\n')
        f.write('#   g: (2S+1)(2L+1) of the ground term (NIST ASD ground states);\n')
        f.write(f'#   checked against CHIANTI 11.0.2 .elvlc for {nchk} of '
                f'{len(lines)} stages (all agree)\n')
        f.write('#   E_1: first excited term of the stage, CHIANTI 11.0.2 .elvlc '
                '(observed, else theoretical); -1 = no .elvlc (not available)\n')
        f.write('#   read by src/milne_recomb_spectrum_mod.f90 (Milne factor '
                'g(X^i)/(2 g(X^(i+1))) of the metal ground continua)\n')
        f.write(f'# generated by tools/fitting/make_ground_term_weights.py on '
                f'{datetime.date.today().isoformat()}\n')
        f.write('\n'.join(lines) + '\n')
    print(f'wrote {out}: {len(lines)} stages, {nchk} checked against CHIANTI')


if __name__ == '__main__':
    main()
