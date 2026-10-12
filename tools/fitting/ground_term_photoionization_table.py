#!/usr/bin/env python3
"""Tabulated R-matrix photoionization cross sections of metal ground terms.

WHAT THIS PRODUCES
------------------
The PHTABLE record of data/atomic/element_<el>.txt, written by
make_element_data.py through ground_term_photoionization_table(): for the
transitions listed in TABLES, the photoionization cross section of the
ground TERM of the ion, tabulated against the photoelectron energy
eps = E - ETH [eV] from the threshold up to the energy E_join where it
meets the VFKY96 outer-shell fit (PHOTO row).  Below E_join the table
replaces the PHOTO fit for every consumer: the transport opacity, the
photoionization rate and its heating (species_mod::species_cross_section),
and the Milne ground-term recombination of the diffuse field (channel 6,
species_mod::species_outer_shell_cross_section).  At and above E_join the
phfit2 prescription is unchanged (PHOTO below E_inn, the PHSHELL sum above).

WHICH TRANSITIONS
-----------------
Ne III (Ne2+ + h nu -> Ne3+ + e), transition 3 of neon: Nahar (2019, New
Astron. 67, 97; Breit-Pauli R-matrix, close-coupling expansion of 58 fine-
structure levels of Ne IV up to 2s2 2p2 3d, radiation damping included),
the NORAD-Atomic-Data file ne3.px.gd.txt (total photoionization cross
sections of the three levels 2s2 2p4 3P2,1,0; downloaded 2026-10-11 from
https://norad.astronomy.osu.edu/ne3/ne3.px.gd.txt, stored unchanged as
data/atomic/norad/ne3.px.gd.txt, md5 below).  Checked against the
publications (references/Nahar_2019NewAstro_67_97.pdf, Nahar_2019IJMS_443_61
.pdf, read): the binding energies of the three levels, 4.71421 / 4.70872 /
4.70611 Ry, are the BPRM column of Nahar (2019, New Astron.) Table 2 (4.714,
4.709, 4.706); the 3P2 cross section reproduces that paper's Fig. 1a (6.5 Mb
background at 5.5 Ry, 2.2 Mb at 10.4 Ry, the resonance groups converging on
Ne IV 2D/2P below 5.1 Ry and the n = 3 groups at 7-10 Ry); Nahar et al. (2019,
IJMS) compare the same calculation with the ALS merged-beam measurement
(8.2 meV resolution, absolute scale +-20%) and attribute most of the
measured structure below 68 eV to these three levels.  The authors state an
accuracy of 10-15% (New Astron., Sec. 5).
Why it replaces VFKY96 here: the VFKY96 fit is the 1993 Opacity Project
LS-coupling cross section smoothed over its resonances, and that smoothing
drops most of the near-threshold 2p3(2D, 2P) n(s,d) resonance strength.  The
Milne integral of the VFKY96 fit gives 0.545 of the ground-term capture
rate of this table at 1e4 K (0.365 at 5e3 K; alpha_1 2.78e-12 against
5.10e-12 cm^3 s^-1 at 1e4 K); Badnell's RR + DR partial rates into the
ground term (adf48/adf09, an independent method) are 0.87-0.91 of this
table over 5e3-2e4 K (md/atomic_stage_expansion/excited_level_ionizing_
photons_2026-10-11/measurement/ground_capture/).

PHYSICS OF THE TABLE
--------------------
(1) Ground term, not ground level.  MoCHII carries one ground state per
    stage.  The table is the statistical average over the fine-structure
    levels, sigma = sum_J (2J+1) sigma_J / sum_J (2J+1), each level taken on
    the PHOTOELECTRON-energy scale, so all three thresholds sit at ETH (the
    NIST ionization potential of the 3P2 ground level).  This is the
    cross section whose Milne inverse with the term weight g = 9 is the sum
    of the captures into the three levels, which is what the ground-term
    weights of channel 6 assume.  The photoelectron scale is also where the
    resonances of the three levels coincide: they converge on the same
    Ne IV thresholds, which the calculation places at the observed energies
    (New Astron., Sec. 3), while its own bound-level energies are 0.7 eV too
    deep (64.14 against 63.42 eV for 3P2), so the published photon-energy
    scale is NOT used.  Validity: the three level cross sections agree on
    this scale (their integrals over eps < 40 eV within 0.5% / 0.8% of 3P2's),
    so the absorption does not depend on how the levels are populated
    (3P2 alone at low density); the thresholds of 3P1 and 3P0 (0.080 and
    0.114 eV below 3P2's) are moved up to ETH, which shifts the photon of a
    capture into those levels by at most 0.114 eV.
(2) Resonances.  The table keeps every node of the source on the 3P2 grid:
    steps of 1e-4 Ry (1.36 meV) below 0.5 Ry above threshold, 3e-4 Ry
    (4.08 meV) to 2 Ry, coarser above, where the source has fewer
    resonances; the 3P1 and 3P0 grids are identical to it up to 6.89 Ry and
    are interpolated linearly onto it above (offsets <= 4e-6 Ry).  The
    resolution is that of the source calculation; nothing is smoothed.
    MoCHII interpolates the table linearly in eps, at the exact packet
    energy, and the channel-6 Milne integral and its sampled spectrum
    integrate the same piecewise-linear function node by node.
(3) Energy range.  The source extends to eps = 363 eV, but above its
    highest target threshold (Ne IV 2s2 2p2 3d 4D7/2, 5.3109 Ry: photon
    135.7 eV) it is a background that falls below the Hartree-Dirac-Slater
    cross section VFKY96 reproduces at high energy (0.57 of it at 427 eV);
    that is outside what a 14-orbital R-matrix basis is built for, and the
    HDS single-electron limit is the better determination there.  The
    table is therefore cut at the first crossing of the two above 135.7 eV,
    E_join, found by bisection on the linearly interpolated table, and the
    last node is the VFKY96 value at E_join, so the cross section is
    continuous there.  Above E_join phfit2 is unchanged: the VFKY96 fit
    (the 2s + 2p L shell) below the K edge E_inn, the subshell sum above.
    No shell is counted twice: the R-matrix total contains the 2s channels
    (2s 2p4 targets) and the 2p channels; the PHOTO fit it replaces is the
    same L-shell total.
(4) Heating.  The table is a total over the exit channels (Ne IV 4S, 2D,
    2P, 2s 2p4 ...); NORAD has no partial cross sections for Ne III.  The
    photoelectron is therefore given E - ETH, the convention of the VFKY96
    L-shell fit it replaces; for an exit into an excited Ne IV term the
    excitation energy (>= 5.1 eV) is counted as heat instead of the line
    it radiates, an overestimate confined to E > ETH + 5.1 eV.  Channel 6
    has the same limitation: above ETH + E_1 (E_1 = 5.11 eV, Ne IV 2D) the
    total counts excited-term photoionization as ground recombination;
    milne_recomb_spectrum_mod reports that share (0.25% at 1e4 K).

Run as a module (make_element_data.py), or alone to print the table summary:
    python3 ground_term_photoionization_table.py
"""
import hashlib
import math
import os

RYDBERG_EV = 13.605693122994        # CODATA 2018 R_inf h c [eV]
HERE = os.path.dirname(os.path.abspath(__file__))
ATOMIC = os.path.join(HERE, "..", "..", "data", "atomic")

#--- (element, transition) -> source description
TABLES = {
    ("ne", 3): dict(
        path=os.path.join(ATOMIC, "norad", "ne3.px.gd.txt"),
        md5="964ecbb6669aef6f1ecabd8e61d7fd87",
        # level -> statistical weight 2J+1 (block key: 2J of the NORAD file)
        levels={4: 5, 2: 3, 0: 1},
        reference="Nahar 2019, New Astron. 67, 97 (BPRM, 58 Ne IV levels); "
                  "NORAD ne3.px.gd.txt, 3P2,1,0 averaged with 2J+1 on the "
                  "photoelectron-energy scale",
        # highest target threshold of the close-coupling expansion [Ry]
        # (Nahar 2019 New Astron. Table 1, level 58, 2s2 2p2 3d 4D7/2, NIST)
        last_target_ry=5.3109,
    ),
}


def vfky96(E, p):
    """VFKY96 eq. (1) [Mb], p = (E_th, E0, s0, ya, P, yw, y0, y1), as
    photo_xsec::sigma_vfky96."""
    if E < p[0]:
        return 0.0
    _, E0, s0, ya, P, yw, y0, y1 = p
    x = E / E0 - y0
    z = math.sqrt(x * x + y1 * y1)
    Q = 5.5 - 0.5 * P
    return s0 * ((x - 1.0) ** 2 + yw ** 2) * z ** (-Q) * (1.0 + math.sqrt(z / ya)) ** (-P)


def read_norad_levels(path):
    """{2J: (binding [Ry], eps [Ry] list, sigma with radiation damping [Mb] list)}
    of a NORAD px.gd file."""
    lines = open(path).read().split("\n")
    out = {}
    j = 0
    while j < len(lines) - 1:
        t = lines[j].split()
        u = lines[j + 1].split()
        if (len(t) == 4 and t[0] == "0" and t[2] == "0" and t[3] == "1"
                and len(u) == 3):
            twoj = int(t[1])
            be = -float(u[0])
            n = int(u[1])
            eps, sig = [], []
            for row in lines[j + 3:j + 3 + n]:
                v = row.split()
                if len(v) != 4:
                    raise ValueError(f"{path}: malformed row in level 2J={twoj}: {row!r}")
                eps.append(float(v[1]))
                sig.append(float(v[3]))
            if len(eps) != n:
                raise ValueError(f"{path}: level 2J={twoj} has {len(eps)} rows, header {n}")
            out[twoj] = (be, eps, sig)
            j += 3 + n
        else:
            j += 1
    return out


def _interp(x, xs, ys):
    """linear interpolation of the sorted table (xs, ys) at x (inside it)"""
    lo, hi = 0, len(xs) - 1
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    w = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] + w * (ys[hi] - ys[lo])


def ground_term_photoionization_table(el, transition, eth, photo):
    """(eps [eV] list, sigma [Mb] list, notes) of the PHTABLE record, or None
    when the transition has no table.  eth: the ETH row [eV]; photo: the
    VFKY96 parameters of the PHOTO row (E_th first), for the join."""
    spec = TABLES.get((el, transition))
    if spec is None:
        return None
    with open(spec["path"], "rb") as fh:
        md5 = hashlib.md5(fh.read()).hexdigest()
    if md5 != spec["md5"]:
        raise ValueError(f"{spec['path']}: md5 {md5}, expected {spec['md5']}")
    lev = read_norad_levels(spec["path"])
    if sorted(lev) != sorted(spec["levels"]):
        raise ValueError(f"{spec['path']}: levels {sorted(lev)}, expected {sorted(spec['levels'])}")
    #--- common grid: the ground level's (largest weight) photoelectron grid
    base = max(spec["levels"], key=lambda k: spec["levels"][k])
    grid = lev[base][1]
    if grid[0] != 0.0 or any(b <= a for a, b in zip(grid, grid[1:])):
        raise ValueError(f"{spec['path']}: grid of 2J={base} not increasing from 0")
    gsum = float(sum(spec["levels"].values()))
    sig = [0.0] * len(grid)
    for twoj, g in spec["levels"].items():
        _, eps_j, sig_j = lev[twoj]
        same = len(eps_j) == len(grid) and all(abs(a - b) <= 1e-12 for a, b in zip(eps_j, grid))
        for i, e in enumerate(grid):
            s = sig_j[i] if (same or abs(eps_j[i] - e) <= 1e-12) else _interp(e, eps_j, sig_j)
            sig[i] += g * s / gsum
    eps_ev = [e * RYDBERG_EV for e in grid]
    #--- the join: first crossing table = VFKY96 above the last target threshold
    e_last = spec["last_target_ry"] * RYDBERG_EV
    diff = [s - vfky96(eth + e, photo) for e, s in zip(eps_ev, sig)]
    i0 = next(i for i, e in enumerate(eps_ev) if e > e_last)
    if diff[i0] <= 0.0:
        raise ValueError(f"{el} {transition}: table already below VFKY96 at the last target threshold")
    ij = next((i for i in range(i0, len(eps_ev) - 1) if diff[i] > 0.0 >= diff[i + 1]), None)
    if ij is None:
        raise ValueError(f"{el} {transition}: table never meets VFKY96 above {e_last:.2f} eV")
    e0, e1 = eps_ev[ij], eps_ev[ij + 1]
    s0, s1 = sig[ij], sig[ij + 1]

    def d(e):
        return s0 + (s1 - s0) * (e - e0) / (e1 - e0) - vfky96(eth + e, photo)
    a, b = e0, e1
    for _ in range(200):
        m = 0.5 * (a + b)
        if d(m) > 0.0:
            a = m
        else:
            b = m
    e_join = 0.5 * (a + b)
    eps_out = eps_ev[:ij + 1] + [e_join]
    sig_out = sig[:ij + 1] + [vfky96(eth + e_join, photo)]
    notes = dict(reference=spec["reference"], md5=md5, source=os.path.relpath(spec["path"], ATOMIC),
                 e_join=eth + e_join, n=len(eps_out), e_last=e_last)
    return eps_out, sig_out, notes


if __name__ == "__main__":
    import re
    from make_element_data import ETH, parse_ph2
    ph2 = parse_ph2()
    for (el, it) in TABLES:
        Z = dict(ne=10)[el]
        ne_ion = Z - (it - 1)
        eth = ETH[(Z, ne_ion)]
        eps, sig, notes = ground_term_photoionization_table(el, it, eth, [eth] + ph2[(Z, ne_ion)])
        print(f"{el} transition {it}: {notes['n']} nodes, E_join = {notes['e_join']:.4f} eV, "
              f"sigma(ETH) = {sig[0]:.4f} Mb, max {max(sig):.2f} Mb")
