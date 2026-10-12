# MoCHII: copied from EXHALE/cooling_data/chianti_cooling.py (2026-07-11); origin tree is read-only reference
"""Minimal, dependency-light CHIANTI reader + collisional line-cooling builder.

Reads the raw CHIANTI ASCII ion files (.elvlc, .wgfa, .scups) directly --- no
ChiantiPy required --- and computes the optically-thin collisional line-cooling
function Lambda_ion(T) [erg cm^3 s^-1] in stationary electron-collision
and spontaneous-radiation balance as electron density tends to zero.
Radiatively stable excited states retain collision-determined populations;
excitation and subsequent cascades contribute their emitted photon energy.
The Maxwellian excitation-rate coefficient is

    q_lu(T) = 8.629e-6 / (g_l sqrt(T)) * Upsilon_lu(T) * exp(-dE_lu/kT)

with effective collision strengths descaled from the .scups file.

Provenance: CHIANTI v11 database at $XUVTOP (here
/home/kiseon/RT_Codes/CHIANTI/dbase). This script is the auditable source for
every cooling coefficient that gets ported into the Fortran (Phase 2).
"""

import os
import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator
from functools import lru_cache
from atomic_excitation_energies import selected_energies

# Physical constants (CODATA / NIST)
K_B_ERG = 1.380649e-16        # erg / K
RY_ERG = 2.1798723611035e-11  # 1 Rydberg in erg
EV_ERG = 1.602176634e-12      # erg / eV
KB_OVER_RY = K_B_ERG / RY_ERG  # Rydberg per Kelvin = 6.33363e-6
COLL_PREF = 8.629e-6          # cm^3 s^-1 K^1/2 (Maxwellian rate prefactor)

DBASE = "/home/kiseon/RT_Codes/CHIANTI/dbase"

# Below these adopted bounds, hold Upsilon constant; Maxwellian rates and
# detailed balance still use actual gas T. Ne I's 1000.36 K and Si I's
# 1003.02 K are inferred from rounded first finite CHIANTI BT nodes, not
# published physical validity limits. Their 1,000--100,000 K use is an
# explicit numerical atomic-data approximation.
# O IV's source interval is 3200..3.2e7 K (Liang, Badnell, and Zhao 2012,
# Section 3, 200*(q+1)^2..2e6*(q+1)^2 K with q=3). Below 3200 K,
# constant Upsilon is the adopted cold approximation, not a source result.
# O V's distributed Stout/CHIANTI finite collision grid starts at 5000 K.
# The final 2014 A104 paper prints a different general sequence interval;
# 5000 K is the identified table boundary, not a claimed published accuracy limit.
BT_SOURCE_TMIN_K = {('c', 4): 3200., ('mg', 3): 1800.,
                   ('ne', 1): 1000.36, ('si', 1): 1003.02,
                   ('si', 5): 5000., ('o', 4): 3200., ('o', 5): 5000.}


def ion_dir(elem, ion):
    """CHIANTI directory for e.g. ('fe', 2) -> .../fe/fe_2."""
    return os.path.join(DBASE, elem, f"{elem}_{ion}")


def read_elvlc(path):
    """Return {level_index: (g, E_obs_cm1, E_th_cm1)}.

    Read the first 87 characters using the CHIANTI v8+ fixed fields, as in
    CHIANTI pro/read/read_elvlc.pro: (i7,a30,a5,i5,a5,f5.1,f15.3,f15.3).
    Additional comma-separated or whitespace-separated columns are metadata.
    """
    levels = {}
    with open(path) as fh:
        for number, line in enumerate(fh, 1):
            if line.strip() in ("-1", ""):
                break
            try:
                idx = int(line[:7])
                jval = float(line[52:57])
                e_obs = float(line[57:72])
                e_th = float(line[72:87])
            except ValueError as error:
                raise ValueError(f"{path}:{number}: invalid CHIANTI energy record") from error
            if idx in levels or idx <= 0 or jval < 0 or not np.all(np.isfinite([jval,e_obs,e_th])):
                raise ValueError(f"{path}:{number}: invalid or duplicate energy level")
            levels[idx] = (2.0*jval+1.0, e_obs, e_th)
    return levels


def read_wgfa(path):
    """Return list of (lower, upper, wavelength_A, gf, A_s^-1)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.rstrip("\n")
            if s.strip() in ("-1", ""):
                break
            toks = s.split()
            if len(toks) < 5:
                continue
            try:
                ll = int(toks[0]); ul = int(toks[1])
                wl = float(toks[2]); gf = float(toks[3]); a = float(toks[4])
            except ValueError:
                continue
            rows.append((ll, ul, wl, gf, a))
    return rows


def read_scups(path):
    """Parse the CHIANTI .scups file (scaled effective collision strengths).

    Returns list of dicts with keys:
      ll, ul       lower / upper level index
      de           transition energy [Rydberg]
      gf           weighted oscillator strength
      lim          high-T limit (-1 if unused)
      ttype        Burgess-Tully transition type (1..6)
      cups         scaling parameter C
      xs, ys       scaled temperature and scaled effective collision strength
    """
    trans = []
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    i = 0
    n = len(lines)
    while i < n:
        head = lines[i].split()
        if len(head) < 8 or head[0] == "-1":
            i += 1
            continue
        try:
            ll = int(head[0]); ul = int(head[1])
            de = float(head[2]); gf = float(head[3]); lim = float(head[4])
            nt = int(head[5]); ttype = int(head[6]); cups = float(head[7])
        except (ValueError, IndexError):
            i += 1
            continue
        xs = np.array(lines[i + 1].split(), dtype=float)
        ys = np.array(lines[i + 2].split(), dtype=float)
        if xs.size == nt and ys.size == nt:
            trans.append(dict(ll=ll, ul=ul, de=de, gf=gf, lim=lim,
                              ttype=ttype, cups=cups, xs=xs, ys=ys))
        i += 3
    return trans


# Ions CHIANTI 11.0.2 lacks (no valence .scups), taken from Cloudy c25.00
# Stout. Every other ion is read from CHIANTI.
STOUT_EXCITATION_IONS = {('cl', 7), ('ar', 2), ('ca', 4)}


def read_excitation_atom(elem, ion):
    """Read raw source levels, collision nodes, and identified radiative data.

    Bound-state-only Si I and Mg III radiation corrections are applied after
    the physical bound-state selection in ``bound_excitation_atom``.
    """
    if (elem.lower(), ion) in {('ar', 2), ('ca', 4)}:
        from stout_fine_structure_atoms import stout_fine_structure_atom
        return stout_fine_structure_atom(elem, ion)
    if elem.lower()=='cl' and ion==7:
        from stout_excitation_cooling import clvii_excitation_atom
        return clvii_excitation_atom()
    stem=os.path.join(ion_dir(elem,ion),f'{elem}_{ion}')
    collisions = read_scups(stem+'.scups')
    for transition in collisions:
        # A Burgess-Tully ordinate of types 1-5 is a scaled collision strength
        # and cannot be negative. CHIANTI 11.0.2 carries a few slightly
        # negative values, all at the last node (the infinite-energy limit,
        # x = 1) of type-2 rows: Mg I (2,12) -0.0095 and (13,15) -0.0007,
        # Fe II (1,9) -0.0013, Fe VII 25 rows between -1e-7 and -1.1e-4. The
        # high-energy limit of a collision strength is nonnegative, so such a
        # node is set to zero, its nearest physical value; the original value
        # is kept in 'negative_ordinates' for the file header. Below 1e5 K the
        # last interval is reached only by transitions with small scaling
        # energies. A negative ordinate anywhere else is refused.
        ys = transition['ys']
        if transition['ttype'] != 6 and np.any(ys < 0):
            negative = np.flatnonzero(ys < 0)
            if np.any(negative != ys.size - 1):
                raise ValueError(f'{elem}_{ion}: negative scaled collision strength '
                                 f'before the last node of {transition["ll"]}-{transition["ul"]}')
            transition['negative_ordinates'] = [(int(k), float(ys[k])) for k in negative]
            transition['ys'] = np.where(ys < 0, 0.0, ys)
    bt_tmin = BT_SOURCE_TMIN_K.get((elem.lower(), ion))
    if bt_tmin is not None:
        for transition in collisions:
            transition['bt_tmin_K'] = bt_tmin
    levels = read_elvlc(stem+'.elvlc')
    radiation = read_wgfa(stem+'.wgfa')
    if (elem.lower(), ion) == ('ne', 1):
        from neon_spontaneous_probabilities import nei_identified_spontaneous_probabilities
        radiation = nei_identified_spontaneous_probabilities(
            levels, radiation, collisions, stem)
    if (elem.lower(), ion) in {('o', 2), ('o', 3)}:
        from electric_dipole_probabilities import observed_energy_electric_dipole_probabilities
        radiation = observed_energy_electric_dipole_probabilities(
            ion, levels, radiation, collisions, stem)
    return levels, radiation, collisions


# Ions whose radiation is oriented by their own checked transforms
# (electric_dipole_probabilities.py, silicon_spontaneous_probabilities.py).
OWN_RADIATION_ORIENTATION = {('o', 2), ('o', 3), ('mg', 3), ('si', 1)}
# Record of the last bound_excitation_atom call for each ion, for file
# headers: number of source levels, and the rows reoriented by
# reversed_electric_dipole_probabilities (count, unmatched rows).
BOUND_ATOM_RECORD = {}


def chianti_configurations(elem, ion):
    """Configuration labels of the CHIANTI .elvlc levels, by level index."""
    configurations = {}
    with open(os.path.join(ion_dir(elem, ion), f'{elem}_{ion}.elvlc')) as fh:
        for line in fh:
            if line.strip() in ('-1', ''):
                break
            configurations[int(line[:7])] = line[7:37]
    return configurations


def bound_excitation_atom(elem, ion):
    """Collision-connected bound states for ordinary spontaneous line emission.

    Above-threshold states require autoionization/radiative branching and are
    excluded from a bound-only atom. States unreachable from the ground by
    collisions and downward cascades have exactly zero stationary population.
    Keep source transition identities and their original scaling parameters.
    """
    raw,radiation,collisions=read_excitation_atom(elem,ion)
    levels=selected_energies(elem,ion,raw,radiation,collisions)
    from ionization_potentials import ELEMENTS,read_ionization_potentials
    z=dict(ELEMENTS).get(elem.lower())
    if z is not None:
        threshold=read_ionization_potentials()[z,z-ion+1]*EV_ERG
        levels={i:r for i,r in levels.items() if r[1]*HC_ERG_CM<threshold}
    neighbors={i:set() for i in levels}
    collisions=[tr for tr in collisions if tr['ll'] in levels and tr['ul'] in levels]
    radiation=[r for r in radiation if r[0] in levels and r[1] in levels]
    if (elem.lower(), ion) == ('mg', 3):
        from electric_dipole_probabilities import bound_mgiii_electric_dipole_probabilities
        stem=os.path.join(ion_dir(elem,ion),f'{elem}_{ion}')
        radiation=bound_mgiii_electric_dipole_probabilities(
            raw, levels, radiation, collisions, stem)
    if (elem.lower(), ion) == ('si', 1):
        from silicon_spontaneous_probabilities import bound_sii_spontaneous_probabilities
        stem=os.path.join(ion_dir(elem,ion),f'{elem}_{ion}')
        radiation=bound_sii_spontaneous_probabilities(
            raw, levels, radiation, collisions, stem)
    key=(elem.lower(),ion)
    BOUND_ATOM_RECORD[key]=dict(source_levels=len(raw),reversed_e1_rows=0,
                                reversed_e1_unmatched=[])
    if key not in OWN_RADIATION_ORIENTATION and os.path.exists(
            os.path.join(ion_dir(elem,ion),f'{elem}_{ion}.elvlc')):
        from electric_dipole_probabilities import reversed_electric_dipole_probabilities
        radiation,count,unmatched=reversed_electric_dipole_probabilities(
            raw,levels,radiation,chianti_configurations(elem,ion))
        BOUND_ATOM_RECORD[key].update(reversed_e1_rows=count,reversed_e1_unmatched=unmatched)
    for tr in collisions:
        i,j=tr['ll'],tr['ul'];neighbors[i].add(j);neighbors[j].add(i)
    for i,j,w,g,a in radiation:
        if a>0:neighbors[j].add(i)
    ground=min(levels,key=lambda i:levels[i][1]);active={ground};pending=[ground]
    while pending:
        for i in neighbors[pending.pop()]-active:
            active.add(i);pending.append(i)
    return ({i:r for i,r in levels.items() if i in active},
            [r for r in radiation if r[0] in active and r[1] in active],
            [tr for tr in collisions if tr['ll'] in active and tr['ul'] in active])


@lru_cache(maxsize=8192)
def _scaled_collision_strength_interpolator(xs, ys, logarithmic):
    """Preserve nonnegative tabulated strengths without artificial zeros.

    Keep the source cubic where its exact extrema are nonnegative. If
    nonnegative nodes produce a negative interior extremum, use PCHIP on
    those same nodes. This is an interpolation approximation, not a repair
    of unresolved physical resonances. Type-6 ordinates are logarithms and
    may be negative. Negative source nodes are retained for diagnosis;
    evaluation outside their nonnegative region is explicitly rejected.
    """
    cubic = CubicSpline(xs, ys)
    if logarithmic or np.any(np.asarray(ys) < 0):
        return cubic
    extrema = cubic.derivative().roots(extrapolate=False)
    extrema = extrema[np.isfinite(extrema)]
    if np.any(cubic(extrema) < 0):
        return PchipInterpolator(xs, ys, extrapolate=False)
    return cubic


def scaled_collision_strength_interpolator(tr):
    """Source-node interpolation shared by direct rates and native fitting."""
    return _scaled_collision_strength_interpolator(
        tuple(tr['xs']), tuple(tr['ys']), tr['ttype'] == 6)


def upsilon(tr, T):
    """Burgess & Tully (1992) descaled effective collision strength Upsilon(T).

    T may be a scalar or array [K]. de is in Rydberg.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    if not np.all(np.isfinite(T)) or np.any(T <= 0):
        raise ValueError('Excitation temperature must be finite and positive')
    if "temperature_grid" in tr:
        from stout_excitation_cooling import effective_collision_strength_temperature_table
        return effective_collision_strength_temperature_table(tr["temperature_grid"],tr["upsilon_grid"],T)
    bt_tmin = tr.get('bt_tmin_K')
    if bt_tmin is not None:
        if not np.isfinite(bt_tmin) or bt_tmin <= 0:
            raise ValueError('Invalid BT source temperature minimum')
        T_source = np.maximum(T, bt_tmin)
    else:
        T_source = T
    de = tr["de"]
    C = tr["cups"]
    tt = tr["ttype"]
    et = (KB_OVER_RY * T_source) / de     # kT_source / dE
    if tt in (1, 4):
        st = 1.0 - np.log(C) / np.log(et + C)
    else:
        st = et / (et + C)
    spl = scaled_collision_strength_interpolator(tr)
    sups = spl(np.clip(st, tr["xs"][0], tr["xs"][-1]))
    if tt == 1:
        ups = sups * np.log(et + np.e)
    elif tt == 2:
        ups = sups
    elif tt == 3:
        ups = sups / (et + 1.0)
    elif tt == 4:
        ups = sups * np.log(et + C)
    elif tt == 5:
        ups = sups / et
    elif tt == 6:
        ups = 10.0 ** sups
    else:
        raise ValueError(f"unknown BT type {tt}")
    if np.any(~np.isfinite(ups)) or np.any(ups < 0):
        raise ValueError(f'Invalid collision strength at requested temperature: '
                         f'{tr["ll"]}->{tr["ul"]}; negative source nodes '
                         'cannot supply physical rates in this region')
    return ups


def cooling_lambda(elem, ion, T, lower_levels=None, upper_levels=None):
    """Optically-thin collisional line-cooling Lambda(T) [erg cm^3 s^-1].

    Sum over .scups transitions of  q_lu * dE_lu, q_lu the Maxwellian
    excitation-rate coefficient out of the lower level. Optionally restrict to a
    set of lower and/or upper levels (e.g. {1} ground, {3,4} the 2p Ly-alpha
    doublet for H I).
    Returns (Lambda_total, per_transition) where per_transition is a list of
    (ll, ul, Lambda_array).
    """
    lev, radiation, trans = bound_excitation_atom(elem, ion)
    T = np.atleast_1d(np.asarray(T, dtype=float))
    total = np.zeros_like(T)
    per = []
    for tr in trans:
        lower, upper = excitation_levels(lev, tr)
        if lower_levels is not None and lower not in lower_levels:
            continue
        if upper_levels is not None and upper not in upper_levels:
            continue
        g_l = lev[lower][0]
        de_erg = transition_energy_erg(lev, tr)
        ups = upsilon(tr, T)
        q_lu = COLL_PREF / (g_l * np.sqrt(T)) * ups * np.exp(-de_erg / (K_B_ERG * T))
        lam = q_lu * de_erg
        total += lam
        per.append((lower, upper, lam))
    return total, per


def black1981_HI_lya(T):
    """Black (1981) H I Ly-alpha collisional-excitation cooling currently in
    Cool_coeff.f90 (coex_rate_HI): 7.5e-19/(1+sqrt(T/1e5)) exp(-118348/T)."""
    T = np.asarray(T, dtype=float)
    return 7.5e-19 / (1.0 + np.sqrt(T / 1.0e5)) * np.exp(-118348.0 / T)


HC_ERG_CM = 1.9864458571489287e-16  # h c [erg cm]; E[cm^-1]*HC_ERG_CM = E[erg]
HC_OVER_K = HC_ERG_CM / K_B_ERG  # h c / k_B [cm K], consistent with excitation energies


def _level_energy_cm1(lev):
    """Map {idx:(g,E_obs,E_th)} -> {idx:(g, E_cm1)} using observed energy when
    available (>0 or ground), else the theoretical value."""
    out = {}
    for idx, (g, e_obs, e_th) in lev.items():
        e = e_obs if (e_obs > 0.0 or idx == 1) else e_th
        out[idx] = (g, e)
    return out


def excitation_levels(levels, transition):
    """Lower/upper physical-energy identities of a symmetric collision pair."""
    i,j=transition['ll'],transition['ul']
    return (i,j) if levels[j][1]>=levels[i][1] else (j,i)


def transition_energy_erg(lev, tr):
    """Energy [erg] of a .scups transition from the selected level energies.

    The .scups header value tr['de'] is the Burgess & Tully SCALING energy (a
    theoretical value used to build the scaled temperature); it is not the
    physical transition energy and must not be used either as the radiated
    energy or in the Boltzmann factor of the excitation rate.  Detailed balance
    ties the excitation rate to the SAME level gap that appears in the
    statistical-equilibrium matrix. For most ions this is the observed .elvlc
    difference (e.g. [O II] 3726/3729: .scups 30661 cm^-1 vs observed
    26811 cm^-1, a factor 0.583 in the cooling at 8000 K).
    N IV uses the checked Stout energy vector throughout its rate model.

    Upsilon is symmetric in the pair. Return the absolute physical gap;
    callers assign excitation and de-excitation with excitation_levels.
    Exactly degenerate pairs have zero photon energy but retain collisions.
    """
    de_cm1 = lev[tr["ul"]][1] - lev[tr["ll"]][1]
    return abs(de_cm1) * HC_ERG_CM


def radiative_cascade_probabilities(levels, transitions, radiation):
    """Radiative end-state probabilities in the collision-connected atom.

    States with no outgoing
    spontaneous rate can retain a finite population as electron density tends
    to zero; their populations must be determined by electron collisions.
    Source levels unreachable from the ground through the supplied collision
    and downward radiative transitions do not belong to this atomic model.
    """
    ground = min(levels, key=lambda i: levels[i][1])
    neighbors = {i: set() for i in levels}
    decay = {i: {} for i in levels}
    for tr in transitions:
        l, u = tr['ll'], tr['ul']
        if l not in levels or u not in levels:
            raise ValueError('Collision row refers to an absent energy level')
        neighbors[l].add(u)
        neighbors[u].add(l)
    for l, u, wavelength, gf, rate in radiation:
        if rate <= 0:
            continue
        if l not in levels or u not in levels or u == l:
            raise ValueError('Radiative row must connect distinct supplied levels')
        decay[u][l] = decay[u].get(l, 0.0) + rate
        neighbors[u].add(l)
    active = {ground}
    pending = [ground]
    while pending:
        i = pending.pop()
        for j in neighbors[i] - active:
            active.add(j)
            pending.append(j)
    closed = sorted(i for i in active if not decay[i])
    # With one radiative end state, total cascade energy is the excitation
    # energy independently of branching. No intermediate level gaps enter
    # this total-cooling calculation. Detailed line emission still requires
    # consistent intermediate energies, which some source models lack.
    if len(closed) == 1:
        return closed, {i: np.ones(1) for i in active}
    for i in active:
        if any(levels[i][1] <= levels[j][1] for j in decay[i]):
            raise ValueError('Multiple closed states require consistent downward radiative energies')
    ordered = sorted(active, key=lambda i: levels[i][1])
    index = {i: k for k, i in enumerate(closed)}
    probability = {}
    for i in ordered:
        value = np.zeros(len(closed))
        if i in index:
            value[index[i]] = 1.0
        else:
            total = sum(decay[i].values())
            for j, rate in decay[i].items():
                value += (rate/total)*probability[j]
        if not np.isclose(value.sum(), 1.0, rtol=1e-12):
            raise ValueError('Incomplete radiative end-state probabilities')
        probability[i] = value
    return closed, probability


def cooling_low_density_limit(elem, ion, T, restrict_to_wgfa=False):
    """Steady collisional line cooling as n_e -> 0 [erg cm^3 s^-1].

    No photoexcitation or ionization of excited levels is included. First
    follow spontaneous cascades to states with no radiative loss. Collisions
    between those states, including cascades following excitation, determine
    their stationary populations. Sum the photon energy of each cascade with
    those populations. If the ground is the only radiatively closed state,
    this reduces to sum_u q_ground,u * (E_u-E_ground).

    C III and N IV contain a radiatively closed excited state in the supplied
    data; assigning every ion to the ground is not their steady low-density
    limit. A small positive rate, including a source placeholder, still drains
    a state in the strict limit. Such a rate may imply a critical density far
    below a practical finite-density calculation; validate that separately.
    """
    levels, radiation, transitions = bound_excitation_atom(elem, ion)
    if restrict_to_wgfa:
        pairs = {(l, u) for l, u, wavelength, gf, rate in radiation if gf != 0.0 or rate != 0.0}
        transitions = [tr for tr in transitions if (tr['ll'], tr['ul']) in pairs]
    return cooling_low_density_from_levels(levels, radiation, transitions, T)


def cooling_low_density_from_levels(levels, radiation, transitions, T):
    """Strict stationary low-density cooling of identified supplied levels.

    Accept source-independent energies, spontaneous probabilities, and
    effective-strength curves. Radiatively closed states retain their
    collision-determined populations as electron density tends to zero.
    """
    closed, probability = radiative_cascade_probabilities(levels, transitions, radiation)
    T = np.atleast_1d(np.asarray(T, dtype=float))
    if not np.all(np.isfinite(T)) or np.any(T <= 0):
        raise ValueError('Cooling temperatures must be finite and positive')
    energies = np.array([levels[i][1]*HC_ERG_CM for i in closed])
    index = {i: k for k, i in enumerate(closed)}
    n = len(closed)
    rates = np.zeros((T.size, n, n))  # rates into destination from source
    photon_power = np.zeros((T.size, n))
    for tr in transitions:
        l, u = excitation_levels(levels, tr)
        if l not in index and u not in index:
            continue
        de = transition_energy_erg(levels, tr)
        strength = upsilon(tr, T)
        common = COLL_PREF*strength/np.sqrt(T)
        directions = [(l, u, common/levels[l][0]*np.exp(-de/(K_B_ERG*T))),
                      (u, l, common/levels[u][0])]
        for source, target, coefficient in directions:
            if source not in index:
                continue
            column = index[source]
            for destination in range(n):
                if destination == column:
                    continue
                transferred = coefficient*probability[target][destination]
                rates[:,destination,column] += transferred
                rates[:,column,column] -= transferred
            released = levels[target][1]*HC_ERG_CM - np.dot(probability[target],energies)
            photon_power[:,column] += coefficient*max(released, 0.0)
    if n == 1:
        return photon_power[:,0]
    matrix = rates.copy()
    scale = np.max(np.abs(matrix), axis=2)
    matrix /= np.where(scale > 0, scale, 1.0)[:,:,None]
    matrix[:,0,:] = 1.0
    rhs = np.zeros((T.size,n)); rhs[:,0] = 1.0
    try:
        populations = np.linalg.solve(matrix, rhs[...,None])[...,0]
    except np.linalg.LinAlgError as error:
        raise ValueError('Nonunique collision balance among radiatively closed states') from error
    if not np.all(np.isfinite(populations)) or np.any(populations < -64*np.finfo(float).eps):
        raise ValueError('Invalid stationary populations of radiatively closed states')
    populations = np.maximum(populations,0.0)
    populations /= populations.sum(axis=1)[:,None]
    residual = np.einsum('tij,tj->ti',rates,populations)
    rate_scale = np.max(np.abs(rates),axis=(1,2))
    if np.any(np.max(np.abs(residual),axis=1) > 1e-10*rate_scale):
        raise ValueError('Closed-state collision balance failed')
    return np.sum(populations*photon_power,axis=1)


# Electron density [cm^-3] of the reference (first) column of the runtime
# cooling tables (nlevel_cooling_mod, NE_LO). The cooling fits are made to the
# stationary state at this density (fit_cooling.py, mocassin_cooling_fits.py).
REFERENCE_DENSITY = 1.0e-6
# Below this cooling coefficient [erg cm^3 s^-1] the linear finite-density
# solve has lost relative precision to underflow of its populations.
STATIONARY_COOLING_FLOOR = 1.0e-280


def cooling_stationary_at_density(levels, radiation, collisions, T, density):
    """Cooling coefficient [erg cm^3 s^-1] of the full stationary level
    balance at a finite electron density: electron collisions and spontaneous
    decay, optically thin, no photoexcitation; the same physics as the
    installed n-level atom (nlevel_mod::nlevel_cooling_at), solved with one
    step of iterative refinement as nlevel_mod::level_populations.

    Unlike cooling_low_density_from_levels (the strict n_e -> 0 limit), a
    level whose only spontaneous decay is far slower than its collisional
    rates at this density is collisionally coupled here, not drained.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    ids = sorted(levels)
    position = {k: i for i, k in enumerate(ids)}
    n = len(ids)
    radiative = np.zeros((n, n))
    power = np.zeros(n)
    for row in radiation:
        lower, upper, rate = row[0], row[1], row[4]
        if rate <= 0:
            continue
        i, j = position[lower], position[upper]
        radiative[i, j] += rate
        radiative[j, j] -= rate
        power[j] += rate*(levels[upper][1]-levels[lower][1])*HC_ERG_CM
    lower_index, upper_index, up, down = [], [], [], []
    for tr in collisions:
        lower, upper = excitation_levels(levels, tr)
        gap = transition_energy_erg(levels, tr)
        common = COLL_PREF*upsilon(tr, T)/np.sqrt(T)
        lower_index.append(position[lower]); upper_index.append(position[upper])
        up.append(density*common/levels[lower][0]*np.exp(-gap/(K_B_ERG*T)))
        down.append(density*common/levels[upper][0])
    lower_index = np.array(lower_index, dtype=int); upper_index = np.array(upper_index, dtype=int)
    up = np.array(up).reshape(len(lower_index), T.size)
    down = np.array(down).reshape(len(lower_index), T.size)
    cooling = np.zeros(T.size)
    for k in range(T.size):
        matrix = radiative.copy()
        np.add.at(matrix, (upper_index, lower_index), up[:, k])
        np.add.at(matrix, (lower_index, lower_index), -up[:, k])
        np.add.at(matrix, (lower_index, upper_index), down[:, k])
        np.add.at(matrix, (upper_index, upper_index), -down[:, k])
        scale = np.max(np.abs(matrix), axis=1)
        matrix /= np.where(scale > 0, scale, 1)[:, None]
        matrix[0, :] = 1
        rhs = np.zeros(n); rhs[0] = 1
        population = np.linalg.solve(matrix, rhs)
        population += np.linalg.solve(matrix, rhs-matrix@population)
        cooling[k] = population@power/density
    return cooling


def log_ground_state_cooling_from_levels(levels, radiation, collisions, T):
    """Natural log of the strict-limit cooling when the ground is the only
    radiatively closed state (sum of ground excitations), evaluated in
    logarithms so Boltzmann tails below the floating-point range survive."""
    from scipy.special import logsumexp
    closed, _ = radiative_cascade_probabilities(levels, collisions, radiation)
    ground = min(levels, key=lambda k: levels[k][1])
    if closed != [ground]:
        raise ValueError('Log ground-state cooling requires one radiative end state')
    T = np.atleast_1d(np.asarray(T, dtype=float))
    terms = []
    for tr in collisions:
        lo, up = excitation_levels(levels, tr)
        if lo != ground:
            continue
        energy = transition_energy_erg(levels, tr)
        if energy <= 0:
            continue
        terms.append(np.log(COLL_PREF*upsilon(tr, T)*energy/levels[ground][0])
                     - .5*np.log(T) - energy/(K_B_ERG*T))
    if not terms:
        raise ValueError('No radiating ground excitation')
    return logsumexp(terms, axis=0)


def reference_log_cooling(levels, radiation, collisions, T):
    """Natural log of the cooling the fits are made to: the stationary state
    at REFERENCE_DENSITY.

    Not the strict n_e -> 0 limit. A source can carry a placeholder
    spontaneous rate (CHIANTI 11.0.2 S V 3s3p 3P0 -> 3s2 1S0, A = 1e-30 s^-1,
    a one-photon J = 0 -> 0 decay that is strictly forbidden); the strict
    limit lets it drain the level and radiate every excitation into it,
    whereas at any physical density (n_e above ~A/q ~ 1e-22 cm^-3) the level
    is collisionally de-excited. The runtime takes the fit as the cooling of
    the n_e = 1e-6 cm^-3 table column and scales it by the table ratio
    (nlevel_cooling_mod::nlevel_cooling_apply), so that column is the
    consistent target.

    Where the finite-density cooling falls below STATIONARY_COOLING_FLOOR
    (Boltzmann tails of highly charged ions at low T), the logarithmic
    strict-limit ground-state sum is used. That is valid there because every
    level then has a negligible population whatever its decay rate; the
    returned dictionary records the relative difference of the two at the
    lowest temperature where both are available, so a discontinuity is seen.
    Returns (log cooling, info).
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    plateau = cooling_stationary_at_density(levels, radiation, collisions, T, REFERENCE_DENSITY)
    valid = plateau > STATIONARY_COOLING_FLOOR
    result = np.full(T.size, -np.inf)
    result[valid] = np.log(plateau[valid])
    info = dict(points_below_floor=int((~valid).sum()), join_relative_difference=None)
    if not np.all(valid):
        log_strict = log_ground_state_cooling_from_levels(levels, radiation, collisions, T)
        result[~valid] = log_strict[~valid]
        if valid.any():
            k = int(np.argmax(valid))
            info['join_relative_difference'] = float(np.expm1(result[k]-log_strict[k]))
    return result, info


def log_cooling_ground_state_limit(elem, ion, T):
    """Natural log of stationary cooling with one radiative end state.

    Use only when every connected excited state relaxes to the ground.
    Logarithmic excitation sums preserve Boltzmann tails below the floating
    point range without adding a rate floor. Source collision strengths and
    energy assumptions are unchanged. This is not a detailed emission model.
    """
    from scipy.special import logsumexp
    stem = os.path.join(ion_dir(elem,ion), f'{elem}_{ion}')
    levels, radiation, collisions = bound_excitation_atom(elem, ion)
    closed, probability = radiative_cascade_probabilities(levels,collisions,radiation)
    ground = min(levels,key=lambda i:levels[i][1])
    if closed != [ground]:
        raise ValueError('Log ground-state cooling requires one radiative end state')
    T = np.atleast_1d(np.asarray(T,dtype=float))
    if np.any(T<=0) or not np.all(np.isfinite(T)):
        raise ValueError('Cooling temperatures must be finite and positive')
    terms = []
    for tr in collisions:
        if excitation_levels(levels,tr)[0] != ground:
            continue
        gap = transition_energy_erg(levels,tr)
        if gap<=0:
            raise ValueError('Ground excitation must have a positive physical gap')
        strength = upsilon(tr,T)
        with np.errstate(divide='ignore'):
            terms.append(np.log(COLL_PREF*gap/levels[ground][0])+
                         np.log(strength)-.5*np.log(T)-gap/(K_B_ERG*T))
    if not terms:
        raise ValueError('No supplied ground-state excitation for logarithmic cooling')
    return logsumexp(np.asarray(terms),axis=0)


def cooling_effective(elem, ion, T, pop="coronal", e_cut_cm1=None,
                      restrict_to_wgfa=False):
    """Effective collisional line-cooling per (n_e * n_ion_total) [erg cm^3 s^-1].

    Sum over .scups transitions of  f_l(T) * q_lu(T) * dE_lu, where f_l is the
    fractional population of the lower level l and q_lu the Maxwellian
    excitation-rate coefficient. Upper levels are assumed to decay radiatively
    (optically thin), so every excitation out of a populated lower level cools.
    dE_lu is the selected level difference (transition_energy_erg).

    pop selects the lower-level population model:
      'coronal'     imposed ground-level population (f_1 = 1). This is a
                    steady low-density limit only if all populated excited
                    states drain radiatively. Use cooling_low_density_limit
                    for the general stationary limit, including closed states.
      'ground_term' Boltzmann among the lowest levels sharing the ground config
                    label (e.g. the Fe II a6D fine-structure multiplet).  This
                    is a HIGH-density model (LTE within the ground term): it is
                    reached only for n_e >> n_crit of the fine structure, and it
                    suppresses the FIR fine-structure cooling by roughly
                    g_1/sum(g_term) when applied at nebular densities.
      'boltzmann'   Boltzmann over every level with E < e_cut_cm1 (the
                    metastable manifold); levels >= e_cut are depleted. This is
                    Huang (2023) Section 2.5 for Fe I/Fe II (cut at the lowest
                    odd-parity level that has an E1 decay to the ground term).

    restrict_to_wgfa keeps only transitions that also appear in the .wgfa file
    (i.e. carry a radiative A / oscillator strength), matching Huang's "1105 of
    4339" Fe II selection.
    """
    lev,rows,trans = bound_excitation_atom(elem,ion)
    T = np.atleast_1d(np.asarray(T, dtype=float))
    wgfa_pairs = None
    if restrict_to_wgfa:
        wgfa_pairs = {(ll,ul) for ll,ul,wl,gf,a in rows if gf != 0 or a != 0}

    # build lower-level population fractions f_l(T) (each an array over T)
    if pop == "coronal":
        lows = [1]
    elif pop == "ground_term":
        # Approximate the ground multiplet by the lowest levels within ~0.4 eV
        # of ground (e.g. the Fe II a6D fine-structure block at 0-977 cm^-1).
        thr_cm1 = 0.4 * EV_ERG / (HC_OVER_K * K_B_ERG)  # 0.4 eV -> ~3226 cm^-1
        lows = [i for i, (g, e) in lev.items() if e <= thr_cm1]
    elif pop == "boltzmann":
        if e_cut_cm1 is None:
            raise ValueError("pop='boltzmann' requires e_cut_cm1")
        lows = [i for i, (g, e) in lev.items() if e < e_cut_cm1]
    else:
        raise ValueError(f"unknown pop model {pop!r}")

    if pop == "coronal":
        f = {1: np.ones_like(T)}
    else:
        Z = np.zeros_like(T)
        for i in lows:
            g, e = lev[i]
            Z += g * np.exp(-e * HC_OVER_K / T)
        f = {i: lev[i][0] * np.exp(-lev[i][1] * HC_OVER_K / T) / Z for i in lows}

    lowset = set(lows)
    total = np.zeros_like(T)
    for tr in trans:
        ll,ul = excitation_levels(lev,tr)
        if ll not in lowset:
            continue
        if wgfa_pairs is not None and (ll,ul) not in wgfa_pairs and (ul,ll) not in wgfa_pairs:
            continue
        g_l = lev[ll][0]
        de_erg = transition_energy_erg(lev, tr)
        if de_erg <= 0.0:
            continue
        ups = upsilon(tr, T)
        q_lu = COLL_PREF / (g_l * np.sqrt(T)) * ups * np.exp(-de_erg / (K_B_ERG * T))
        total += f[ll] * q_lu * de_erg
    return total


def free_free_cooling(T, Z=1.0):
    """Huang (2023) Eq. 9 free-free cooling per (n_e * n_ion) [erg cm^3 s^-1]:
    Lambda_ff = 1.9095e-25 * Z^2 * T4^0.55,  T4 = T/1e4 K."""
    T = np.asarray(T, dtype=float)
    return 1.9095e-25 * Z**2 * (T / 1.0e4) ** 0.55


EIGHTPI_RT3 = 8.0 * np.pi / np.sqrt(3.0)  # 14.5104, Van Regemorter prefactor
RY_EV = 13.605693                          # Rydberg [eV] = hydrogen ionization I_H


def van_regemorter_clu_compact(T, dE_eV, f):
    """Compact Van Regemorter form quoted by Huang (2023) Eq. 10:

        C_lu(e) = 2.16 * alpha^-1.68 * exp(-alpha) * T^-3/2 * f ,  alpha = dE/kT.

    WARNING: cross-checking this exact expression against CHIANTI for the Mg II
    and Ca II resonance doublets shows it underestimates the electron-impact
    excitation rate by 3-8x (worse at low T; see cooling_data docs). Huang
    applies it only to the ~1025 permitted Fe I lines that have no CHIANTI
    collision data. Kept for that Fe I use and for documentation; for strong
    resonance lines prefer van_regemorter_rate().
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    alpha = (dE_eV * EV_ERG) / (K_B_ERG * T)
    return 2.16 * alpha ** (-1.68) * np.exp(-alpha) * T ** (-1.5) * f


def van_regemorter_rate(T, dE_eV, gf, g_l, gbar=0.2):
    """Rigorous Van Regemorter (1962) excitation rate coefficient [cm^3 s^-1]:

        Upsilon = (8 pi/sqrt 3) * gf * (I_H/dE) * gbar
        C_lu    = 8.629e-6 / (g_l sqrt(T)) * Upsilon * exp(-dE/kT)

    gf = g_l * f_absorption is the weighted oscillator strength; gbar is the
    effective (Maxwellian-averaged) Gaunt factor. gbar ~ 0.2 for ions: this
    reproduces the CHIANTI effective collision strengths of the Mg II and Ca II
    resonance doublets to within ~20% over 5000-30000 K (back-solved gbar =
    0.19-0.24 for Mg II, 0.14-0.22 for Ca II). Used for Na I (isoelectronic with
    Mg II, no CHIANTI data) and, optionally, Fe I.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    alpha = (dE_eV * EV_ERG) / (K_B_ERG * T)
    ups = EIGHTPI_RT3 * gf * (RY_EV / dE_eV) * gbar
    return COLL_PREF / (g_l * np.sqrt(T)) * ups * np.exp(-alpha)


# Na I resonance D doublet (3s 2S1/2 -> 3p 2P_J); NIST A/f values (Kramida 2020).
#   (dE_eV, f_absorption, label).  Lower level 3s 2S1/2 has g_l = 2.
NAI_GL = 2.0
NaI_D_LINES = [
    (2.10230, 0.320, "D1 5895.92A 3s2S1/2-3p2P1/2"),
    (2.10440, 0.641, "D2 5889.95A 3s2S1/2-3p2P3/2"),
]


def cooling_NaI(T, gbar=0.2):
    """Na I D-doublet collisional line cooling per (n_e * n_NaI) [erg cm^3 s^-1].

    Coronal/optically-thin limit: every electron-impact excitation out of the
    3s ground level is followed by a radiative decay, so Lambda = sum_D C_lu*dE.
    Na I is absent from CHIANTI v11.0.2, but it is isoelectronic with Mg II
    (Na-like 3s 2S -> 3p 2P resonance doublet), so the collision rate uses the
    Van Regemorter approximation with the Gaunt factor gbar calibrated against
    CHIANTI Mg II/Ca II (~0.2). This yields Upsilon(3s-3p) ~ 36 near 1e4 K,
    consistent with published Na I R-matrix collision strengths; the anomalous
    two-electron Mg I 1S-1P line (gf=3.4 but Upsilon~2) is NOT a valid analog.
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    total = np.zeros_like(T)
    for dE_eV, f, _ in NaI_D_LINES:
        gf = NAI_GL * f
        total += van_regemorter_rate(T, dE_eV, gf, NAI_GL, gbar) * (dE_eV * EV_ERG)
    return total


if __name__ == "__main__":
    T = np.logspace(3.0, 5.0, 60)
    lya, per = cooling_lambda("h", 1, T, lower_levels={1}, upper_levels={3, 4})
    tot, _ = cooling_lambda("h", 1, T, lower_levels={1})
    blk = black1981_HI_lya(T)
    print(f"{'T[K]':>10} {'CHIANTI_Lya':>13} {'CHIANTI_tot':>13} "
          f"{'Black':>12} {'Black/CHIANTI_Lya':>18}")
    for k in range(0, len(T), 4):
        print(f"{T[k]:10.3e} {lya[k]:13.4e} {tot[k]:13.4e} "
              f"{blk[k]:12.4e} {blk[k]/lya[k]:18.2f}")
    # Cross-check the descaled Upsilon(1s-2p) against the pre-extracted file.
    d = ion_dir("h", 1)
    trans = read_scups(os.path.join(d, "h_1.scups"))
    t13 = [t for t in trans if t["ll"] == 1 and t["ul"] == 3][0]
    t14 = [t for t in trans if t["ll"] == 1 and t["ul"] == 4][0]
    for Tc in (1.0e4, 2.0e4):
        u = upsilon(t13, Tc)[0] + upsilon(t14, Tc)[0]
        print(f"Upsilon(1s-2p) at T={Tc:.0e} K : {u:.4f}")
