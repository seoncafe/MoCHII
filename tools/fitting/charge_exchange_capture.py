"""Classical polarization capture bound for ground-state H and positive ions.

Ferland et al. 1997, ApJ 481 L115, Eq. (2), DOI 10.1086/310655.
Their Fe-specific Eq. (3) is rounded to 1.92e-9*q cm^3/s. Use the actual
reduced mass for other elements. Every trajectory crossing the centrifugal
barrier is assumed to transfer an electron. This is a capture-model bound,
not an independently validated electron-transfer coefficient. The user
adopted this upper-limit approximation for the ten q>=5 added boundaries.
Element generation selects it there; production activation remains separate.

Applicability: isolated ion + H(1s), a long-range induced-dipole potential,
and classical trajectories. No published temperature interval establishes
reaction accuracy for q >= 5. Diagnostics below use E=k_B*T only to show the
characteristic barrier radius and angular momentum, not to define validity.
"""
import argparse
import json
import math
from pathlib import Path

# SI constants, CODATA 2022; H(1s) static polarizability volume is 4.5*a0^3.
BOHR_M=5.29177210544e-11
ELEMENTARY_CHARGE_C=1.602176634e-19
EPSILON_0=8.8541878188e-12
ATOMIC_MASS_KG=1.66053906892e-27
BOLTZMANN_J_K=1.380649e-23
HBAR_J_S=1.054571817e-34
HYDROGEN_MASS_U=1.00782503223
# Representative natural atomic masses; electron-mass corrections are omitted.
ION_MASS_U={'c':12.0107,'n':14.0067,'o':15.9994,'ne':20.1797,
            'mg':24.305,'si':28.0855,'s':32.065,'ar':39.948,
            'fe':55.845,'cl':35.453,'ca':40.078}
LANGEVIN_BOUNDARIES=(('n',5),('s',5),('s',6),('ar',5),('ar',6),
                    ('fe',5),('fe',6),('cl',5),('cl',6),('ca',5),
                    ('c',5),('c',6),('n',6),('o',5),('o',6),
                    ('ne',5),('ne',6),('mg',5),('mg',6),
                    ('si',5),('si',6),('ca',6))


def hydrogen_ion_reduced_mass(ion_mass_u):
    """Reduced collision mass [kg] for H(1s) and the specified ion mass."""
    if not math.isfinite(ion_mass_u) or ion_mass_u<=0:
        raise ValueError('Ion mass must be finite and positive')
    return HYDROGEN_MASS_U*ion_mass_u/(HYDROGEN_MASS_U+ion_mass_u)*ATOMIC_MASS_KG


def polarization_potential_coefficient(charge):
    """C4 [J m^4] in U(R)=-C4/R^4, using the H(1s) polarizability."""
    if isinstance(charge,bool) or not isinstance(charge,int) or charge<=0:
        raise ValueError('Ion charge must be a positive integer')
    alpha_volume=4.5*BOHR_M**3
    return alpha_volume*(charge*ELEMENTARY_CHARGE_C)**2/(8*math.pi*EPSILON_0)


def langevin_hydrogen_capture_bound(charge,ion_mass_u):
    """Classical unit-transfer capture bound [cm^3/s], independent of T.

    k=2*pi*sqrt(2*C4/mu). A physical rate additionally requires a reaction
    probability, product distribution, and assessment of the potential model.
    """
    mu=hydrogen_ion_reduced_mass(ion_mass_u)
    return 2*math.pi*math.sqrt(2*polarization_potential_coefficient(charge)/mu)*1e6


def polarization_capture_conditions(charge,ion_mass_u,temperature_K):
    """Characteristic barrier diagnostics at relative energy E=k_B*T.

    Long-range and classical conditions require R/a0 >> 1 and L/hbar >> 1;
    these are not sufficient to establish unit electron-transfer probability.
    """
    if not math.isfinite(temperature_K) or temperature_K<=0:
        raise ValueError('Temperature must be finite and positive')
    c4=polarization_potential_coefficient(charge)
    mu=hydrogen_ion_reduced_mass(ion_mass_u);energy=BOLTZMANN_J_K*temperature_K
    return {'temperature_K':temperature_K,'barrier_radius_a0':(c4/energy)**.25/BOHR_M,
            'critical_angular_momentum_hbar':2*math.sqrt(mu)*(c4*energy)**.25/HBAR_J_S}


def langevin_recombination_rates():
    """Selected unit-probability H(1s) capture rates; no certified T interval."""
    return {key:langevin_hydrogen_capture_bound(key[1],ION_MASS_U[key[0]])
            for key in LANGEVIN_BOUNDARIES}


def avoided_crossing_mean_heat(charge):
    """Mean kinetic release [eV] at a 9.5-bohr avoided crossing, q>4.

    The exit channel has ion charges +1 and +(q-1), hence its Coulomb
    energy is (q-1)/9.5 hartree, rounded to 2.86*(q-1) eV. KF99 Section 2
    prints 2.86*q despite defining q as the initial charge; that expression
    is inconsistent with this exit-channel potential. Cloudy c25 uses q-1.
    Ground-state reactants and low-density radiative relaxation are assumed.
    Product branching, cascade reabsorption, collisional redistribution, the
    reverse reaction and a physical temperature accuracy range are unresolved.
    """
    if isinstance(charge,bool) or not isinstance(charge,int) or charge<=4:
        raise ValueError('KF99 high-charge heat requires integer q>4')
    return 2.86*(charge-1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    rows=[]
    for element,charge in LANGEVIN_BOUNDARIES:
        mass=ION_MASS_U[element];rate=langevin_hydrogen_capture_bound(charge,mass)
        rows.append({'element':element,'recombining_charge':charge,'ion_mass_u':mass,
                     'capture_bound_cm3_s':rate,'ratio_to_rounded_Fe_rule':rate/(1.92e-9*charge),
            'mean_kinetic_heat_eV':avoided_crossing_mean_heat(charge),
            'conditions':[polarization_capture_conditions(charge,mass,t)
                                   for t in (1000.,10000.,100000.)]})
    result={'source_doi':'10.1086/310655','status':'user-selected capture upper-limit approximation; production activation separate',
            'temperature_accuracy_interval':None,'unit_transfer_probability':'assumed, not verified',
            'rows':rows}
    text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    else:print(text,end='')

if __name__=='__main__':main()
