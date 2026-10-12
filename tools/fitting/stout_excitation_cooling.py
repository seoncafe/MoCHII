"""Published-temperature outer-shell Cl VII excitation cooling.

Only the 32 bound outer-shell levels are used. Core-excited Stout levels
above the physical ionization threshold require autoionization branching
and are excluded from this ordinary bound-level emission model. The tables
start at 9800 K. The adopted cold approximation holds Upsilon at its
9800 K value, while rates retain their actual Maxwellian temperature.
The numerical use interval is 1000--100000 K; low-temperature resonance
accuracy is not established by the source.
Ground-state reactants, Maxwellian electrons and n_e -> 0 are assumed.
"""
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from chianti_cooling import COLL_PREF,HC_ERG_CM,K_B_ERG

STEM=Path.home()/'CLOUDY/c25.00/data/stout/cl/cl_7/cl_7'


def clvii_outer_atom():
    """Read actual outer-shell energies, A values and effective strengths."""
    levels={};radiation=[];collisions=[];temperatures=None
    for line in Path(str(STEM)+'.nrg').read_text().splitlines()[1:]:
        if line.startswith('*'):break
        f=line.split(maxsplit=3)
        if len(f)>=3 and f[0].isdigit() and int(f[0])<=32:
            levels[int(f[0])]=(float(f[2]),float(f[1]))
    for line in Path(str(STEM)+'.tp').read_text().splitlines()[1:]:
        f=line.split()
        if f and f[0]=='A' and int(f[1])<=32 and int(f[2])<=32:
            lower,upper,rate=int(f[1]),int(f[2]),float(f[3])
            if rate<=0 or levels[upper][1]<=levels[lower][1]:raise ValueError('Invalid Cl VII radiative pair')
            radiation.append((lower,upper,rate))
    for line in Path(str(STEM)+'.coll').read_text().splitlines()[1:]:
        f=line.split()
        if f and f[0]=='TEMP':temperatures=np.array(list(map(float,f[1:])))
        if f and f[:2]==['CS','ELECTRON'] and int(f[2])<=32 and int(f[3])<=32:
            lower,upper=int(f[2]),int(f[3]);strength=np.array(list(map(float,f[4:])))
            if levels[upper][1]<=levels[lower][1] or np.any(strength<=0):raise ValueError('Invalid Cl VII collision pair')
            collisions.append((lower,upper,strength))
    if set(levels)!=set(range(1,33)) or temperatures is None or min(temperatures)!=9800:
        raise ValueError('Cl VII source changed or incomplete')
    if np.any(np.diff(temperatures)<=0) or any(u.shape!=temperatures.shape for l,h,u in collisions):
        raise ValueError('Invalid Cl VII collision grid')
    decays={i:[] for i in levels}
    for l,u,a in radiation:decays[u].append(l)
    if [i for i in levels if not decays[i]]!=[1]:raise ValueError('Cl VII outer atom requires one radiatively stable ground state')
    return levels,radiation,temperatures,collisions


def effective_collision_strength_temperature_table(grid, strength, temperature):
    """Log-PCHIP Upsilon on supplied nodes, constant below the first node.

    The Cl VII caller supplies 9800--9.8e7 K nodes. Other callers supply
    their own temperature tables; this routine does not assign the Cl VII
    interval to them. Cold continuation is the adopted approximation, not
    a published accuracy claim. Maxwellian rates use the actual temperature.
    """
    grid=np.asarray(grid,dtype=float)
    strength=np.asarray(strength,dtype=float)
    temperature=np.atleast_1d(np.asarray(temperature,dtype=float))
    if (grid.ndim!=1 or len(grid)<2 or strength.shape!=grid.shape or
        not np.all(np.isfinite(grid)) or np.any(grid<=0) or np.any(np.diff(grid)<=0) or
        not np.all(np.isfinite(strength)) or np.any(strength<=0)):
        raise ValueError('Require ordered positive temperature nodes and positive strengths')
    # A printed native endpoint can differ by a few ulps after exp(log T).
    # Permit only representational rounding, not a physical hot continuation.
    ceiling=np.nextafter(grid[-1],np.inf)
    for _ in range(7):ceiling=np.nextafter(ceiling,np.inf)
    if not np.all(np.isfinite(temperature)) or np.any(temperature<=0) or np.any(temperature>ceiling):
        raise ValueError(f'Temperature table requires finite 0 < T <= {grid[-1]:.17g} K')
    return np.exp(PchipInterpolator(np.log(grid),np.log(strength),extrapolate=False)(
        np.log(np.clip(temperature,grid[0],grid[-1]))))


def clvii_excitation_atom():
    """Bound Stout atom in the shared excitation representation.

    Type 2 supplies an invertible coordinate for native Chebyshev storage;
    source strengths are log-PCHIP tables, not Burgess-Tully source splines.
    Its lower coordinate is the 9800 K source node. Existing native clipping
    then gives constant Upsilon below 9800 K without changing Maxwellian T.
    """
    from chianti_cooling import RY_ERG,K_B_ERG
    levels,rad,grid,coll=clvii_outer_atom()
    raw={i:(g,e,e) for i,(g,e) in levels.items()}
    radiation=[(l,u,1e8/(levels[u][1]-levels[l][1]),0.,a) for l,u,a in rad]
    transitions=[]
    for l,u,strength in coll:
        de=(levels[u][1]-levels[l][1])*HC_ERG_CM/RY_ERG
        et=K_B_ERG*grid/(de*RY_ERG)
        transitions.append(dict(ll=l,ul=u,de=de,gf=0.,lim=-1.,ttype=2,cups=1.,
            xs=et/(1+et),ys=strength,
            temperature_grid=grid,upsilon_grid=strength))
    return raw,radiation,transitions


def clvii_outer_cooling(temperature):
    """Coronal cooling [erg cm^3 s^-1] with the adopted cold Upsilon limit."""
    levels,rad,grid,collisions=clvii_outer_atom()
    temperature=np.atleast_1d(np.asarray(temperature,dtype=float))
    value=np.zeros_like(temperature)
    for lower,upper,strength in collisions:
        if lower!=1:continue
        ups=effective_collision_strength_temperature_table(grid,strength,temperature)
        gap=(levels[upper][1]-levels[lower][1])*HC_ERG_CM
        value+=COLL_PREF*ups*gap/(levels[1][0]*np.sqrt(temperature))*np.exp(-gap/(K_B_ERG*temperature))
    return value
