# Fe V--VII recombination source audit

## Physical judgment

The available files do not establish a physically reliable low-temperature
dielectronic recombination (DR) prescription for all of Fe V, Fe VI, and Fe VII.
The production generator's Fe-stage guard must remain in place. In particular,
using the old CHIANTI DR fits at nebular temperatures would effectively omit
low-temperature resonances, and entering zero DR would assert a false physical
result. The Cloudy GuPC fit covers Fe VII to Fe VI only and is labeled as a
private communication. It cannot fill Fe V or Fe VI.

## Charge indexing and file contents

CHIANTI `fe_s` RR/DR data specify recombination of Fe `s` into Fe `s-1`.
For atomic number 26, the recombining ions Fe V, Fe VI, and Fe VII have 22,
21, and 20 bound electrons respectively. In the production generator, these
would be transitions `i=4,5,6`, with `Nrec=26-i` and CHIANTI stage `i+1`.
The local Badnell RR and DR tables contain no `M=1` record for these three
recombining ions. Cloudy's `Fe_Gu_c[0]` and `Fe_Gu_E[0]` are labeled Fe+6,
and the loop stores them at `ion+5` with `ion=0`, where Cloudy's zero-based
index 5 means Fe+6 recombining to Fe+5. Thus the fit is Fe VII to Fe VI.

| Recombination | CHIANTI RR fit | CHIANTI DR fit | Cloudy GuPC |
| --- | --- | --- | --- |
| Fe V to Fe IV | Type 3, interpolated from several older sources | Type 1, file cites Mazzotta et al. (1998) | No row |
| Fe VI to Fe V | Type 3, same source description | Type 1, same citation | No row |
| Fe VII to Fe VI | Type 3, same source description | Type 1, same citation | First row, Gu private communication |

The local CHIANTI files (`VERSION` = `11.0.2`) describe RR as Ken Dere's 2007
interpolation from Shull and Van Steenberg (1982), Landini and Monsignori Fossi
(1990, 1991), and Mazzitelli and Mattioli (2002). They cite Mazzotta et al.
(1998) for DR. These are *file provenance statements*, not independently
verified claims about the final journal papers. The RR type-3 form is
`A (T/10^4 K)^(-eta)` in `cm^3 s^-1`. DR type 1 is
`T^(-3/2) sum_j c_j exp(-E_j/T)`, with `E_j` in K and `c_j` in
`cm^3 s^-1 K^(3/2)`. The Cloudy GuPC form uses the latter expression with
temperature and `E_j` in eV; its coefficient units therefore contain
`eV^(3/2)` rather than `K^(3/2)`. Directly copying the Cloudy coefficients
into a Kelvin-based DR record would be incorrect.

## Measured candidate rates

The accompanying audit script evaluated the local files directly. Values are
in `cm^3 s^-1` at `10^4 K` and are source-fit values, not recommended
production coefficients.

| Recombination | CHIANTI RR | CHIANTI DR | Cloudy GuPC DR |
| --- | ---: | ---: | ---: |
| Fe V to Fe IV | 7.8000e-12 | 6.0480e-27 | unavailable |
| Fe VI to Fe V | 1.5100e-11 | 3.8722e-35 | unavailable |
| Fe VII to Fe VI | 2.6200e-11 | 1.0797e-30 | 5.8929e-11 |

The full source-fit curves on a 301-point logarithmic grid from `10^3` to
`10^6 K` are saved in [the JSON](fe_recombination_expansion/fe_recombination_curves.json)
and [the vector PDF](fe_recombination_expansion/fe_recombination_curves.pdf).
RR and DR are shown separately; the GuPC curve appears only for Fe VII.
All curves are marked as having unverified temperature validity.

At `10^3 K`, the three CHIANTI DR rates are respectively `1.2465e-194`,
`1.7762e-279`, and `1.4340e-235 cm^3 s^-1`; Cloudy GuPC gives
`3.7786e-10 cm^3 s^-1` for Fe VII. This large source dependence shows why
the CHIANTI high-temperature fits alone cannot settle the nebular DR choice.
Neither file states a demonstrated temperature validity interval for the
selected Fe fits. Numerical evaluation outside a documented interval is not
physical validation.

| Comparison condition | CHIANTI file evaluation | Cloudy GuPC evaluation |
| --- | --- | --- |
| Recombining ion | Fe V, VI, VII | Fe VII only |
| Electron temperature | Same 301 nodes from `10^3` to `10^6 K` | Same nodes, converted to eV |
| Density and radiation field | Rate coefficient files contain no dependence | Local expression contains no dependence |
| RR contribution | Evaluated separately | Not part of GuPC DR fit |
| Boundary conditions, gravity, species chemistry | No plasma model evaluated | No plasma model evaluated |

This is a coefficient comparison; it produces no temperature, velocity,
density, or composition profiles. Such profiles must be compared once complete
model inputs and a validated Fe prescription exist.

## Other local codes and source access

MOCASSIN 3.0 `source/update_mod.f90` reads `data/dielectronic.dat` into its
`diRec` array at all temperatures, then copies the complete `aldroPequi`
high-temperature array into `diRec` when `TeUsed > 60000 K`. The local
`dielectronic.dat` has no Fe (`Z=26`) row, while `aldrovandi.dat` has Fe
entries for 22, 21, and 20 electrons. The caller initializes `diRec` to zero
before the call. At lower temperatures, `calcAlpha` can fill a missing Fe DR
entry from the corresponding O ion: in its `lgBadnell` branch this requires O
to be enabled, whereas the other branch applies the substitution whenever
`ion < 7`. This can affect Fe V and VI (`ion=5,6`), but excludes Fe VII
(`ion=7`). The final rate then adds `max(0, diRec)`. Thus a nonzero
low-temperature Fe V or VI contribution from this path can be an O-based
estimate; it is not an independent Fe DR calculation. The MAPPINGS V
source was searched for Fe recombination data, but its rate provenance and
applicable Fe stages were not established in this audit; it is not a verified
replacement.

NASA ADS API lookup was attempted first, as required for literature checks,
but DNS resolution of `api.adsabs.harvard.edu` failed in this environment.
Publisher final versions were therefore not retrieved. The cited papers'
specific Fe coefficients, methods, and validity intervals remain unverified.
The Cloudy private-communication fit has no accessible final paper indicated
by its source label. Source review is required before promoting any candidate.

## Reproduction and integration boundary

Run `python tools/fitting/fe_recombination_expansion.py` for the coefficient
table, `MPLCONFIGDIR=/tmp/mochii-fe-recombination-mpl XDG_CACHE_HOME=/tmp/mochii-fe-recombination-cache python tools/fitting/plot_fe_recombination_expansion.py`
for the full curves, and
`python -m unittest tests/fitting/test_fe_recombination_expansion.py` for
source-index and expression checks. The scripts read external files and do not
edit the production generator or atomic products. A future
production choice must establish published low-temperature DR for all three
transitions, verify its validity range and units, and then remove the explicit
Fe-stage guard in `tools/fitting/make_element_data.py`. That work also depends
on the separate Fe ionization thresholds and charge-exchange boundary review.
