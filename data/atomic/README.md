# data/atomic

Atomic rate and excitation files used by the production calculation. Cooling
fits produced by `tools/fitting` identify their source, fit form, fit interval,
and measured approximation error. Source-node atoms instead preserve their
collision tables, interpolation choice, numerical limits, and source-accuracy
qualifications. An exponential fit is not required for absolute atom cooling.

## Contents (state of 2026-10-11)

| path | content | read by |
|---|---|---|
| `element_<el>.txt` (11: c, n, o, ne, s, ar, mg, fe, si, cl, ca) | stages I-VII of each element (77 stages, 66 transitions): ionization potentials, phfit2 photoionization rows (and the Ne III `PHTABLE`), collisional ionization, recombination, charge exchange with hydrogen and its heat (see Element files) | `species_mod` |
| `ground_term_weights.txt` | ground-term statistical weight and first excited term of every metal stage (NIST ground terms; `tools/fitting/make_ground_term_weights.py`) | `milne_recomb_spectrum_mod` (diffuse channel 6) |
| `badnell_rr.dat`, `badnell_dr.dat` | Badnell radiative and dielectronic recombination fit parameters (2023 release) | the recombination generators |
| `nlevel_h_1.txt`, `cooling_h_1.txt` | the H I atom and its cooling fit | `cooling_mod`, `lines_mod` |
| `sh95_{hi,heii}_case{A,B}.txt`, `hei_porter_caseB.txt` | recombination-line emissivity tables | `lines_mod` |
| `he1_case_b.dat` | the Porter He I grid from Cloudy, input of `tools/fitting/make_hei_lines.py`; not read at run time | -- |
| `h_2s_2p_mixing.txt` | H(2s) -> H(2p) collisional rates by collider and A_2s | two-photon quenching |
| `hei_metastable.txt` | He I 2^3S rate coefficients and the 10833 A line constants | `hei_metastable_mod` |
| `excitation/chianti/` | the default metal atoms (below) | `nlevel_mod`, `nlevel_cooling_mod` |
| `excitation/mocassin_v2/` | the MOCASSIN v2 atoms (below) | the same, with `metal_excitation_model = 'mocassin_v2'` |
| `norad/ne3.px.gd.txt` | Nahar (2019) Ne III ground-term photoionization cross sections (NORAD), source of the Ne III `PHTABLE` | the generators only |
| `source/si_i/` | the NIST MCHF Si I file the Si I strengths come from (`tools/fitting/silicon_source_records.py`) | the generators only |

## Excitation models

`par%metal_excitation_model` selects `excitation/<model>/` for the level
energies, radiative probabilities, electron collisions, matching cooling fits
and hydrogen-collision files of every metal stage
(`nlevel_<el>_<stage>.txt`, `cooling_<el>_<stage>.txt`, `hcoll_<el>_<stage>.txt`,
and `manifest.json` with the SHA-256 of every file). Every stage with an atom
in the selected model cools and emits; there is no switch that drops some.

- `'chianti'` (default): full CHIANTI 11.0.2 atoms (every bound,
  collision-connected level, every radiative and collision row) for 67
  stages, Stout atoms for Cl VII, Ar II and Ca IV: 70 atoms (6629 levels;
  Fe III 322, Si VII 460), 70 cooling fits, hydrogen collisions for C I,
  C II, O I. No line cooling for C VII (bare), Ar I, Cl I, Ca I, Fe I (no
  data) and Cl V, Ca III (CHIANTI has only excitation-autoionization data).
  Report: `md/atomic_stage_expansion/chianti_model_2026-10-11/README.md`.
- `'mocassin_v2'`: the MOCASSIN v2 atoms copied exactly (55 atoms, 571
  levels, 49 cooling fits and their listing `cooling_fit_parameters.txt`), plus the 141 CHIANTI forbidden A-values of Fe II
  levels 14-35, for parity runs only.

The cooling fits of both models target the stationary balance at
n_e = 1e-6 cm^-3, the reference density the (T, n_e) tables are divided by;
four fits were refit to it on 2026-10-11 (chianti `ar_7`, `s_5`, `o_4`,
mocassin_v2 `n_2`), and the others differ from the strict n_e -> 0 limit by
less than 1e-4, so their headers still name that limit. The chianti fits
carry 2-21 terms (mocassin_v2 up to 69); every chianti header states a
maximum fit error at or below 0.98% over 1e3-1e5 K (N V; read from the
headers, not re-measured here).

Beyond a collision table's temperature range Upsilon is held at the end
node (TEMPERATURE_TABLE rows); Burgess-Tully rows hold the scaled variable
and descale and form the Maxwellian at the true T. `UPSILON_CONSTANT_UPPER_K`
records are read and ignored. The `'extended'` model, `par%fe_levels_full`
and `par%additional_line_cooling`, the root metal `nlevel_*`/`cooling_*`/
`hcoll_*` files and `cooling_fit_parameters.txt` were removed on 2026-10-11
(`md/ATOMIC_PHYSICS_RESTORATION_RESULTS_2026-10-11.md`); an input naming
them is refused with the replacement.

## Collisional line-cooling fits

The cooling used in the thermal balance is stored as closed-form fits, one
file for each ion in the selected model
(`excitation/<model>/cooling_<element>_<stage>.txt`, e.g.
`excitation/chianti/cooling_o_3.txt` for O III), and `cooling_h_1.txt` for
H I at this directory's root. (`cooling_fit_parameters.txt`, a single
listing of the former root fits written by
`tools/fitting/make_cooling_table.py`, was removed with them.)

The fitted quantity is the rate coefficient

```
Lambda_ion(T) = T^(-1/2) * sum_{i=1}^{N} A_i * exp(-T_i / T)
```

| symbol | meaning | unit |
|---|---|---|
| `Lambda_ion` | line-cooling rate coefficient of the ion | erg cm^3 s^-1 |
| `A_i` | amplitude of term *i* | erg cm^3 s^-1 K^(1/2) |
| `T_i` | excitation temperature of term *i* | K |
| `T` | electron temperature | K |
| `N` | number of terms (2 to 21 in `chianti`, up to 69 in `mocassin_v2`; one line for each) | — |

The volumetric cooling rate of the ion is `n_e * n_ion * Lambda_ion(T)`
in erg cm^-3 s^-1; MoCHII sums it over the ions of the active registry.

File layout of `cooling_<ion>.txt`: comment header (lines starting with
`#`), then `N`, then `N` lines of `A_i  T_i`.

The `chianti` fits use CHIANTI v11.0.2 level energies and effective
collision strengths with Burgess-Tully descaling over 10^3-10^5 K (the three
Stout atoms their own tables). The stationary balance at the reference
density n_e = 1e-6 cm^-3 includes collisions and spontaneous radiative
cascades; it equals the strict n_e -> 0 limit to 1e-4 for all but four
fits, which were refit to it (Excitation models, above).
When only the ground state is radiatively stable, its excitation power
is the total cooling. C III and Si III contain another radiatively stable
state in the source atom; their stationary populations must also be solved.
Their cooling fits were corrected on October 9, 2026. Older results using
these fits are stale until recalculated. Previous products are preserved.
(That sentence dates from the root fits of October 9; the `chianti` C III
and Si III fits were generated from full atoms on October 11.)

Observed level energies are used when available, otherwise theoretical
energies are used. The same gaps enter excitation Boltzmann factors and
cascade energy accounting. Exponential coefficients describe a fit to the
total cooling and can be signed; they are not individual emission rates.
The total must be positive and accurate over the stated range;
`cooling_eval` holds a signed fit at zero where it would go negative.

The fits omit finite-density collisional de-excitation. At run time,
`par%cooling_model = 'local_ne'` (the default) applies a density correction
from a (T, n_e) table built from `nlevel_<el>_<stage>.txt` of the same
model. That selected atom also
supplies line diagnostics. `par%cooling_model = 'low_density'` uses the
fits directly, within their temperature and low-density validity limits.
Reduced-atom completeness and its low-density reference require separate
validation; matching filenames does not establish identical atomic models.

The fits hold over their stated range, 10^3-10^5 K, and must not be
extrapolated below it: their exponentials do not follow the true excitation
temperatures of the far-infrared fine-structure lines ([C I] 609.1 um is
23.6 K; the `chianti` C I fit has one term with T_i = 8.8e-4 K and the next
at 4123 K, the root fit it replaced had its lowest at 2054 K), so in cold
neutral gas the electron-impact line cooling from these fits is not the
atom's.  That gas is also where hydrogen, not electrons, is the collider: the
C I, C II and O I atoms carry atomic-hydrogen impact rates for their ground
terms (`excitation/<model>/hcoll_<el>_<stage>.txt`, below), solved in the same statistical equilibrium
as the electron rates.  The thermal balance adds, for those ions, the
increment of the radiated power over the electron-only populations
(`species_mod::metal_cooling_H`), and the line output uses the same atoms.

## Atomic-hydrogen impact on fine-structure lines

`hcoll_<element>_<stage>.txt` sits beside `nlevel_<element>_<stage>.txt` in
`excitation/<model>/` (both models carry them for C I, C II and O I) and is read
with it by `nlevel_mod::nlevel_load` whenever it exists; its rates enter the
level equilibrium, the cooling and the lines under
`par%neutral_hydrogen_excitation` (default on; off only for parity with
electron-only codes).  Generated by
`tools/fitting/make_hcoll_fine_structure.py`.

| file | transitions | source |
|---|---|---|
| `hcoll_c_1.txt` | C I 3P0-3P1 (609 um), 3P0-3P2, 3P1-3P2 (370 um) | Abrahamsson, Krems & Dalgarno (2007), ApJ 654, 1171, Table 2 (5-1000 K) |
| `hcoll_c_2.txt` | C II 2P1/2-2P3/2 (158 um) | Barinovs et al. (2005), ApJ 620, 537: de-excitation 7.6e-10 (T/100)^0.14 cm^3 s^-1 (20-2000 K) |
| `hcoll_o_1.txt` | O I 3P2-3P1 (63 um), 3P2-3P0, 3P1-3P0 (145 um) | Abrahamsson et al. (2007), Table 4 (30-1000 K) |

File layout: comment header, then `NHCOL n`, then `n` rows

```
lo hi form Tmin Tmax c_1 ... c_10
```

with `lo`, `hi` the level indices of the n-level file (`lo` the lower level).

| form | rate | coefficients |
|---|---|---|
| 1 | Abrahamsson et al. (2007) eq. (1), EXCITATION lo -> hi: `ln(k/1e-11 cm^3 s^-1) = sum_{i=0}^{8} c_(i+2) T^(-i c_1)` | `c_1` = b, `c_2 ...` = a_i |
| 2 | DE-excitation hi -> lo: `q = c_1 (T/100)^c_2` cm^3 s^-1 | `c_1`, `c_2` |

The de-excitation rate is evaluated at T clamped to `[Tmin, Tmax]` (the range
of the calculation; no extrapolation), and excitation follows by detailed
balance at the actual T with the level energies and statistical weights of
the n-level file.  The negative signs that the Abrahamsson PDF text layer
drops were restored from the printed glyphs; the generator checks the fits
against the paper's Tables 1 and 3 to <= 1.0% (the one exception, O I k1 at
50 K, +9.1%, is the paper's own fit at the edge of its range).

## Element files

`element_<element>.txt` holds, one `TRANSITION it` block for each ionization
stage `it -> it+1`, the photoionization, collisional-ionization,
recombination and charge-exchange data (`tools/fitting/make_element_data.py`;
the header of each file documents every row).  The photoionization rows
follow Verner's `phfit2.f`, parsed from it directly:

| row | content |
|---|---|
| `ETH` | ionization potential [eV] (NIST); the cross section is zero below it |
| `PHOTO E_th E0 sigma0 ya P yw y0 y1` | VFKY96 outer-shell fit (phfit2 `PH2`), used below `E_inn`; present when `E_inn > 0` |
| `PHINNER E_inn` | phfit2 `einn`, the threshold of the first inner subshell [eV]; 0 for the elements absent from the VFKY96 table (Cl), 1e30 for ions with fewer than three electrons |
| `PHSHELL l E_th E0 sigma0 ya P yw` | Verner & Yakovlev (1995) fit of one occupied subshell (phfit2 `PH1`), one row for each subshell; at `E >= E_inn` the cross section is the sum over the subshells with `E >= E_th` |
| `PHTABLE n`, then `n` rows `eps sigma` | tabulated R-matrix cross section of the ground term [Mb] against the photoelectron energy `eps = E - ETH` [eV], linear in `eps`; replaces `PHOTO` below its last node `E_join = ETH + eps_n`, where it equals `PHOTO` (Ne III only: Nahar 2019 BPRM, `norad/ne3.px.gd.txt`; `tools/fitting/ground_term_photoionization_table.py`) |

`ETH`, the PHOTO cutoff, and the occupied-valence PHSHELL cutoff use the
identified NIST ASD snapshot dated 2022-11-15 in Cloudy c25.00. Inner PH1
edges and Voronov CI fit energies remain independent source parameters.
The terminal photon-band check refuses a band (`eion_max`) above the
ionization potential of the last included stage (VII) of an element present,
since the stage it would feed is not carried; with stages I-VII the lowest
such potential is Cl VII, 114.2 eV (Ar VII 124.4, Fe VII 125.0 eV), so the
default 13.6-100 eV band is accepted.

The photoheating coefficient is `sum_s sigma_s (1 - E_th,s/E)` over the
absorbing subshells (on a `PHTABLE`, `sigma (1 - ETH/E)`).  The reader refuses a
transition with no `PHSHELL` rows (a file written before they existed), with
`E_inn > 0` and no `PHOTO` row, or with a `PHTABLE` that is short of its count,
does not start at `eps = 0` or increase, carries a negative or nonfinite value,
reaches `E_inn`, or does not end on the `PHOTO` fit (relative 1e-8).

`norad/ne3.px.gd.txt` is the NORAD-Atomic-Data file of Nahar (2019, New Astron.
67, 97), downloaded unchanged on 2026-10-11 from
https://norad.astronomy.osu.edu/ne3/ne3.px.gd.txt (md5
964ecbb6669aef6f1ecabd8e61d7fd87): total photoionization cross sections, with
radiation damping, of Ne III 2s2 2p4 3P2, 3P1, 3P0.  It is the source of the
Ne III `PHTABLE` and of `tests/ne3_ground_term_photoionization/reference.txt`.
Inside the default 13.6-100 eV band the inner edges that matter are Mg I/II
2p (54.9/65.7 eV), Ca I/II 3p (34.4/40.9 eV), Fe I-III 3p (66.0-87.1 eV) and
the Cl I-V 3s subshells (25.3, 36.9, 50.2, 64.7, 80.0 eV; Cl VI's 3s is its
outer shell).

Charge exchange with hydrogen is `CXI`/`CXR form c_1..c_6` (forms 1-7 in the
file header).  Form 8 is a table: `CXR 8` followed by

```
CXRTAB n T_1 ... T_n k_1 ... k_n
```

(`n <= 1024`; T in K, k in cm^3 s^-1), interpolated log-log in T and held at
the end values outside the table.  N2+ + H and O2+ + H use it for the
Barragan et al. (2006, ApJ 636, 544) Table 1 k1 and Table 2 k4 rates
(26 points, 1e2-1e5 K).
With stages I-VII and all eleven elements, 45 reactions are carried
(7 ionization, 38 recombination; census in `tests/charge_exchange/`).
N4+ + H and Si2+ + H take the updated KF96 database distributed with
Cloudy (Ferland et al. 2013, RMxAA 49, 137; Cloudy c25.00
`data/ctrecombdata.dat`) in place of the printed KF96 Table 1 fits.
The recombination rows for q = 5, 6 (22 of them, `CXR 1`) are the Langevin
capture upper limit; they act only with `par%langevin_charge_exchange`
(default off).
The 13 existing KF96 higher-boundary fits now use individual published
temperature ranges, with at most 653 nodes and measured interpolation error
below 1e-5. Nitrogen's first-boundary fit uses 100--50000 K; its reverse
activation factor takes the actual temperature. Endpoint continuation outside
a published interval does not establish physical accuracy there.
Four formerly omitted Cl3/4 and Ca3/4 recombination reactions now use KF96
total rates with KF99 mean kinetic heat. Ca uses Table 1 fit ranges; Cl uses
1000--30000 K inside the Table 2 sampled span, with its Landau-Zener
limitations. Ground-state reactants and low-density radiative relaxation
are assumed. No excited-product reverse rate is inferred. See
[the rate audit](../../md/atomic_stage_expansion/charge_exchange_expansion.md)
and [the correction validation](../../md/atomic_stage_expansion/kf96_total_recombination_validation.md).

Form 9 resolves a transition by the product term of the ion (sulfur,
Zhao et al. 2005 Table III): `CXI 9` and `CXR 9`, then

```
CXSTATE nterm nT g_react dIP
CXSTATET T_1 ... T_nT
CXSTATEK label lev_lo lev_hi g E k_1 ... k_nT      (one row per term)
```

with `lev_lo..lev_hi` the term's levels in `nlevel_<el>_<stage+1>.txt`, `g`
its weight, `E` its energy [eV] and `k` the forward rate into it.  The
forward rate is the sum over terms; the reverse is each channel's detailed
balance, `k_j g_react / (2 g_j) exp[-(dIP - E_j)/kT]`, weighted by the
term's population from that n-level atom at the local (T, n_e).  The
first-stage reverse rates of Mg, Si and Fe are the detailed balance of the
total forward rate (the ion's excited terms in LTE), an upper bound on the
reverse from the ground term, and each file says so.  The He-H charge exchange of the H/He balance is not
in these files: its rates are in `src/ion_balance_mod.f90`, closed forms except
the non-radiative He+ + H rate (Loreau et al. 2014 by reciprocity), which is
tabulated into `src/heii_h_charge_exchange_data.inc` by
`tools/fitting/make_heii_h_charge_exchange.py`.

The corrected C III and Si III serialized cooling fits were measured on
801 temperatures from 1000 to 100000 K: maximum relative errors are
0.127901% and 0.438816%, respectively. Earlier quoted fit errors for other
ions are stale until their serialized fits are rechecked against the adopted
stationary source. Such checks must specify the atomic model and temperature
interval; small abundance does not establish physical validity. See
[the measured cooling validation](../../md/atomic_stage_expansion/cooling_validation.md).

## Wavelengths

Every wavelength written by MoCHII is a **vacuum** wavelength. The SH95 H I and
He II tables carry vacuum values directly, and the n-level wavelengths follow
from observed level energy differences, which are vacuum by construction. The
one table whose source uses air wavelengths is the Porter He I set (indexed by
air wavelength in the Cloudy file it is parsed from); `make_hei_lines.py`
converts those to vacuum on the way out with the inverse of the IAU standard
dispersion of air, so the emitted line table has one convention throughout. The
He I triplet line is therefore 10833 A here, not the 10830 A of the air-based
literature, and the identifiers that name it follow: the `LINE10833` keyword of
`hei_metastable.txt`, the `tau10833_map` reader, and the `HeI10833` label of
the Porter table.

## Other files

| file | content |
|---|---|
| `excitation/<model>/nlevel_<el>_<stage>.txt` | n-level atom: levels, A-values, effective collision strengths (`NLEVEL_VERSION 2` source nodes in both models; the Chebyshev form, read for `nlevel_h_1.txt` at the root) |
| `element_<el>.txt` | photoionization, recombination, collisional-ionization, and charge-exchange data of one element (see above) |
| `excitation/<model>/hcoll_<el>_<stage>.txt` | atomic-hydrogen impact rates for the ground-term fine structure (see above) |
| `badnell_rr.dat`, `badnell_dr.dat` | Badnell radiative and dielectronic recombination fit parameters |
| `sh95_*.txt`, `hei_porter_caseB.txt` | recombination-line emissivity tables |
| `he1_case_b.dat` | the Porter He I case-B grid copied from Cloudy c23.01, the input `tools/fitting/make_hei_lines.py` parses into `hei_porter_caseB.txt`; not read at run time |
| `h_2s_2p_mixing.txt` | H(2s) -> H(2p) collisional rate coefficients by collider (e, p, He+, He2+) and A_2s, for the two-photon quenching |
| `hei_metastable.txt` | He I 2^3S metastable rate coefficients and the 10833 A line constants; the `PENNING A c d1 d2 d3` row is the Garcia Munoz (2025, A&A 698, A199) Penning + associative ionization rate `Q_31 = A exp(c/T + d1 lnT + d2 (lnT)^2 + d3 (lnT)^3)` cm^3 s^-1 (`A` = 1e-9); the `UPSILON role lower n 2S+1 l g_lower dE_eV Upsilon(1..n)` rows are the Bray, Burgess, Fursa & Tully (2000, A&AS 146, 481) Tables 2 and 3 effective collision strengths, 1^1S -> n^3L and 2^3S -> n^(2S+1)L for n <= 5 (41 transitions, role `Q13`/`FEED`/`Q31A`/`Q31B`/`SINGLET`/`TRIPLET`), at the `UPSGRID` nodes log T = 3.75-5.75 as printed plus 3.00-3.50 from the CHIANTI he_1 temperature dependence anchored at 10^3.75 K; `q = BETA Upsilon exp(-dE/kT)/(g_l sqrt T)` with the exact Maxwellian constant `BETA` = 8.629132e-6 and observed term energies (`E23S` = E(2^3S)); the 2^3S -> 1^1S electron de-excitation is not a row of its own but the detailed-balance inverse of the `Q13` row; the recombination source is not in the file (alpha_3 = f_3S alpha_B(He II), hei_cascade_mod x recomb_mod; the Oklopcic & Hirata 2018 `ALPHA3` row is gone and a file carrying it is refused); `CION3 A p T0` is the Black (1981) electron-impact ionization of 2^3S, `A T^p exp(-T0/T)` |

Not in this directory: `data/gauntff_vh14.dat` (the van Hoof et al. 2014
free-free Gaunt factors, `par%gaunt_vh14`), `data/gamma{HI,HeI,HeII}.dat`
(the free-bound + free-free emission tables of the nebular continuum) and
`data/HeI2phot.dat` (the He I two-photon shape) sit one level up; the
hydrogenic recombination of H II and He III
(alpha_B, the 2s share f_2s in case A and case B, and the He II case-A
Lyman-line shares) is generated as Fortran DATA into
`src/hydrogenic_recombination_data.inc` by
`tools/fitting/make_hydrogenic_recombination.py`, from one zero-density
nl-resolved cascade on the Storey & Hummer (1991) hydrogenic routines, so
that the elemental rate functions of `recomb_mod` need no setup call.  It
replaced `heii_caseB_n2.txt` (SH95, four anchors, case B) on 2026-09-26.

## History of the October 9-11 preparation

Kept as a record. Stage counts ("43 metal stages", "production remains at ..."), "staged" or "prepared" data and root file locations in these subsections describe the tree of October 9-11, before the restoration of 2026-10-11; the present state is the sections above.

### Excitation corrections and source validity

Eleven production low-density cooling files were corrected after excluding
above-threshold states from ordinary bound excitation and preserving
nonnegative collision strengths between nonnegative source nodes. Rounded
source fits were checked at 801 temperatures from 1,000 to 100,000 K; their
maximum error is 0.9485%, within the stated 1% coefficient target. Linked
native evaluation differed from the rounded coefficients by at most
9.053e-12 relatively. Original coefficients are preserved in the validation
directories. These are corrections to existing stages, not activation of
the current 77-stage target; the earlier 64-stage plan is historical.

The earlier corrected five-level C I native atom passed 3,204 linked checks
with a maximum source cooling difference of 0.14966% and a rounded transition
fit error of 0.19306%. Those historical products are preserved. It has since
been replaced by the complete 42-level bound source after full-source and
registered-cooling verification; the earlier five-level check did not establish
source convergence.
S I's missing low-density cooling fit is now installed with 21 signed terms.
Its maximum native/source difference is 0.496965% over 801 temperatures; its
native evaluation matches its rounded primary coefficients within 2.818e-12.
Its five-level source-node atom preserves all sampled source-zero intervals.
S III now uses the full 53-level bound collision-connected source (355 radiative
rows and 1,378 collision rows), replacing the stale six-level fitted atom.
Their actual 401-by-161 registered tables passed grid/interior cooling and
application checks. Old S III and all failed products are preserved under
`md/atomic_stage_expansion/validation_2026-10-10/`. These changes do not activate
the expanded element roster or validate a spatial nebular calculation.

Ca IV now includes the published two-level ground fine-structure electron
pair from Pelan and Berrington (1995), Table 3, as preserved in Stout. The
801-point fit error is 0.261151%; its actual registered table has maximum
sampled interior/application differences of 0.167021%/0.163713%. Constant
first-node Upsilon below 10,000 K is the adopted approximation; the source
table ends at 1,000,000 K. Cooling from higher excited states remains
unestablished by that pair. The installed serialized coefficients agree with
the compiled reader within 3.800e-12.

Fe IV's missing cooling now uses the full 37-state collision-connected bound
source, with 391 radiative and 645 electron-collision rows. Its coronal limit
includes two radiatively stable states. The dense-grid fit error is
0.258290%, compiled coefficient reading agrees within 1.539e-12, and the
actual registered table passes its direct-solve checks. The aggregate listing
now contains 31 ions and 182 terms, preserving every primary coefficient
exactly. These additions retain the 43-stage production roster.

Source-node files begin with `NLEVEL_VERSION 2` and retain original collision
nodes and interpolation choices. Existing Chebyshev files remain readable.
`DENSITY_RATIO 1` selects fitted cooling multiplied by the calculated density
ratio for a complete or independently source-converged bound atom. It prevents
fit or interpolation errors from creating a fictitious unsuppressed component.
The absence of that marker retains the older partial-atom approximation.
Source convergence for reduced atoms is currently established only between
1,000 and 100,000 K; the larger numerical table does not establish physical
accuracy outside that interval. Fe II's negative high-temperature source
endpoint also remains outside the tested region.
See [interpolation correction](../../md/atomic_stage_expansion/excitation_interpolation_correction.md)
and [storage design](../../md/atomic_stage_expansion/excitation_storage_design.md).

### Reviewed Fe recombination and staged N IV

The existing Fe IV to Fe III recombination records in `element_fe.txt`
use Zhang et al. (2026), ApJS 285:50, final Tables 6/7: ground-level RR
and zero-density S1 DR. Their domain is 180--9e7 K; extrapolation outside
it is uncertified. These fits omit density suppression and capture from
metastable parents. The compact DR fit differs from author ADF09 data by
up to 4.718% at inspected nodes within the domain, exceeding the paper's
stated 2% fit accuracy. Prior results using the old Fe IV low-temperature
recombination records are stale until recalculated for reuse.

The generator also implements ground-level q=4--6 Fe rates, each valid
over 20*q^2--1e7*q^2 K, with the same stated limitations. Fe V can be
generated; higher Fe stages still require selected charge-transfer inputs.
Production remains at 43 metal stages. Numerical evidence is in
[`final_fe_recombination_validation.md`](../../md/atomic_stage_expansion/final_fe_recombination_validation.md).

N IV cooling and 29-level data were generated and checked in staging only.
They use the checked Cloudy Stout N IV energies, observed for levels 1--29
and BSR above, with the original CHIANTI radiative and collision indices.
Their numerical fit interval 1,000--100,000 K does not establish source
resonance accuracy or convergence of the complete target-state expansion.
See [`niv_staging_validation.md`](../../md/atomic_stage_expansion/niv_staging_validation.md).

### Adopted added high-charge hydrogen transfer

The selected element generator uses the user-adopted reduced-mass Langevin
capture upper limit for N q=5; S, Ar, Fe and Cl q=5--6; and Ca q=5.
Classical unit electron-transfer probability is assumed for H(1s); no physical
temperature-accuracy interval is established. Mean kinetic release is the
9.5-bohr avoided-crossing estimate 2.86(q-1) eV using the exit charges +1
and +(q-1), explicitly correcting the inconsistent literal 2.86q prose in
KF99. Product branching, reverse rates and diffuse cascade reabsorption remain
unresolved. These data are staged; production still declares 43 metal stages.
See [adoption validation](../../md/atomic_stage_expansion/langevin_adoption_validation.md).

### Full bound C I, C II, and N I excitation atoms

The source-node atoms now use 42, 53, and 50 levels, respectively, after
independent full-source population/cooling and actual registered-table
validation. `DENSITY_RATIO 1` applies their computed density dependence to
the fitted low-density cooling. The previous five-level atoms could not
resolve that density dependence; their retained low-density fit residual
produced artificial high-density cooling. The low-density fits and existing
C I/C II H(1s) collision tables were preserved. Their original ground-term
level identities remain unchanged, and actual published-rate, detailed-balance,
and thermal/emission checks passed with the full atoms. Source convergence
is measured over 1,000--100,000 K; it does not establish accuracy of every
source outside that interval. See
[physical comparison and validation](../../md/atomic_stage_expansion/existing_compact_atom_correction.md).

N II and S II now use the full identified bound sources with 58 and 70
levels and density-ratio cooling after independent native and registered
checks. Original five-level atoms are preserved in the same installation
record. Their low-density fits were retained. O II and O III now use full
196- and 177-level atoms with source-node collisions. Higher E1 probabilities
preserve line strengths at the selected physical photon energies and emitting
statistical weights; unidentified same-parity rates remain explicit effective
rate approximations. Native independent balance and energy closure passed.
Full registered tables passed with maximum application errors of 0.197768%
and 0.143563%. These numerical errors do not establish atomic accuracy. See
[oxygen source and validation](../../md/atomic_stage_expansion/oxygen_valence_convergence.md).

### C IV and Si V source-temperature cooling

C IV and Si V now have 24-state and 169-state source-node atoms and
low-density cooling fits. Below their source minima of 3,200 and 5,000 K,
Upsilon is held constant while Maxwellian factors retain actual temperature.
The fit range is 1,000--100,000 K; rounded numerical errors are 0.344146%
and 0.929616%, not physical uncertainty estimates. Installation and native
validation are recorded in
[the cooling products](../../md/atomic_stage_expansion/validation_2026-10-10/bt_floor_cooling_products/README.md).

The O III optical transitions 4 to 2 and 4 to 3 use the final Storey and
Zeippen (2000), MNRAS 312, 813, Table 3 totals 6.791e-3 and 2.046e-2 s^-1,
matching the directly inspected MOCASSIN v2 rows. The full O III atom also
uses the E1 energy treatment described above. Its low-density cooling fit was preserved because the
coronal emitted-energy sum is unchanged within numerical precision.

### Ar II ground fine-structure cooling

Ar II now uses the published Pelan and Berrington (1995), A&AS 110, 209,
Table 3 two-state Upsilon source. The original 11 nodes are preserved,
with constant first-node Upsilon below 1,580 K. The registered cooling
temperature grid ends at 158,000 K and rejects higher temperatures.
Its low-density fit covers 1,000--100,000 K with numerical error 0.202152%.
The paper excludes Ar II from its better-than-10-percent accuracy estimate.
Higher Ar II excitation is not included. See
[installation and validation](../../md/atomic_stage_expansion/validation_2026-10-10/source_temperature_registry/README.md).

### Full Mg III excitation and source-energy correction

Mg III now uses 167 bound states, 3,364 radiative rows, and 13,861 original
collision pairs with 207,915 nodes. NIST shows that the observed energy
assigned to source ID167 belongs to a different J=1 7d state. The J=0 7p
state therefore uses its own theoretical 634430 cm^-1 energy; external
CHIANTI bytes remain unchanged. Higher E1 rates preserve source strengths
with physical frequency cubed and emitting statistical weight. The first27
A values remain unchanged, and 112 unidentified multipole rates are explicit
effective-rate approximations. Upsilon is held at 1,800 K below that source
minimum while Maxwellian factors use actual temperature.

The 1,000--100,000 K fit has maximum numerical error 0.210705%. Full native
profiles, independent collision/radiation energy closure, and the registered
table passed. Maximum checked table application error is 0.000534832%.
These are numerical representation errors, not atomic accuracy estimates.
The previous 27-state candidate loses up to 17.1870% of full-source cooling
and is not installed. See
[source identity, approximations, and validation](../../md/atomic_stage_expansion/mgiii_radiative_consistency.md).


### Ne I spontaneous-channel completion

Ne I now has a complete 31-state source atom with 144 spontaneous rows,
465 collision pairs, and 6,975 original CHIANTI nodes. NIST-dated Stout
rates A(3,15)=9.00e5 and A(5,15)=6.82e7 s^-1 restore two omitted decays;
all 142 original CHIANTI A values remain unchanged. The ground is the only
radiatively closed state. Upsilon is held constant below an inferred common
finite-node floor of 1000.36 K while actual T enters Maxwellian factors.
This floor and the 1,000--100,000 K adopted numerical interval do not certify
published source accuracy. The CHIANTI conversion Methods remain unverified.

Source/native cooling and complete populations passed 401-temperature,
four-density checks; independent energy balances passed. The installed
fit's 801-temperature maximum error is 0.0317686%; actual GNU coefficient
reader error is 7.78e-12. The 64,561-cell Intel density table's sampled
application error is 0.0041481%. Installation preserves exact validated
bytes, and a production debug reader opened the installed atom successfully.
See [the root review](../../md/atomic_stage_expansion/validation_2026-10-10/ne_i_radiative_completion/main_review.json).

### Prepared O IV full-bound cooling data

`nlevel_o_4.txt` and `cooling_o_4.txt` contain the adopted 147-level source-data excitation model and 20-term low-density fit. The full atom supplies finite-density suppression through its density ratio. Imported spontaneous rates remain atomic-model estimates; detailed A replacement provenance and absolute accuracy are unverified. Constant Upsilon below the published 3,200 K collision boundary is the adopted cold approximation; rates retain actual gas temperature. Numerical use is 1,000--100,000 K. The complete model and installed reader checks are recorded in `md/atomic_stage_expansion/validation_2026-10-10/o_iv_full_completion/installation/README.md`.

Oxygen's production ELEMENT record still declares three stages. These prepared O IV data do not independently add a fourth oxygen stage or activate the current 77-metal-stage target.

`nlevel_o_5.txt` and `cooling_o_5.txt` now also provide the qualified 166-level O V source-data model and 19-term fit. Constant Upsilon below the distributed 5,000 K table boundary is the adopted approximation, with actual gas temperature in Maxwellian rates. The table boundary differs from the general interval in the final 2014 paper and does not certify source accuracy. Imported A values are estimates; no unused gf/wavelength field enters cooling. See `md/atomic_stage_expansion/validation_2026-10-10/o_v_full_completion/installation/README.md`. The preserved root listing now contains 39 species including H I and 381 terms; runtime metal listings are model-specific as described above. Production oxygen remains at three stages.


### Supplied MOCASSIN models with absolute cooling

23 previously absent finite-level atoms were installed with the validated
MOCASSIN v2 source-table conversion. Existing atoms and rate files were not
replaced. The source headers retain atomic-accuracy qualifications, adopted
cold constant-strength behavior, and positive log-log PCHIP interpolation.
These atoms can calculate absolute cooling without an exponential coefficient
file. Cl V is already within the active chlorine chain; its actual consumer
passed 264 samples in both existing density modes. The other new atoms await
expanded stage activation. Mg I now uses the approved constant-strength
model in both coherent selections; only Ca VII remains unresolved. See
[installation and source limits](../../md/atomic_stage_expansion/validation_2026-10-11/missing_atom_installation/README.md).
