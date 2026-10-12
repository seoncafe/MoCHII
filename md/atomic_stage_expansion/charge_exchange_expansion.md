# Charge-transfer review for added metal boundaries

The nine rates below are published recombination fits for reactions
`X^q+ + H -> X^(q-1)+ + H+`. They have been checked against Table 1 of the
final journal version of Kingdon and Ferland (1996, *ApJS* 106, 205; PDF in
`md/heii_h_charge_exchange/references/`). The code reader uses the first four
rows of each Cloudy element block, indexed by reactant charge `q=1..4`.
`charge_exchange_expansion.py` now supplies the nine published KF96 total
rates and separately verified KF99 Table 1 mean kinetic heat for selective
staging. The generator writes these boundaries with `CXHEAT` and `CXRTAB`.
Production remains at 43 stages until the remaining physics is complete.
KF96 Table 1 lists rate coefficients, temperature
ranges, and references in its seven columns. It does **not** print an energy
defect or product-state column. Thus the Cloudy energy defects cannot be
claimed to be values printed in KF96 Table 1. Cloudy's source directly uses
the stored defect as gas heat per recombination; it does not identify the
populated product levels. The original final
[KF99 publisher PDF](references/Kingdon_Ferland_1999_ApJ_516_L107.pdf) was
obtained; Section 2 and Table 1 were read and visually checked. These published
means assume ground-state reactants and low-density radiative relaxation, and
omit diffuse cascade reabsorption and product collisional redistribution.
Individual product channels remain unverified. Zero heat or ground products
would not represent this published mean-energy approximation.

| Added boundary | Cloudy row within element (zero based) | KF96/Cloudy rate at 10,000 K (cm^3 s^-1) | Cloudy gas heat per event (eV) | Fit range (K) |
|---|---:|---:|---:|---:|
| C 4 | 3 | 2.18455e-9 | +11.30 | 10–100,000 |
| N 3 | 2 | 3.32815e-9 | +6.40 | 10–100,000 |
| O 3 | 2 | 4.14226e-9 | +5.00 | 1,000–50,000 |
| O 4 | 3 | 2.60181e-10 | +8.47 | 1,000–30,000 |
| Ne 3 | 2 | 5.65491e-9 | +5.82 | 5,000–50,000 |
| Ne 4 | 3 | 6.59560e-9 | +8.60 | 1,000–30,000 |
| S 4 | 3 | 6.49855e-9 | +8.60 | 1,000–30,000 |
| Ar 4 | 3 | 6.49938e-9 | +8.60 | 1,000–30,000 |
| Mg 3 | 2 | 6.49889e-9 | +5.73 | 1,000–30,000 |

The rates are direct evaluations of `a * 1e-9 * T4^b *
(1 + c*exp(d*T4))`, with `T4=T/10,000 K`, using the parsed Cloudy rows.
The displayed `+` sign means heat added to gas by the reaction as written,
following Cloudy's `SumCTHeat_v` at `atmdat_char_tran.cpp:595–608`. In the
MoCHII convention it would be `dE_rec`, the numeric value in a `CXHEAT` row.
Cloudy multiplies its stored column 7 directly by the recombination rate and
the reactant densities; it is not a threshold or a product excitation energy.
The defects `5.73` and `8.60` eV recur across unlike ions and nearly equal
`2.86*(q-1)` eV for charges three and four. Cloudy's explicit high-charge
heating rule is `2.86*ion` eV. KF99 adopts these repeated values through a
mean avoided-crossing approximation. They are not product-level energies.
The exit Coulomb potential involves charges `+1` and `+(q-1)`, giving
`27.2114*(q-1)/9.5` eV at 9.5 bohr. KF99's higher-charge prose uses `2.86q`,
inconsistent with that definition and its repeated table values. Cloudy's
product-charge indexing is not changed on that basis.
Temperatures outside each table's stated interval are clamped in the audit
function. This reproduces Cloudy's numerical convention; it does not establish
physical accuracy outside the fitted interval. KF96 specifically cautions
against distant extrapolation. Some rows, including O 4, are flagged as
uncertain in the published table. Rate and mean-heat records now support
staging; they do not establish complete cooling or reverse channels.

N4 now uses the original KF96 fit for selective staging; the source comparison remains important: KF96 Table 1 and both local MOCASSIN versions list
`(3.28, 0.52, -0.52, -0.19)`, whereas Cloudy c25.00 lists
`(2.95, 0.55, -0.39, -1.07)` with a 1,000–1,000,000 K interval. Direct
evaluation at 10,000 K gives 1.870e-9 and 2.555e-9 cm^3 s^-1, respectively.
The provenance of Cloudy's changed fit is not established. Fe4 has a Cloudy
fit but no matching KF96 Table 1 row. KF99 gives N4 and Fe4 mean heat
(11.00 and 10.00 eV). The follow-up identified Fe4's KF96 Table 2
Landau-Zener total, including its radiative-transfer caveat; Fe stage
ground-level recombination has now passed final-source and native checks. See
[further validation](kf96_total_recombination_validation.md).
Both now belong to the reviewed staging set. The other added boundaries are N 5; S 5–6; Ar 5–6; Fe 5–6;
Cl 5–6; and Ca 5. Cloudy `HCTRecom` uses a configurable default
`1.92e-9 * q cm^3 s^-1` when `q >= 5`, regardless of species or temperature.
This is a capture estimate and does not establish the electron-transfer
probability, product channels, or heat. The user adopted the reduced-mass unit-transfer capture upper limit for the
ten q >= 5 target boundaries. The separate mean heat uses the physically
consistent 9.5-bohr exit-channel estimate 2.86(q-1) eV, matching Cloudy
indexing while explicitly identifying the inconsistent literal KF99 prose.
See [adoption validation](langevin_adoption_validation.md). Local MOCASSIN
tables stop at charge four for most elements, and the inspected MAPPINGS
heating routine consumes tabulated energy defects without identifying these
missing product channels.

KF96 Section 2 explicitly allows excited recombined products, while its
Table 1 prints total rate fits rather than resolved product populations. For
state-resolved refinements beyond the adopted KF99 mean heat, determine significant
product channels and branching, compare their excitation energies with the
physical ionization-potential difference, and decide whether photons escape
or reprocess locally. Rate availability and heat availability must remain
separate decisions. The PH1 fit edge must not silently replace a physical
ionization potential in this calculation.

## Source record and checks

| Source | SHA-256 | Direct use |
|---|---|---|
| `md/heii_h_charge_exchange/references/Kingdon_Ferland_1996_ApJS_106_205.pdf` | `9de90ffa22b42fe59aea23b42c4fd2daa37a5efb0ff89d80355d5e11a227cf26` | Final journal Table 1, Section 2, and fit caveats |
| `references/Kingdon_Ferland_1999_ApJ_516_L107.pdf` | `2c681fe32532baddee5d19a67e472cefc51f6a0135de01ab82eb9aa182e78def` | Final publisher Table 1 and mean heat assumptions |
| `/home/kiseon/CLOUDY/c25.00/data/ctrecombdata.dat` | `4ff7f4aeaf24dd5e0c6ab450694e480ef16136fab443c4c68b7acd799c5cbe45` | Coefficients, fit intervals, and stored defects |
| `/home/kiseon/CLOUDY/c25.00/source/atmdat_char_tran.cpp` | `22a1443eed27765841679a705f53120bb25fbbf8d35757942e2096c73ce605c9` | `HCTRecom` indexing, clamping, and charge-five default |
| `/home/kiseon/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/source/update_mod.f90` | `54072877e18dd5d4b9cb386c95ebd24dcad4240c094eeb7d8733113e9a737330` | Independent table transcription and N 4 discrepancy |
| `/home/kiseon/RT_Codes/MOCASSIN/MOCASSIN-3.0/source/update_mod.f90` | `4c7fbbafe43a38836cc75e7a272071c760c12ce72cf9d636ef19f4a4908662a5` | Same N 4 coefficients; Ar 4 exponent differs from KF96 |

The MOCASSIN versions also put Ca-labelled rows into `chex(18,...)`, which
overwrites Ar rows. This is a confirmed defect in preserved external source;
its Ca entries are not used here. KF96 and Cloudy agree on Ar 4's `b=0.85`,
whereas MOCASSIN lists `b=2.12`. The ADS API was attempted first for a
citation check, but DNS resolution was unavailable. The local final journal
PDF was read directly with `pdftotext -layout`; no abstract or preprint was
substituted. The original rate audit did not regenerate atomic products.
The subsequent [KF99 mean heat validation](kf99_mean_heat_validation.md)
generated separate staging files and left production files unchanged.

## Existing KF96 rate intervals

The original `CX_HIGH` form-6 evaluator clamped all rows to 5,000–50,000 K.
KF96 Table 1 gives the following reaction-specific intervals. The published
fits and Cloudy coefficients match for twelve of these thirteen rows. Silicon
charge two is a documented exception; its published interval also differs
from Cloudy's revised interval. The rate reader returns the published fit and
interval for all thirteen rows without changing production data.

| Existing boundary | KF96 Table 1 interval (K) | Cloudy interval (K) | Coefficients match? |
|---|---:|---:|---|
| C 2 | 5,000–50,000 | 5,000–50,000 | Yes |
| C 3 | 1,000–100,000 | 1,000–100,000 | Yes |
| Ne 2 | 5,000–50,000 | 5,000–50,000 | Yes |
| S 2, S 3 | 1,000–30,000 | 1,000–30,000 | Yes |
| Ar 2, Ar 3 | 1,000–30,000 | 1,000–30,000 | Yes |
| Mg 2 | 1,000–30,000 | 1,000–30,000 | Yes |
| Fe 2, Fe 3 | 1,000–100,000 | 1,000–100,000 | Yes |
| Si 2 | 10–1,000,000 | 500–100,000 | No |
| Si 3 | 1,000–30,000 | 1,000–30,000 | Yes |
| Si 4 | 1,000–50,000 | 1,000–50,000 | Yes |

For Si 2, KF96 Table 1 and current `CX_HIGH` have
`(a,b,c,d)=(1.23,0.24,3.17,0.00418)`, while Cloudy has
`(6.77,0.0736,-0.43,-0.11)`. The origin of that revision has not been
confirmed from a primary publication. Direct evaluation at 10,000 K gives
`5.14543e-9` and `4.16214e-9 cm^3 s^-1`, respectively. N 2 and O 2 use the later Barragan
rates and are excluded from this KF96 comparison. The published fit limits
do not establish accuracy outside those intervals; clamping is a numerical
choice that needs separate physical review.

`tools/fitting/kf96_rate_table.py` prepares adaptive log-temperature nodes
for the existing form-8 `CXRTAB` log-log interpolation. The element generator
now writes these tables for the thirteen existing KF96 higher-boundary rows
after comparing the original `CX_HIGH` coefficients with the published fit.
It refines each
segment until its quarter, middle, and three-quarter checks differ from the
published fit by at most `5e-6`. An independent 2,001-point check over each
published interval, after six-decimal scientific output rounding, found the
largest relative error `5.629e-6` for Fe 3; the maximum table length was 653
nodes for C 2. The distinct private Fortran `MAX_CXR_TAB=1024` can load
these tables, while `MAX_CX_TAB=40` remains for heating tables. Beyond each published interval, form 8 holds
the endpoint rate. This is a finite continuation, not a claim that the KF96
fit is physically valid at lower or higher temperatures.

KF96 Table 1 N+ + H has `(1.01e-3,-0.29,-0.92,-8.38)` in units of
`1e-9 cm^3 s^-1` over 100–50,000 K; its Table 3 reverse N + H+ has
`(4.55e-3,-0.29,-0.92,-8.38)` over the same range with
`Delta E/k=1.086e4 K`. These match the current first-transition `CXR` and
`CXI` coefficients, respectively. Table 3 Mg + H+ has
`(9.76e-3,3.14,55.54,-1.12)` over 5,000–50,000 K, matching current `CXI`.
The first-transition Mg+ + H form-7 fit and Fe first-transition rates have
no matching KF96 Table 1 or Table 3 row; their validity must be checked
against their own cited sources. KF96 Table 3's Fe++ row concerns the second
charge boundary and does not validate Fe + H+.

After conversion to `unittest.TestCase`, direct execution of
`python3 -m unittest tests/fitting/test_charge_exchange_expansion.py -v`
ran three focused tests successfully. `python3 tests/fitting/run_checks.py
--result /tmp/kf96_full_fitting_checks.json` included them in its discovery
and passed 25 tests overall after the separate adaptive-table module's three
tests were added. The third test generates all current elements in a temporary
directory and checks the thirteen serialized KF96 tables. Both test modules
passed by direct `unittest` execution.
No third-party test package is
needed.

The subsequent KF99 mean heat checks bring the fitting test total to 28.
The nine added tables passed 18,027 compiled samples with a maximum
relative interpolation error of `5.1541e-6`; all 27 compiled mean heat
samples matched the adopted published values. The full curve figure,
approximation limits, and controller records are linked in
[mean heat validation](kf99_mean_heat_validation.md).
