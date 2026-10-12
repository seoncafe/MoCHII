# MOCASSIN v2 atomic coverage absent or inactive in MoCHII

Current state verified on October 11, 2026 (KST); initial inventories below
were measured on October 10. The revised first-phase goal is 77 metal
stages (66 adjacent transitions), covering stages I--VII of C, N, O, Ne, S,
Ar, Mg, Fe, Si, Cl, and Ca. H and He add five stages, for 82 total.
The earlier 64-stage target is superseded.

## Scope and evidence

The reference is the local `mocassin.2.02.73.2` source and its actual
`fileNames.dat` mapping. The [source inventory](source_inventory.json)
contains 56 available line-data files and 21 missing files within this metal
roster. These are declared data models, not 56 successfully executed models.
The [initial MoCHII comparison](mochii_coverage.json) was measured before missing-atom installation from
`ELEMENT`, `NLEV`, and cooling-file records and includes file hashes.
There was no matched spatial plasma calculation or full regression run.

MoCHII's installed production records still declare 43 metal stages and
32 adjacent transitions. Stage storage allows seven stages, and the current
metal-transition capacity is 66, derived from eleven elements and seven
stages. Capacity and generated coefficients do not establish activation.
`species_setup` reads cooling data only through each active element's declared
stage count; `lines_write` likewise iterates over active stages. The
`local_ne` cooling option uses the installed atom through
`nlevel_cooling_add`, `nlevel_cooling_solve`, and `nlevel_cooling_apply`.

## Missing installed atom and cooling files

Before the missing-atom installation, twenty-five MOCASSIN data models had
neither the corresponding installed MoCHII atom nor its cooling fit. This
table is the preserved initial snapshot. The current selected reference
directory contains 55 atoms with 571 levels and 49 primary cooling fits.
Mg I uses the approved five-level constant-strength model; only Ca VII
remains absent. Fits are not mandatory for the absolute-cooling consumer. Some have preparation files or validation
routines elsewhere; this table does not claim that no callable routine or
source data exists. An inactive ion also requires production stage activation.

| Ion | MOCASSIN declared levels | Active MoCHII stage |
| --- | ---: | --- |
| N IV | 10 | No |
| N V | 15 | No |
| O VI | 15 | No |
| Ne IV | 13 | No |
| Ne V | 15 | No |
| Ne VI | 10 | No |
| Ne VII | 10 | No |
| S V | 5 | No |
| S VI | 3 | No |
| Ar V | 5 | No |
| Ar VI | 2 | No |
| Ar VII | 5 | No |
| Mg I | 5 | Yes |
| Mg IV | 2 | No |
| Mg V | 9 | No |
| Mg VI | 13 | No |
| Mg VII | 15 | No |
| Fe V | 5 | No |
| Fe VI | 19 | No |
| Fe VII | 9 | No |
| Si VI | 2 | No |
| Si VII | 9 | No |
| Cl V | 2 | Yes |
| Ca VI | 5 | No |
| Ca VII | 27 | No |

## Smaller installed explicit level models

At the initial checkpoint, seven installed MoCHII atoms had fewer explicit
levels than the corresponding MOCASSIN file. The table below preserves that
checkpoint; all seven selected reference atoms now have the source level
counts after the approved paired replacements. This establishes a gap in explicit level-population and line
emission coverage, but level counts alone do not establish state identity or
collision-pair completeness. Identify configurations, terms, J values,
physical energies, and transition endpoints before importing additional data.

| Ion | MOCASSIN declared levels | MoCHII installed levels |
| --- | ---: | ---: |
| C III | 10 | 5 |
| N III | 10 | 5 |
| O I | 7 | 5 |
| Ne III | 9 | 5 |
| S IV | 12 | 5 |
| Fe II | 142 | 13 |
| Si II | 12 | 5 |

The low-density cooling fits and explicit finite-density atoms are separate
products. `chianti_cooling.cooling_low_density_limit` uses the bound excitation
atom and spontaneous cascades, rather than simply the installed truncated
`NLEV` model. Therefore a smaller explicit atom does not establish that the
corresponding high-level contribution is absent from the low-density cooling
fit. The existing Fe II fit, for example, identifies its CHIANTI low-density
source; its actual source-to-product coverage must be checked before claiming
that all transitions above level 13 are absent from total cooling. A present
fit is likewise not proof of finite-density completeness or current source
accuracy.

## Installed but inactive oxygen models

O IV (147 levels) and O V (166 levels) are installed in the explicitly
selected `extended` model and passed installed native and coefficient-reader
checks. The default `mocassin_v2` model has ten levels for each ion. Their current source-data models retain
explicitly unverified atomic accuracy. They exceed the corresponding v2
10-level file sizes, but the production oxygen chain still ends at O III.
They are therefore prepared physics, not currently invoked production ions.

## Reference-code defects and unavailable line data

The inspected v2 `constants_mod` sets `nForLevels=17` and
`nForLevelsLarge=142`. `emission_mod.equilibrium` stops when the atom exceeds
the supplied line-emission array. Consequently the supplied Fe VI 19-level
and Ca VII 27-level files exceed the default capacity when those ions are
invoked. Their declared levels remain appropriate intended coverage targets;
this audit does not claim successful execution with the unmodified defaults.

The Fe II packed radiation reader in `hydro_mod` stores `(lower, upper, A)`
records into `A(lower, upper)`, whereas the equilibrium equations consume
`A(upper, lower)`. The measured source file has 1,144 positive raw rates and
zero positive rates in the downward triangle under that assignment. This is
a source-level reader defect, not a measured complete-run cooling result.
Preserve the external source for comparison and construct the physically
oriented MoCHII model. Four positive Ca VII rates also have nonpositive energy
gaps in the supplied level ordering; do not reproduce negative photon energies.

The 21 absent MOCASSIN line files are not evidence of absent ionization stages,
and do not impose an invented cooling model as a parity requirement. In
particular C VII is fully stripped and has no bound-level line emission.
MoCHII already has several line atoms beyond the supplied v2 files, including
Ne I, Mg III, Si I, and Si V. Keep those physically justified models.

## Implementation order under the revised goal

1. Prepare all 77 stages and 66 adjacent transitions, complete and validate
   their required ionization, recombination, and charge-transfer data, increase
   the transition capacity, and treat the fully stripped carbon terminal stage
   physically. Do not activate an incomplete rate chain.
2. Complete the 25 missing installed models, distinguishing actual v2
   temperature support from extended-table approximations; retain validated
   larger MoCHII atoms. Account explicitly for the two v2 capacity defects.
3. Expand the seven explicit level models using verified physical state
   identities and collision/radiation endpoints. Preserve valid low-density
   cooling products unless their affected source path requires regeneration.
4. Validate the actual installed reader, populations, energy accounting,
   finite-density cooling, line emissivities, and production invocation before
   declaring completion. Any subsequent plasma comparison requires matched
   input conditions and full relevant profiles.

The existence of Auger tables or recombination-pumping expressions alone is
not a parity gap. Establish invocation and nonzero data before identifying a
process as calculated by v2 and absent from MoCHII. The present result is an
atomic coverage comparison, not an exhaustive comparison of both codes.


## Absolute atomic cooling integration

The production cooling consumer now supports an available atom without a
cooling fit, including Ne VI and Fe V. The focused linked check passed 528
samples and eight zero-density cases, with maximum matrix comparison error
0.124251%. See [the validation report](../validation_2026-10-11/atomic_cooling_without_fit/README.md).
An exponential fit is therefore not mandatory for these atomic calculations.
Ca VII installation and the expanded production roster remain unfinished.


## Installed missing atoms

The [historical coverage inventory](mochii_coverage_after_missing_installation.json)
measured two missing supplied atom files after 23 exact candidate files were
installed. Mg I has since been installed. The current
[selected-model inventory](../validation_2026-10-11/selected_stage_coverage/current_reference.json)
identifies Ca VII as the only missing supplied model. Existing atoms and rate chains were retained. Cl V's actual
production no-fit consumer passed 264 scoped samples with 0.00600158% maximum
matrix comparison difference. See
[the root-reviewed installation](../validation_2026-10-11/missing_atom_installation/README.md).
The earlier 25 missing exponential fits and seven smaller explicit atoms
are historical counts, not current remaining work. Expanded production activation and source accuracy remain unfinished.
