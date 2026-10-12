# Mg III full-bound radiative consistency candidate

**Physical source correction required:** official NIST comparison subsequently
identified a misassigned observed energy at source level 167. The candidate
and numerical profiles described below are preserved conditional diagnostics
and are physically stale for production acceptance. The J=0 `7p 1S` state
must use its own theoretical 634430 cm^-1 energy instead of the observed
628148.375 cm^-1 energy of a J=1 `7d` state. This is a source-identity defect,
not an atomic uncertainty inferred from a small cooling difference.

This is a private physical and numerical preparation, not an accepted full-ion
cooling model. The source is the local CHIANTI Mg III atom. The selected bound
set has 167 collision-connected levels, 3,364 positive spontaneous rows,
13,861 electron-collision pairs, and 207,915 original Burgess–Tully nodes.
The first 27 levels reproduce the existing restricted atom: 117 radiative rows
and 351 collision rows with 5,265 nodes. Their serialized level, radiative,
and collision records are identical in the private first-27 and full
dipole-strength candidate files.

## Radiative rule and its physical limits

The source `mg_3.wgfa` attributes transitions within the first 27 levels to
Jönsson et al. (2014) and the remaining transitions to Liang and Badnell
(2010). For the 3,135 higher-level rows with positive `gf`, opposite parity,
and allowed electric-dipole angular momentum, the frequency inferred from
source `A`, `gf`, and the source upper-level weight agrees with the source
theoretical energy gap to within 0.386% (maximum). The 16 rows whose source
index order is opposite to the selected physical energy order all belong to
this class; their maximum corresponding frequency discrepancy is 0.180%.
These measurements support one explicit electric-dipole transformation for
the higher-level rows. The tabulated `gf` is a weighted oscillator strength,
proportional to transition frequency times electric-dipole line strength. It
cannot remain fixed when the photon frequency changes. With source theoretical
gap \(\Delta E_{\rm src}\) and selected physical gap
\(\Delta E_{\rm phys}\), preserving dipole strength gives:

\[
 A_{\rm phys}=A_{\rm src}
 \left(\frac{\Delta E_{\rm phys}}{\Delta E_{\rm src}}\right)^3
 \frac{g_{\rm src,upper}}{g_{\rm phys,upper}}.
\]

Here the selected physical energies identify the actual emitting level and
photon frequency. The source `A` and theoretical gap define the adopted dipole
strength; source `gf` corroborates the source frequency convention. This
transformation does not claim that the calculated wavefunctions
or fine-structure mixing are exact. For 68 of the 140 levels above 27, an
observed energy is available; the other 72 retain a theoretical energy. Some
observed/theoretical differences reach 6,281.625 cm^-1. For example, the
selected gap of source levels 127 and 143 is 8.688 cm^-1 while their theoretical
gap is 5,818 cm^-1. Such near degeneracy makes the frequency correction large,
but it is not a reason to keep a spontaneous rate for a photon with the wrong
energy. It does require verification that the observational assignment and
source wavefunction label identify the same pair before production adoption.

The 117 Jönsson rows within the first 27 levels retain their source A values.
Their 110 E1 rows do not follow the same simple source A/gf theoretical-frequency
relation: the median inferred-frequency discrepancy is 5.85%, and some values
are much larger. Therefore the Liang–Badnell rescaling is not applied to this
separate source. The remaining 112 higher-level positive A rows have `gf=0`:
109 connect same-parity configurations and three connect opposite-parity
configurations. Their multipole identities and line strengths are not present
in the local table. The private candidate retains each original A as an
**unidentified-multipole effective-rate approximation** and marks that choice
at the atom header. All 112 rows retain the physical selected energy for
emitted photons. Their summed A contribution reaches at most 0.412% of one
emitting level's total candidate A, but that small numerical fraction is not
a physical verification of the approximation. No higher-level A source was
silently substituted.

Across the 3,135 corrected E1 rows, the candidate/source A ratio ranges from
`3.330e-9` to `462.08`; 146 rows lie outside a factor of two. Summed A for a
physical emitting level changes by a factor from `0.6451` to `1.3800`. These
are direct measurements from the private candidate; they are sensitivity
information, not an accuracy estimate. Source radiative rows are not merely
reversed at read time: their A depends on photon frequency and emitting
statistical weight, unlike symmetric Upsilon.

The original collision nodes, scaling gaps, types, interpolation choices, and
adopted 1,800 K source temperature floor are unchanged. Collision-strength
accuracy and high-state convergence remain separate physical questions. The
source ranges printed by the current writer are 1,800–18,000,000 K for the
identified BT collision calculation; below 1,800 K the adopted Upsilon is held
at its lowest published-temperature value while Maxwellian rates use the
actual temperature. The present numerical use under discussion is
1,000–100,000 K.

## Private files and measured preflight

The files under
[`validation_2026-10-10/mgiii_radiative_candidate`](validation_2026-10-10/mgiii_radiative_candidate/)
contain the unaltered-source serialized full atom, the first-27 atom, and the
dipole-strength-corrected full candidate. `build_candidate.py` rechecks raw
source hashes before and after writing, reads the preserved private full atom,
writes the corrected candidate, and records
class counts and rate changes in `result.json`. The source hashes are:

| Raw input | SHA-256 |
| --- | --- |
| `mg_3.elvlc` | `07e61cf56c782806986ca9a19c2197f2cfa85adaf568be4b0ab711850084851e` |
| `mg_3.wgfa` | `bb00af4bf3369c4facb78b068e29df4d7dfb9c0a042ae2f21ca377bed44ab841` |
| `mg_3.scups` | `9aba89f5d1f6ca158a88547420377d988df1363e98d0cbd0c948fb444a302402` |

The first private candidate incorrectly held the oscillator strength `gf`
fixed, resulting in a frequency-squared A transformation. The initial atom,
script, results, native strength sample, and PDF are preserved under
[`rejected_fixed_oscillator_strength`](validation_2026-10-10/mgiii_radiative_candidate/rejected_fixed_oscillator_strength/)
with a hash manifest and explicit rejection note. The current candidate uses
the frequency-cubed rule above. `verify_candidate.py` independently compares
\(A g_{\rm emitter}/\nu^3\) between source and candidate for all 3,135
corrected E1 rows; the maximum relative discrepancy is `6.66e-16`. It also
checks unchanged level and collision records, all 117 first-27 radiative
rows, all 351 first-27 collision rows, 112 unidentified-multipole A values,
and the three original CHIANTI source hashes. The corrected atom SHA-256 is
`46e71f59c1c885a46ec4304f4a26bc1549e716a5490c6118d4cd242b16711873`.
The current source writer rejects an uncorrected reversed radiative A. Thus
the preserved private atom is a diagnostic input, not a product that the
current writer can or should recreate from raw Mg III radiative rows. Its
header records the writer hash at the time it was generated; the controller
separately freezes the current reader and writer hashes.

The 112 higher positive-A rows with `gf = 0` lack source multipole labels and
line strengths. An angular-momentum and parity audit of the source level
assignments finds 65 rows for which E2 is the lowest allowed multipole, 44
rows allowing M1 and E2, and three opposite-parity rows with
`|Delta J| = 3` for which E3 is the lowest allowed multipole. These are
selection rules, not identifications of the calculated transition type.
Consequently, the candidate retains all 112 source A values without a
frequency correction. The row-level audit is in
[`multipole_selection_audit.json`](validation_2026-10-10/mgiii_radiative_candidate/multipole_selection_audit.json).

The private native driver linked the current `src/nlevel_mod.o` without
modifying the source or object. It read all 13,861 candidate collision rows at
10,000 K and reproduced the Python source Upsilon values with maximum
relative error `1.55e-15`. In the subsequent full-profile test, 69,305
Upsilon samples at five temperatures agreed to `2.66e-15` relative. An older
optimized driver linked before `BT_TMIN_K` support rejected the new atom's
metadata as an unknown section; the current-object driver was used. The
[static diagnostic figure](validation_2026-10-10/mgiii_radiative_candidate/radiative_correction_diagnostic.pdf)
was visually inspected; it shows the full range of E1 rate corrections and
the unidentified-multipole A fractions by emitting level.

| Input condition | Full candidate | First-27 comparison |
| --- | --- | --- |
| Level energies | Selected observed when available, otherwise theoretical | Exact first-27 subset |
| Collision strengths | Original 13,861 BT source pairs | Exact 351-pair subset |
| Spontaneous rates | 117 original Jönsson, 3,135 corrected higher E1, 112 original unidentified multipoles | 117 original Jönsson |
| Temperature interval tested | 1,000–100,000 K; adopted constant Upsilon below 1,800 K | Same |
| Electron densities tested | `1e-6`, `1e2`, `1e6`, `1e10 cm^-3` | Same |
| Radiation field and other colliders | None in this statistical-equilibrium test | Same |

## Full-profile numerical result and remaining physical questions

The approved private calculation used 801 log-spaced temperatures from 1,000
to 100,000 K and four electron densities. A finite persistent controller
froze 26 source, code, object, and input hashes. The first attempt stopped at
a mistyped source path before producing numerical profiles. The recovery
controller completed the pilot, compiled native profile, independent Python
profile, and comparison commands. Its final artifact check referenced
`compare_step.json`, while the successful comparison wrote
`comparison_step.json`; it therefore recorded a failed final step despite
the command's zero exit status and the result's `passed = true`. The original
controller state and products remain intact. A separate
[`audit_completed_results.py`](validation_2026-10-10/mgiii_radiative_candidate/recovery_path_fix/audit_completed_results.py)
verified the preserved source/input hashes, seven product hashes, file shapes,
finite nonnegative cooling, population normalization, and comparison result
without rerunning a numerical profile. Its
[`audit_completed_results.json`](validation_2026-10-10/mgiii_radiative_candidate/recovery_path_fix/audit_completed_results.json)
records this distinction.

Native cooling agreed with independent Python statistical equilibrium to at
most `1.13e-14` relative for the full atom and `1.55e-15` for the first 27
levels. The maximum absolute population difference was `2.11e-15`. The
smallest positive direct cooling value was `1.82e-284 erg cm^3 s^-1`. The
first-27/full cooling difference reached `17.20%`, `17.20%`, `17.21%`, and
`17.24%` at the four densities, so the 27-level truncation fails the `0.3%`
criterion **within this corrected candidate model**. The complete
[cooling profiles](validation_2026-10-10/mgiii_radiative_candidate/full_and_first27_cooling_profiles.pdf)
and [relative-difference profiles](validation_2026-10-10/mgiii_radiative_candidate/first27_truncation_differences.pdf)
were visually inspected. The differences rise toward the upper temperature
boundary. The uncorrected full atom is not a coherent physical reference for
cutoff acceptance. Native and direct agreement establishes implementation
consistency; it does not establish accuracy of the underlying radiation or
collision source.

The high-lying observed level assignments, especially near-degenerate
6p–6d/7p–7d pairs, and the 112 unidentified-multipole effective rates remain
physical uncertainties. A final publisher PDF for Liang and
Badnell (2010), *Astronomy & Astrophysics* **518**, A64,
DOI `10.1051/0004-6361/201014170`, was not obtained in this work; an earlier
downloaded `/tmp` file is HTML and is not evidence from the final paper. The
Jönsson 2014 source verification is a separate outstanding check. The
rejected-candidate calculation changed no production atomic data.

## Corrected level-167 identity and current private numerical validation

The preceding 167-level candidate is retained as a **physically rejected
diagnostic**. Its selected level-167 energy, 628148.375 cm^-1, came from a
misassigned CHIANTI observed-energy field. CHIANTI labels that level
`2s2 2p5 7p 1S`, `J = 0`, with theoretical energy 634430 cm^-1. The
[official NIST level response and row-level audit](validation_2026-10-10/mgiii_radiative_candidate/source_identity/README.md)
identify the 628149 cm^-1 level as `7d`, `J = 1`, which has different parity
and statistical weight. NIST lists 111 levels; 67 of the other 68 higher
CHIANTI observed-energy levels match NIST energy within 1 cm^-1, `J`, and
outer orbital. This supports their level identities but does not establish
wavefunction accuracy or the correctness of all transition probabilities.
CHIANTI level 127 has no observed energy, so its selected 620693 cm^-1
theoretical energy and its 8.688 cm^-1 gap to level 143 are not an
independently observed interval.

The [corrected private source-node atom](validation_2026-10-10/mgiii_corrected_level_identity/README.md)
selects the same source level 167 at its theoretical 634430 cm^-1 energy.
Its SHA-256 is
`ca27317e9a76235bcd45243c61ab4659fcfb562928c917f25606c4dd4197764f`.
All other selected energies, 13,861 collision records and their 207,915
nodes, and the first-27 atom remain unchanged. Twenty-four serialized
radiative rows involving level 167 changed; an independent reconstruction
from raw `.wgfa` rates agreed exactly with the new A values. The source
writer's focused Mg III test passed with the corrected atom.

For the corrected full atom, the optimized native calculation and independent
Python statistical equilibrium agreed over 801 temperatures from 1,000 to
100,000 K and four densities from `1e-6` to `1e10 cm^-3`. Their maximum
relative cooling difference was `2.58e-14` and maximum absolute population
difference was `2.00e-15`. The maximum stationarity residual normalized to
total gross collisional flow was `3.49e-15`; collisional energy input minus
radiative output, normalized to gross collisional energy flow, reached
`4.14e-14`. The [complete profiles](validation_2026-10-10/mgiii_corrected_level_identity/full_and_first27_cooling_profiles.pdf)
and [truncation differences](validation_2026-10-10/mgiii_corrected_level_identity/first27_truncation_differences.pdf)
were visually inspected. The first 27 levels differ from the corrected full
atom by up to `17.19%` in cooling and remain an unacceptable truncation for
the `0.3%` criterion. The maximum cooling difference between the rejected
and corrected full candidates was `0.0687%` on this grid; that small effect
does not validate the wrong level identity. Numerical closure does not verify
the high-level radiative source or collision calculation. The complete
private cooling-product controller is tracked separately under
[`mgiii_cooling_products`](validation_2026-10-10/mgiii_cooling_products/).

## Production installation after main review

The corrected full167 atom and fit are now installed in `data/atomic`.
Main review inspected the complete source and fit PDFs and verified all
40 and 44 frozen input hashes in the completed controllers before later
header-only generator comments changed those current source files. Numerical
records in the installed products match the validated private products
exactly; headers additionally state the ID167 correction. The actual
production-module native driver opened the installed atom at 10,000 K and
four densities and reproduced the validated profile exactly. No full source
profile or registered table was regenerated for installation.
[Installation manifest](validation_2026-10-10/mgiii_corrected_level_identity/installation_manifest.json)
and [main review](validation_2026-10-10/mgiii_corrected_level_identity/main_review.json)
record installed hashes and the limited file-open check. Production has
34 ordinary metal electron-excitation fits among 43 stages/32 transitions;
this installation does not activate the staged 64-stage roster.
