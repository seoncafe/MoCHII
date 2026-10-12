# Remaining excitation-energy source assessment

The 16 remaining CHIANTI ions cannot all use the current observed-when-present
energy vector as a complete radiative atom. O IV, Ne V, S VII, Ar V, Mg IV,
and Fe VII contain positive Einstein-A rows whose indexed upper level falls
below the indexed lower level under that vector. Changing collision-pair
orientation is physically valid because effective collision strength is
symmetric in the pair. Changing an Einstein-A row's upper state is a separate
atomic-data decision and must not happen implicitly.

The direct raw-file audit is
[`source_energy_audit.json`](validation_2026-10-10/remaining_energy_audit/source_energy_audit.json).
It records SHA-256 values for each `.elvlc`, `.wgfa`, and `.scups` file, all
nonpositive indexed gaps, and the Cloudy c25.00 NIST ASD 2022-11-15
ionization threshold used here. Counts below include positive-A rows only.
All numbers were recomputed from the files on 2026-10-10; no cooling product
was generated.

| Ion | Levels | Bound levels under mixed energy | Observed values | Initial consecutive observed levels | Reversed/equal A: mixed | Reversed/equal A: theory |
|---|---:|---:|---:|---:|---:|---:|
| C V | 577 | 127 | 168 | 39 | 0 | 0 |
| N V | 923 | 63 | 37 | 18 | 0 | 0 |
| N VI | 577 | 127 | 47 | 25 | 0 | 0 |
| O IV | 204 | 148 | 133 | 79 | 4 | 0 |
| Ne IV | 22 | 22 | 22 | 22 | 0 | 1 equal |
| Ne V | 304 | 304 | 153 | 33 | 40 | 0 |
| S VI | 161 | 32 | 56 | 32 | 0 | 0 |
| S VII | 209 | 177 | 95 | 41 | 4 | 0 |
| Ar V | 220 | 197 | 43 | 5 | 6 | 0 |
| Ar VI | 202 | 182 | 37 | 15 | 0 | 0 |
| Ar VII | 196 | 196 | 72 | 38 | 0 | 0 |
| Mg IV | 195 | 157 | 105 | 63 | 39 | 0 |
| Fe V | 34 | 34 | 34 | 34 | 0 | 0 |
| Fe VI | 96 | 96 | 80 | 70 | 0 | 0 |
| Fe VII | 189 | 189 | 113 | 17 | 21 | 16 |
| Ca VI | 5 | 5 | 5 | 5 | 0 | 2 equal |

The source model must first retain levels strictly below the physical
ionization threshold. C V and N VI each have 127 such levels; N V has 63;
S VI has 32. The N V `.elvlc` footer explicitly calls its higher levels
"above IP" and cites Goryaev et al. (2017) and Kato et al. (1997) for
doubly excited states. The C V and N VI footers make the same distinction.
The S VI footer separates outer-shell and inner-shell calculations; its 129
levels above the measured threshold cannot enter a bound-only radiative
cooling model without autoionization and competing channels. Other ions with
above-threshold levels are listed in the JSON. This threshold is a physical
model boundary, not a cutoff chosen to improve numerical fit quality.

For C V, N V, N VI, Ne IV, S VI, Ar VI, Ar VII, Fe V, Fe VI, and Ca VI, the
mixed observed/theoretical vector has no nonpositive radiative gap among the
bound levels. It is an internally directional starting vector, not an
accuracy certificate. Fe V and Ca VI show why wholesale theoretical energies
are inappropriate: seven Fe V collision pairs reverse, and two Ca VI
positive-A pairs become exactly degenerate. Collision inversions in otherwise
valid vectors should be oriented by the chosen energies while retaining the
same effective collision strength and source scaling energy.

S VI needs a fine-structure correction before using its 32 bound levels.
CHIANTI records exactly equal observed energies for bound pairs 11/12,
18/19, 20/21, 27/28, 29/30, and 31/32, although each pair has a nonzero
`.scups` scaling gap. The corresponding theoretical separations are 10, 5,
3, 3, 2, and 2 cm^-1. Cloudy Stout `s_6.nrg` (SHA-256
`c6cc5e64dd43d947f3da433967a4cfba2e845c41988c17c604274c2a0d7760d4`)
puts 0.001 cm^-1 between each pair; this tiny inserted separation has no
stated experimental support. Moreover, Stout indices 20/21 and 31/32 reverse
the J/statistical-weight identity of the corresponding CHIANTI indices.
Stout energies cannot be copied by index. A physically motivated candidate is
to retain each observed pair's statistical-weight term centroid and apply
its CHIANTI calculated fine-structure separation with the CHIANTI J
assignment. For an unresolved group, the selected energy is
`E_obs + E_theory(i) - sum_j[g_j E_theory(j)] / sum_j[g_j]`. This requires an
explicit source note that the splitting is theoretical and a radiative and
collision-gap audit after construction. Two additional zero collision gaps,
78/80 and 83/84, occur above the S VI ionization threshold and disappear
from the bound model.

For O IV, Ne V, S VII, Ar V, and Mg IV, the full theoretical vector has
positive gaps for every retained positive-A row. It is the conservative
source-consistent candidate for a complete atom. A low-level observed/theory
splice can preserve precise lines, but its cut must follow checked level
identities and term groups. Exhaustive cut scans found the largest index cuts
with no positive-A reversal: O IV 106, Ne V 91, S VII 176, Ar V 96, and Mg IV
58. These numbers are *diagnostic bounds*, not selected physical cutoffs.
For example, Ne V 91/92 splits a `2s2 2p 4s` term, and Ar V 96/97 splits a
`3s2 3p 4d` manifold. The Mg IV 58/59 boundary separates `3d` and `4s`
configurations, but discards observed 4s levels. The O IV 106/107 boundary
separates `4d` and `4f`; S VII 176/177 separates one observed `7d` level from
other 7d levels. The maximum index alone cannot justify any of these splices.
Using a complete theoretical vector is therefore the immediately testable
coherent choice for those five ions; retaining observed low levels requires
further identity and wavelength checks, especially at low temperature.

A configuration-based candidate now retains observed energies only in the
ground electronic configuration and uses original CHIANTI theoretical
energies for all other bound levels. The exact vectors and source hashes are
in [`ground_configuration_splice.json`](validation_2026-10-10/remaining_energy_audit/ground_configuration_splice.json).
The observed indices are O IV 1--2, Ne V 1--5, S VII 1, Ar V 1--5, and
Mg IV 1--2. All five vectors have zero nonpositive radiative gaps. O IV,
Ne V, S VII, and Ar V also have zero nonpositive collision gaps. Mg IV has
one **equal** collision gap: theoretical levels 135 and 136 both have
788126 cm^-1, despite a source `.scups` row with nonzero scaling energy.
Thus the simple ground-configuration Mg IV vector is incomplete.

The more physically informative candidate preserves the complete observed
low valence complex: all `n=2` configurations for O IV (levels 1--15),
Ne V (1--20), and Mg IV (1--3), and the `n=3` complex for S VII and Ar V.
Above it, levels use the original CHIANTI theoretical values, except where
the source gives identified fine-structure information. The exact vectors,
source hashes, and row checks are in
[`valence_complex_resolved.json`](validation_2026-10-10/remaining_energy_audit/valence_complex_resolved.json).
All five have zero nonpositive radiative and zero equal collision gaps.
O IV, Ne V, and Ar V respectively have 2, 6, and 6 collision rows whose
original index direction reverses; their effective collision strengths
remain the same after physical orientation. S VII has no reversed collision
row after its identified 23/24 correction. Mg IV has 149 reversed collision
rows because its measured `4f` manifold is retained while other high
configurations follow the theoretical target. This mixed high-level ordering
requires full population and cooling validation; the number of reversed
collision indices is not itself an energy defect.

S VII levels 23 and 24 have the same CHIANTI observed value of 1644609
cm^-1. Cloudy Stout `s_7.nrg` identifies the same `3D1` and `1D2` levels
as 1644599 and 1644674 cm^-1. Substituting these identified energies
removes the equal collision gap without changing an Einstein-A direction.
Mg IV CHIANTI levels 135/136 have identical theoretical energies but
distinct observed energies, 806911.500 and 806909.188 cm^-1. Their J
identities and approximately equal energies are independently present as
Stout levels 109/108. The candidate retains all identified observed `4f`
energies so that the pair is not isolated from its configuration. Its
observed 133/134 pair is unresolved; the original CHIANTI theoretical
5 cm^-1 fine-structure separation is centered on their common observed
statistical-weight term centroid using the same formula as S VI. The resulting
selected splitting estimates the unresolved fine structure, rather than
inventing an arbitrary 0.001 cm^-1 gap. The selected atom still needs
convergence and source-sensitive cooling checks.

The negative `.wgfa` wavelength on each of the mixed-vector reversed rows
for those five ions flags an approximate wavelength in CHIANTI's own line-list
documentation. The O IV and Mg IV `.wgfa` file notes state that the compilers
already removed other transitions with observed upper energies below lower
energies. Thus a negative wavelength does not authorize changing the Einstein
A direction. It does support using the theoretical target vector while
checking any later observed-energy splice against every radiative row.

The CHIANTI Fe VII mixture needs its own source reconstruction. Its mixed vector reverses 21
positive-A rows, and its theoretical vector still reverses 16. In the mixed
vector, level 117 has 554703.9 cm^-1 and level 126 has 548228 cm^-1, yet
`.wgfa` labels 117 to 126 with A = 33500 s^-1 and a **positive** wavelength
of 15441.869 Angstrom, whose inverse is approximately 6476 cm^-1. The
observed gap is -6475.9 cm^-1. This is a direct inconsistency between the
indexed A direction and the measured level ordering, even though the
wavelength magnitude agrees. The `.elvlc` footer says its theoretical values
come mostly from Li et al. (2018) for 134 levels, with Tayal and Zatsarinny
(2014) filling others, plus adjustments for levels 25--27 and 29. Li et al.'s
published journal article indeed calculates 134 fine-structure levels, not
this full 189-level vector. Fe VII levels 1--125 contain no reversed A under
the mixed vector, but truncation there must still pass complete temperature
and density convergence tests; it cannot be declared a complete atom.
Cloudy Stout `fe_7.nrg` has different high-level positions and labels, so it
cannot be substituted by index.

Cloudy Stout c25.00 provides a separate coherent full Fe VII atom. Its
`fe_7.nrg`, `.tp`, and `.coll` SHA-256 values are
`42a060bca5c3dab121ca5c5cc4102235fc28976a1b1635962cafc3799f8c4415`,
`d3c426196550ccfe986b55c765d5407b27e1b1be6f5340678504a7dd0a417f74`,
and `2fa169e62921cafda0490659a95d8a2d11183a3e14ad4635939879f321e0a72d`.
Direct parsing finds 189 strictly increasing bound energies, 17,766 distinct
pairs in both radiation and collision files, 17,531 positive Einstein A
values, and no nonpositive radiative or collision gap. The collision file
gives 13 temperature nodes from 9800 to 9.8e7 K. These files attribute the
Fe6+ ICFT atom to Witthoeft and Badnell (2008, A&A 481, 543--551). The
[author repository record](https://strathprints.strath.ac.uk/19746/)
describes an 89-term, 189-level scattering target. The final typeset journal
article has now been obtained from the author's public server and checked;
see [final-paper verification](fevii_final_paper_verification.md). All 189
energy rows and 57 printed radiative rows were compared directly. The first
nine Stout energies match the observed column and the others the calculated
column at the printed precision. Every deposited collision-strength value
has not been independently authenticated against an original author archive.

`tools/fitting/stout_excitation_atoms.py` reads those three files as one atom
and checks hashes, row counts, bound energies, rates, and the temperature
grid. At 1000--9800 K it holds only Upsilon at its 9800 K source value;
the Maxwellian rate still uses actual temperature. This explicit
approximation has unknown resonance accuracy below 9800 K. The Fe VII Stout
model does not reassign any CHIANTI Einstein-A row.

Local MOCASSIN v2 `data/fevii.dat` has only nine Fe VII levels, whose
energies match the first nine Stout levels. Its invoked `hydro_mod.f90`
reader imports the packaged energy and A tables, and `emission_mod.f90`
uses indexed energy gaps for the Boltzmann factor and photon energy. The
nine-level model does not establish a policy for the conflicting CHIANTI
189-level mixture.

Before accepting any atom, use one selected energy for each retained level in
Boltzmann factors, detailed balance, emitted photon gaps, and the finite-density
rate matrix. Keep `.scups` `de` as its Burgess-Tully scaling energy. Reject a
nonpositive radiative gap, and account explicitly for retained or excluded
collision-connected states. Use the bound-level subset in both the direct
source-cooling calculation and finite-density serialization, then compare full
cooling and population curves. Fe VII can use the complete Stout atom,
subject to convergence and the stated cold-strength approximation. No public
interface change is needed to
make the source audit or collision orientation correction; a change to the
meaning or layout of existing atom files would require a separate decision.

Source access: The NASA ADS API initially failed DNS resolution in the
restricted environment. An authorized network retry verified Fe VII bibcode
`2008A&A...481..543W`, DOI `10.1051/0004-6361:20079073`, journal title,
and year. The attempted final A&A publisher PDF returned HTTP 403, its ADS
publisher gateway returned HTTP 403, and the ADS PDF timed out. An additional
authorized request obtained the final journal PDF from
`https://amdpp.phys.strath.ac.uk/tamoc/PUB/2008/aa_481_543.pdf`; the
physical methods and printed numerical tables have now been checked. The
[published Li et al. 2018
journal page](https://academic.oup.com/mnras/article/479/1/1260/5035830)
was inspected for its 134-level method and scope. The original CHIANTI file
footers were read for the other paper attributions. Their final journal
methods and numerical data were not independently checked here, so conclusions
about their accuracy, level mixing, or resonance placement remain unverified.
