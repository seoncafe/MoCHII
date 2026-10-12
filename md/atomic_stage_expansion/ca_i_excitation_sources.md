# Ca I excitation source and level-model requirements

Date: October 10, 2026. This audit does not install a Ca I atom or infer
electron collision strengths from radiation data.

## Published calculation and numerical access

Zatsarinny, Parker, and Bartschat, [*Physical Review A* **99**, 012706
(2019)](https://doi.org/10.1103/PhysRevA.99.012706), calculate electron
collisions among 39 physical Ca I target states from threshold to 100 eV.
Their BSR-483 calculation adds 444 continuum pseudostates to the BSR-39
model. The [APS version-of-record abstract](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.99.012706)
confirms the article identity and coverage. The typeset final-journal text,
[publicly posted by the author](https://www.researchgate.net/publication/330330912_Electron-impact_excitation_and_ionization_of_atomic_calcium_at_intermediate_energies),
was readable on the hosting page. Its Section II describes term-dependent
MCHF target orbitals and close-coupling models; the calculations resolve
orbital angular momentum and total spin. The listed example targets are LS
terms such as `4s4p 3P` and `3d4s 3D`, rather than individual J levels.
Section III compares BSR-39 with BSR-483 and reports that continuum coupling
can materially change excitation cross sections, including the resonance
transition. The paper estimates about 5% uncertainty for the ground-to-
`4s4p 1P` resonance transition but says independent tests are needed for
reliable uncertainty estimates for all other transitions. That 5% must not be
applied to every matrix entry.

Section IV states that the entire numerical dataset is **available in
electronic form upon request**. The paper presents only selected curves, so
the 39-state collision matrix cannot be reconstructed from the article.
The APS article landing page exposes a PDF link and an author-accepted
manuscript link but no supplemental numerical file. The APS publisher PDF
request returned HTTP 401; the legacy `link.aps.org/pdf` URL redirected to
the abstract HTML. The author-posted typeset PDF download returned HTTP 403,
though its text was available in the hosting page. A private shell attempt
to reach the ADS API failed DNS resolution, and the ADS article-page request
returned HTTP 429; the DOI and final citation were instead checked against
the APS page. The accepted manuscript and arXiv copy were not used as substitutes
for the final-journal methods. No publisher PDF or numerical data file was
archived because neither was actually retrieved. No author was contacted.
The actionable access path is an authorized request to the surviving authors
for the complete BSR-483 state-to-state cross sections, target-state index
and energies, units, incident-energy grid, and redistribution terms.

## Local radiation and identity constraints

The local Cloudy c25.00 Stout Ca I files are
`/home/kiseon/CLOUDY/c25.00/data/stout/ca/ca_1/ca_1.{nrg,tp,coll}`.
The reproducible [inventory script](validation_2026-10-10/ca_i_excitation_sources/audit_stout_ca1.py)
and [result](validation_2026-10-10/ca_i_excitation_sources/stout_ca1_inventory.json)
record source hashes and counts. There are 780 J-resolved energy levels,
136 radiative `A` records, and zero tabulated electron-collision records;
`ca_1.coll` says `Use gbar`. Only 58 `A` records have both indices at most
39. Of the first 39 *Stout J levels*, ten have no radiative exit among the
listed records and 23 have no radiative path to the ground. Those first 39
J levels are **not** the 39 LS targets of the 2019 calculation.

For example, Stout index 1 is `4s2 1S0`; indices 2–4 are the three distinct
J levels of `4s4p 3P`, with energies 15157.901, 15210.063, and
15315.943 cm^-1. The 2019 cross section for the `4s2 1S` to `4s4p 3P`
term transition is therefore not a direct rate for any one of indices 2–4.
Assigning its total cross section to each J level, dividing it statistically,
or equating the paper's state 2 with Stout level 2 would invent fine-structure
physics. The Stout `A` file is also too sparse to serve as a complete
radiative decay network for the paper's 39 terms. The local Stout provenance
is NIST energy and radiative records dated 2013-10-01; it is not a source
of BSR-483 electron excitation data.

## Physically defined candidate after source delivery

1. Obtain the published BSR-483 numerical cross sections with their 39 LS
   target identifiers and incident-energy convention. Verify that they are
   total excitation cross sections, rather than apparent emission cross
   sections that include cascades, and match the final article's selected
   curves and target thresholds. Preserve the received original and record
   its origin and redistribution permission.
2. Build either a consistent **term-resolved** Ca I model using complete
   term-to-term radiative decay data, or obtain a source calculation that
   resolves the necessary J transitions. A mixed LS collision/J radiation
   model needs a published recoupling or statistical-equilibrium derivation
   with a stated validity interval; no such derivation is established here.
   Restricting to fewer terms also requires an explicit cooling-convergence
   test and a radiatively connected retained set.
3. Maxwellian integrate each verified cross section from its physical
   threshold, then calculate deexcitation by detailed balance using the same
   level degeneracies and energies. Evaluate whether the 100 eV upper limit
   and any high-energy continuation adequately bound the integral across the
   intended 1,000–100,000 K range; its adequacy cannot be assumed from the
   source endpoint alone. Check positivity, near-threshold resolution,
   term-index symmetry, and agreement with any published temperature rates.
4. Validate a source-consistent stationary population and emitted cooling
   over the temperature-density domain before registration. Keep Ca I
   absent from production electron cooling until the collision matrix and
   radiative closure are established. A Stout `gbar` calculation cannot
   replace the missing BSR-483 numerical matrix.

The present block is **numerical source access plus a coherent level/radiation
mapping**, not absence of a published electron-collision calculation. The
2019 paper supplies a strong physical candidate, but no importable Ca I
excitation atom has been demonstrated from the accessible records.
