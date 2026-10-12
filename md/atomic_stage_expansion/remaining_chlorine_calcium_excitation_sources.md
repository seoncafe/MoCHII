# October 11 update to the missing-atom snapshot

The following two-missing-atom inventory was measured before the approved
five-level Mg I model was installed. Only Ca VII remains absent from the
default `mocassin_v2` selection. The source audits below retain their dated
physical qualifications. See the [current result](../ATOMIC_IONIZATION_STAGE_EXPANSION_RESULTS_2026-10-11.md).

# Historical missing-atom installation status

A fresh direct file inventory after installing 23 absent MOCASSIN source
models finds two missing atoms, Mg I and Ca VII; 25 exponential fits remain
absent and seven installed models remain smaller than their supplied v2
counterparts. The production roster remains 43 metal stages and 32 adjacent
transitions. Cl V's installed two-level atom now contributes absolute cooling
without a fit; 264 actual linked samples passed with maximum independent
matrix difference 0.00600158%. This supersedes earlier missing-Cl-V-file
statements below, but not its published-data accuracy qualifications or
higher-level incompleteness. See [the current measured inventory](phase1_77_mocassin_v2/mochii_coverage_after_missing_installation.json)
and [the root-reviewed installation](validation_2026-10-11/missing_atom_installation/README.md).

# Electron-excitation source audit: Cl I, Cl V, and Ca I

Date: October 10, 2026. This audit identifies source physics and actual code paths; it does not install an atomic model. The local numeric checks are reproducible with [`check_local_sources.py`](validation_2026-10-10/neutral_chlorine_excitation/check_local_sources.py), which writes [`result.json`](validation_2026-10-10/neutral_chlorine_excitation/result.json).

## Cl I: the Stout electron label does not identify an electron rate

Cloudy c25.00 `data/stout/cl/cl_1/cl_1.coll` contains `RATE ELECTRON 1 2` with nine identical values of `4.7e-8 cm^3 s^-1` on 30, 100, 250, 500, 1000, 1500, 2000, 2500, and 3000 K. Its two levels have `g=(4,2)` and an observed separation of 882.351 cm^-1; `cl_1.tp` gives `A_21=0.0123 s^-1`. In `source/atmdat_chianti.cpp`, the `RATE ELECTRON` token selects the electron collider and stores the quantity as a rate rather than a collision strength. `source/species2.cpp::StoutCollRate` interprets it as the **downward deexcitation coefficient**. The level-balance code derives the upward coefficient using `g_upper/g_lower exp(-DeltaE/kT)` at the actual temperature. Its endpoint rule holds the first or last tabulated coefficient constant outside the Stout temperature grid. These statements describe Cloudy's implementation, not validation of the input physics.

The cited final article is [Hollenbach and McKee (1989), *Astrophysical Journal* **342**, 306, DOI 10.1086/167595](https://doi.org/10.1086/167595). The original final journal scan was obtained from [NASA ADS](https://ui.adsabs.harvard.edu/link_gateway/1989ApJ...342..306H/ADS_PDF) and inspected visually, including printed p. 333, Appendix C, Table 8 and footnote f. The Cl I row is `4.7(-8)` in a column headed gamma-e, but footnote f explicitly says, “Proton rates are substituted for electron rates for neutral target atoms.” This footnote applies to Table 8 itself, which contains Cl I. Consequently the Stout `ELECTRON` token assigns a proton approximation to electron collisions; converting it to electron Upsilon would preserve that physical error. The paper lists a constant coefficient for this two-level Cl I transition. Its explicit 30% accuracy statement over 30–3000 K concerns **power-law fits for three-level systems**, so it does not establish a Cl I electron-rate validity interval. The nine constant Stout entries are a tabulation of the one paper value, not nine independently calculated temperatures.

As a proton **conditional** two-level model, the table value would imply `q_12^p(T)=4.7e-8 (2/4) exp[-(882.351 hc/k)/T] cm^3 s^-1` if the proton Maxwellian detailed-balance relation is appropriate. This requires proton density and an explicit proton collision term; it is not a substitute for electron cooling in the present electron-only source-node format. The private diagnostic records the numerical consequences of the mislabeled electron conversion solely to make the error auditable. No Cl I electron rate was found in the local Stout file or the ADS searches performed here; this is a bounded search result, not a proof of absence from the literature.

## Cl V: a coherent two-level electron source is locally invoked

MOCASSIN v2 `data/clv.dat` exists, unlike its listed `data/cai.dat`. Its first line declares five comments, followed by `2 15 0`: two levels, 15 electron-temperature nodes, and source Upsilon rather than direct rate coefficients. Its single `1 2` Upsilon sequence spans 2,000–200,000 K. The same file gives `g=(2,4)`, `DeltaE=1490.8 cm^-1`, and `A_21=2.98e-2 s^-1`. MOCASSIN `source/hydro_mod.f90::makeElements` reads `data/fileNames.dat`, opens the actual file, and sets `lgDataAvailable(17,5)` when present. `source/emission_mod.f90::forLines` then calls `equilibrium` for that ion. The equilibrium routine interpolates Upsilon as a spline in log10 temperature and computes the upward and downward electron coefficients using the proper lower and upper weights. Thus this is an invoked calculation, not merely a listed filename.

The file's dated header calls the collision calculation “Saraph & Storey 1996 ... in prep.” ADS identifies the final publication as [Saraph and Storey (1999), *Astronomy and Astrophysics Supplement Series* **134**, 369, DOI 10.1051/aas:1999437](https://doi.org/10.1051/aas:1999437), Iron Project XXX. The published abstract explicitly identifies the ground-term Cl V electron fine-structure transition and a close-coupling treatment sensitive to resonances. The 15 tabulated values and their exact relationship to the final journal table remain to be verified against its methods and numbers; the publisher PDF and ADS PDF routes returned HTTP 403 and 404, respectively. The MOCASSIN header's original A attribution is unpublished, but independent Cloudy Stout/NIST `cl_5.tp` gives the same `2.98e-2 s^-1`, and `cl_5.nrg` gives `g=(2,4)` with a 1490.5 cm^-1 gap. This supports a coherent **two-level** source candidate within the tabulated 2,000–200,000 K interval. It does not establish accuracy outside the interval or complete higher-level Cl V cooling. CHIANTI `cl_5.easplups` is explicitly an excitation-autoionization parameter file, not a bound-bound Upsilon source. Cloudy Stout `cl_5.coll` has no collision rows; its generic g-bar/default-strength branch is an approximation, not independent Cl V data.

## Ca I: published electron data exist, but direct machine-readable adoption is unresolved

ADS confirms [Samson and Berrington (2001), *Atomic Data and Nuclear Data Tables* **77**, 87–95, DOI 10.1006/adnd.2000.0847](https://doi.org/10.1006/adnd.2000.0847). Its final published abstract describes an LS-coupled, 22-state R-matrix calculation with thermally averaged effective collision strengths from the `4s^2 1S` ground state to the next ten states at 1,000–10,000 K. This provides **ground-state-only term-resolved** electron data, not a complete matrix among excited fine-structure levels. A restricted ground-excitation model would need explicit treatment of the omitted excited-state collisions and a physically justified mapping from LS terms to any fine-structure radiative model; assigning equal or statistical fine-level strengths without a source is not verified here.

ADS also confirms [Zatsarinny, Parker, and Bartschat (2019), *Physical Review A* **99**, 012706, DOI 10.1103/PhysRevA.99.012706](https://doi.org/10.1103/PhysRevA.99.012706). Its final journal abstract describes B-spline R-matrix cross sections for all transitions among 39 physical Ca I states, from threshold to 100 eV, with 444 pseudostates used to test continuum coupling. That is the stronger potential complete Ca I source if the level-resolved numerical cross sections and radiative data can be obtained and Maxwellian integrated. The final paper's machine-readable full cross-section set was not found in this pass; the abstract alone does not provide numerical rows, temperature accuracy, or fine-level mapping. Elsevier's 2001 article route returned HTTP 403. An accepted manuscript for the 2019 article is indexed online, but it was not substituted for the published methods or tables.

The actual local MOCASSIN v2 and v3 `fileNames.dat` list `data/cai.dat`, but both installed data directories lack that file, and `makeElements` sets `lgDataAvailable` false when the open fails; `forLines` therefore skips its Ca I equilibrium calculation. Cloudy Stout `ca_1.nrg` and `ca_1.tp` contain many level energies and radiative rates, while `ca_1.coll` has no electron collision rows and says only “Use gbar.” Cloudy's source code then uses its g-bar or default collision-strength approximation for unprovided transitions. TORUS has `torusdata/CaI.atm` with `CBX` bound-bound entries, but `modelatom_mod.F90::collisionRate` implements `CBB` only for H I, He I, and He II and raises a fatal error for its default atom case; it is not a working Ca I electron-excitation calculation. CHIANTI's local `ca_1.diparams` is direct-ionization data, not excitation. The `Wood_Codes/ionize` routines contain a variable `cl5` in carbon line-cooling expressions, not a chlorine-five atomic dataset; searches of their Fortran and data files identified no invoked Ca I or Cl V atom. MAPPINGS output contains Cl V lines, but inspection here did not identify a standalone Ca I or Cl V level-resolved electron source or establish the exact internal line-data mapping; MAPPINGS output alone is insufficient to transfer a dataset.

## Physical decision points

The subsequent [MAPPINGS invoked-source audit](validation_2026-10-10/mappings_remaining_excitation/README.md)
located actual constant-strength two-level Cl I and Cl V calculations in
MAPPINGS V 5.2.1. Its Cl I ground fine-structure weights are reversed relative
to the identified physical levels. Its legacy resonance file is not invoked
by the inspected reader. These results refine the initial MAPPINGS search;
they do not establish a verified temperature-dependent Cl I electron source.
The constant strengths lack individual source attribution and accuracy ranges.

The [Cl V final-source access investigation](validation_2026-10-10/cl_v_final_source_completion/README.md)
compares ten MOCASSIN nodes with search-index excerpts of the final table,
but the complete final Methods and higher-temperature table remain unread.
Those excerpts are explicitly provisional evidence. The
[Osorio Ca I model investigation](validation_2026-10-10/ca_i_model_source_audit/README.md)
also found no accessible numerical atom tied to that publication.

1. Exclude the Stout Cl I `RATE ELECTRON` row from electron cooling. Investigate a genuine Cl I electron fine-structure calculation; if modeling proton excitation, add a proton collider term with its own density and documented range.
2. Cl V has a narrow, physically coherent two-level candidate. Verify the original final Saraph–Storey table and interpolation validity before production import, then compare full temperature and density cooling profiles with that same two-level source.
3. Ca I requires numerical access to the 2001 table or preferably the 2019 complete cross sections, plus source-consistent radiative data and a declared term/fine-structure mapping. A file listing, generic g-bar estimate, or rates for Ca II and higher calcium stages cannot supply those data.
