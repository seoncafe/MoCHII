# October 11 update to the missing-atom snapshot

The two-missing-atom inventory below predates the approved Mg I installation.
The default `mocassin_v2` selection now contains 55 atoms and 49 cooling fits;
only Ca VII lacks an accepted source atom. Its mixed radiative model and four
nonpositive photon gaps remain unresolved. Production still declares 43 metal
stages and 32 adjacent transitions. See the [current result](../ATOMIC_IONIZATION_STAGE_EXPANSION_RESULTS_2026-10-11.md).

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

# Remaining physical requirements after oxygen cooling installation

The user has revised the first-phase goal to **77 metal stages and 66 adjacent
transitions**, with explicit level, collision-strength, and line-cooling
coverage of the supplied MOCASSIN v2 models. See the
[current measured comparison](phase1_77_mocassin_v2/README.md).
The 64-stage inventories below describe earlier preparation. Requirements
for line atoms absent from v2 must be distinguished from first-phase parity;
their absence does not invalidate physically justified models already installed.
A 165-stage expansion is outside the approved implementation scope.
Terminal ionization must therefore be assessed as a bounded charge-truncation
approximation, with its physical assumptions and validity conditions recorded;
the scope decision alone does not validate that approximation.

The approved expansion is not complete. Current production `ELEMENT` records
declare 43 metal stages and 32 adjacent transitions. Direct file inspection
finds 36 ordinary electron-excitation fits among those stages. The 64-stage,
53-transition rate roster remains staged. Source-specific corrections and
their scoped validation are complete for the newly installed full O II,
O III, Mg III, and Ne I atoms. No numerical calculation controller remains
running from these installations.

## Remaining excitation inputs

| Ion | Established requirement before a full production cooling model |
| --- | --- |
| Ar I | Resolve two extraordinary effective collision strengths in the 35-term ADAS LS source; establish radiation provenance and data-use permission. Its metastable LS term is a model limitation, not automatically a source error. |
| Ca I | Obtain the published 39-LS-state electron data and an identified, compatible spontaneous-rate model. Stout's 780 fine-structure states cannot be joined by index to 39 LS states. |
| Ca III | Resolve source-use permission and a coherent physical-energy/radiation model for the private ADAS distorted-wave data; do not treat placeholder A values as physical spontaneous probabilities. |
| Cl I | Obtain electron-impact data. A proton substitute in the inspected source does not establish neutral electron excitation. |
| Cl V | Verify the final source Methods and supported temperature range of the two-state collision model. Numerical reproduction of MOCASSIN's extended table does not establish source accuracy up to 100,000 K. |
| Fe I | Obtain the complete final Conroy energy/collision tables and establish every retained collision state's identity in the chosen radiation source; integer indices and dominant configuration labels alone are insufficient. |
| Mg I | Resolve the CHIANTI negative infinite-temperature endpoints and low-temperature conversion discrepancy. Published LS strengths extend only to 20,000 K, and cannot be inserted directly into a fine-structure atom. |

The full [Si I original-radiation audit](validation_2026-10-10/si_i_electric_dipole_source_audit/README.md)
has now recovered the official NIST original calculation and identified all
185 represented bound E1 pairs, including 15 doubled coefficients and the
14/15 inversion. Root independently reproduced the audit and verified the
original E1, E2, M1, and M2 strength-to-rate conventions. The corrected coherent
35-state model has completed full numerical validation and is installed;
see the [installation review](validation_2026-10-10/si_i_cooling_products/INSTALLATION_REVIEW.md).
This resolves the identified radiative defects and the missing cooling file,
but does not verify absolute atomic accuracy or the unread 2005 Methods.
The remaining seven missing excitation inputs are listed above.

A new [Cl I literature investigation](validation_2026-10-10/cl_i_electron_source_completion/README.md)
identifies Wang et al. (2013) as the appropriate neutral fine-structure
calculation. Its official accepted manuscript is preserved in the central
references directory, explicitly distinguished from the inaccessible final
article. Numerical cross sections are still unavailable. The subsequent
rate-fit interval starts near 5,802 K and does not justify a continuation
to 1,000 K; its final table's coefficient signs remain unverified.

The root session inspected the exact Si I/Ar I source rows, independently
remeasured the Ar I inventory and Fe I source counts and published energy
residuals, and recorded its review in
[the evidence record](validation_2026-10-10/remaining_source_identity_main_review.json).
Detailed direct-source measurements are in the
[Ar I audit](validation_2026-10-10/ar_i_ls_source_identity/README.md),
[Fe I audit](validation_2026-10-10/fe_i_radiative_identity/README.md),
[Si I audit](validation_2026-10-10/si_i_radiative_identity/README.md),
[Ne I/Mg I assessment](source_validity_ne_i_mg_i.md), and
[chlorine/calcium assessment](remaining_chlorine_calcium_excitation_sources.md).
The subsequent [Si I/Mg I source investigation](validation_2026-10-10/si_mg_source_completion/README.md)
was an earlier checkpoint before Si I installation. Its unresolved absolute
Si I accuracy caveat remains valid; Mg I's author-provided matrices still stop
at 20,000 K. The
[Fe I/Ca I access record](validation_2026-10-10/fe_ca_source_completion/README.md)
distinguishes unavailable downloads from absence of data.
These requirements cannot be satisfied by clipping a source node, inferring
a missing A as zero, or guessing a corrected exponent. New source-supported
repairs remain within the approved work; choosing an additional physical
approximation requires an explicit validity statement and scientific review.

## Ionization and thermal closure

The direct-photon boundary at 100 eV closes outer photoionization for the
proposed terminal stages. It does not close electron-impact ionization at
finite temperature. The existing terminal diagnostic is a restricted
electron-only sensitivity calculation, not an accepted truncation criterion.
Private diagnostics may evaluate omitted charge states without changing the
approved production roster. Adding one diagnostic stage for every element
would not establish closure against all higher charge states. The accepted
64-stage boundary must follow a defined physical model and source-valid
reverse rates; no new temperature limit or error tolerance has been adopted.
An [original-node reverse-rate diagnostic](validation_2026-10-10/terminal_native_node_bounds/README.md)
now separates total recombination fit discrepancies from the analytic bound
for exactly one omitted state in an adjacent electron-only chain. Its sparse
original nodes and excluded processes prevent its use as a complete
production truncation criterion.

The nine energetically reachable nonadjacent photoionization products fit
the proposed 64-stage roster. Their production charge balance and thermal
accounting remain incomplete. Published Kaastra--Mewe electron-count
probabilities, mean Auger energy, and fluorescence yields do not provide a
joint final-charge/electron/photon/excitation distribution. Combining those
means with the current PH1 edges and successive ionization energies can
violate event energy conservation. A tested standalone nonadjacent solver is
neither linked nor invoked by production.
The [charge and energy bounds](validation_2026-10-10/auger_energy_completion/README.md)
now explicitly locate three unresolved shell channels at 100 eV. The other
six support a charge-marginal diagnostic lookup, not a complete thermal model.

The concrete required input information and acceptance checks are recorded
in [Auger implementation requirements](auger_implementation_requirements.md)
and [terminal reverse-rate source validity](terminal_reverse_rate_source_validity.md).
New public atomic layouts, a different ion roster, or changed meanings of
thermal outputs are material changes. They are not implied by routine
continuation of the approved 64-stage plan. No such change has been made.

## Scope relative to MOCASSIN v2

Shared process classes include adjacent photoionization/recombination,
available charge exchange, and electron-excitation/spontaneous line cooling.
Larger atoms and corrected coefficients improve these same processes.
The inspected local v2 data do not invoke Ne I collision cooling or Mg III
line cooling; MoCHII now has validated numerical models for both.

Current MoCHII also invokes metal electron-impact ionization and sums active
inner-shell cross sections into photoionization rates. The inspected v2 metal
balance zeroes its collisional-ionization coefficient and uses outer-shell
photoionization rates. MoCHII's broader shell rate is still assigned to an
adjacent stage, so this does not establish complete inner-shell ionization.
The adopted high-charge Langevin capture upper limit and cold constant-Upsilon
continuation remain explicit approximations. Full Auger and terminal closure
must be reported separately as unfinished physics.

No input-matched, converged MOCASSIN/MoCHII nebular comparison has been run.
Such a comparison must match spectra, boundaries, density and geometry,
species/reactions, and enabled energy terms, and inspect complete spatial
profiles. The source and atom validation completed here cannot substitute
for that comparison or establish whole-code parity. See the
[invoked-path scope report](mocassin_v2_physics_scope.md).

## Resumed 64-stage implementation: O IV source review

The user has resumed the 64-metal-stage objective. O IV full-atom numerical
validation is the next prepared calculation, but root first detected a
source-frequency inconsistency requiring radiative provenance review. The
current selector returns 147 states, 2,319 radiative rows, and 10,731 collision
pairs. For 343 higher E1 pairs, stored gf and A do not imply a frequency
within 1% of any available observed/theoretical/wavelength frequency. No
correction is inferred from this diagnostic alone. See the
[direct source review](validation_2026-10-10/o_iv_main_radiative_review/README.md).
Production remains at 43 metal stages and 36 ordinary cooling fits.

## Latest O IV/O V checkpoint

The complete O IV source-data cooling model is now installed after qualified root adoption and actual installed reader checks. The earlier requirement to resolve every unused gf/A discrepancy before installing a collision/radiation model was too strong: the actual n-level path uses only identified energies and imported A. No incorrect A has been established; detailed source replacements and atomic accuracy remain explicitly uncertain. This uncertainty is not certification of the source, but it does not itself demonstrate a physically wrong cooling equation. See the [installation review](validation_2026-10-10/o_iv_full_completion/installation/README.md).

O V completed its full-bound numerical validation, with physical source review still pending. The latest [21-ion inventory](validation_2026-10-10/upper_stage_readiness/README.md) is a historical pre-O-IV-installation snapshot; its roster difference still has 21 ions, while O IV is now the first with installed data. The other 20 require their remaining source-specific checks and installation. The seven original-roster missing excitation sources, nonadjacent Auger energy/charge accounting, terminal-stage validity, and final 64-stage activation remain outstanding. Ca VI's CHIANTI collision footer gives logT=4.5--6.5; Fe V's stated fit interval ends at 50,000 K. Neither can be treated as a verified 1,000--100,000 K model merely because an extrapolated private fit was generated.

The latest oxygen checkpoint supersedes the preceding O V pending-installation statement: both O IV and O V are adopted and installed as explicitly qualified imported-rate cooling models. The [direct inventory](validation_2026-10-10/upper_stage_readiness/post_oxygen_installation.json) leaves 19 of 21 additional-ion atoms/fits uninstalled. Atomic source accuracy remains unverified, and numerical passes do not establish a coupled thermal/charge model. Production still has 43 stages/32 transitions; the original seven missing excitation sources and Auger/terminal closure requirements remain outstanding.
