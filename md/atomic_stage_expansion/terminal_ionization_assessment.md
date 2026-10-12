# Terminal electron-impact ionization assessment

## Physical judgment

The 64-stage roster is a defensible **photon-threshold staging choice** below
100 eV, but it is not a physically closed ionization network at unrestricted
temperature or in the presence of inner-shell multiple ionization. Electron
impact can ionize every listed terminal ion at finite temperature. The current
runtime has no outgoing transition from those ions. The 100 eV photon cutoff
alone cannot justify setting those rates to zero. More stages, or a stated and
enforced physical validity range for a deliberate truncation approximation,
are required before calling the network generally complete. The limited
calculation below does not decide the allowable range because the terminal
hydrogen charge-transfer rates and nonadjacent Auger products remain unresolved.

## What was evaluated

The source diagnostic reads the staged 64-ion element files under
`validation_2026-10-10/langevin64_exit_heat/products/atomic`, not the
currently installed 43-stage `data/atomic` files. It uses their 53 internal
Voronov CI and RR+DR transitions, and adds one outgoing transition per element
from the local Voronov `cfit.f`, Badnell M=1 RR/DR, or CHIANTI `.rrparams` and
`.drparams` APIs. The next-stage ionization energy is from the identified NIST
snapshot used by `make_element_data.py`. The 53 internal CI rates receive the
exact stage constants in `recomb_mod::ci_dere_ratio` for `dere_hybrid`. That
function has no terminal-stage entry for any of the 11 elements, so it would
return unity at all 11 additional transitions; this is an implementation
return value, not a validated Dere correction. Both truncated and extended
calculations use this same convention.

All 801 sampled temperatures from 1,000 to 100,000 K use stationary adjacent
electron processes, zero photoionization, fully ionized hydrogen, and no H
charge transfer. Electron density cancels from the stage ratios. The energy
diagnostic uses the physical ionization energy of each transition, as
`species_mod::metal_ionization_cooling` does for its CI term. The plots show
every ion stage across the full temperature range. There is no spatial domain,
gravitational potential, gas velocity, or density profile in this rate test;
no claim about a nebular solution follows from it.

| Condition | Truncated calculation | One-stage extension |
|---|---|---|
| Temperature | 1,000–100,000 K, 801 logarithmic points | Identical |
| Radiation | No photons or diffuse field | Identical |
| Hydrogen | Fully ionized; no H charge transfer | Identical |
| Electron density | Positive, cancels from adjacent equilibrium | Identical |
| Atomic rates | Staged 64-ion files, internal HYBRID factors | Same rates plus local source fits at each terminal boundary; terminal factor 1 |
| Thermal quantity | Internal CI ionization-energy loss | Same sum including terminal CI loss |
| Auger ionization | Excluded | Excluded |

The measured maxima over this temperature interval are:

| Element | Next-stage maximum fraction | Maximum added CI energy fraction |
|---|---:|---:|
| C | 3.152e-20 | 4.129e-20 |
| N | 4.067e-33 | 1.017e-32 |
| O | 5.541e-8 | 4.461e-7 |
| Ne | 1.911e-8 | 2.276e-7 |
| S | 2.635e-18 | 1.609e-18 |
| Ar | 9.814e-13 | 8.855e-12 |
| Mg | 2.377e-4 | 6.152e-3 |
| Fe | 4.141e-9 | 1.236e-7 |
| Si | 9.521e-7 | 1.064e-6 |
| Cl | 1.930e-11 | 4.483e-12 |
| Ca | 3.974e-9 | 3.611e-8 |

These are measured source-fit diagnostics, not an acceptance tolerance. Mg IV
has the largest added population and energy loss in the specified interval.
At higher temperatures these omitted transitions grow; this finite interval
cannot establish an unrestricted maximum. The internal-ion populations were
renormalized after adding each stage, so the plotted difference includes the
effect on all stages.

[Ion-fraction profiles](validation_2026-10-10/terminal_ionization/ion_fractions.pdf),
[CI energy profiles](validation_2026-10-10/terminal_ionization/ionization_energy.pdf),
[all 801-point fraction and log-fraction arrays](validation_2026-10-10/terminal_ionization/full_profiles.npz),
and [measured results with source hashes](validation_2026-10-10/terminal_ionization/result.json)
are produced by `tools/fitting/validate_terminal_ionization.py`. Both PDF
pages were visually inspected for legibility. The fraction figure displays
logarithms down to -100; the array file retains the smaller values without
the display limit. The energy figure plots only strictly positive values on
its logarithmic axis, with no artificial positive floor. No production
objects or installed atomic files were changed or rebuilt.

## Source and implementation checks

The actual runtime calls `voronov_ci * ci_dere_ratio` in
`species_mod::species_ne_prepare`, `species_fractions`, and
`metal_ionization_cooling` only for `it = 1 .. nstage-1`. The cooling term is
`n_e n_i C_i E_th` for those transitions. `make_element_data.py` also requires
an explicit `CX_HIGH` and `CX_HEAT` entry before adding a boundary to an
element file; none of the 11 terminal boundaries has such an entry. The
polarization-capture rates for selected *internal* high-charge transitions
cannot be extended to terminal transitions as measured charge-transfer rates:
their electron-transfer probability, product branching, and reliable
temperature interval are not established. A zero terminal CX rate is therefore
not a verified physical value.

The Fe VIII recombination source has been checked against the final journal
paper: Schmidt et al. (2008), *A&A* **492**, 265–275,
DOI [10.1051/0004-6361:200810834](https://doi.org/10.1051/0004-6361:200810834),
[author-hosted final journal PDF](https://amdpp.phys.strath.ac.uk/tamoc/PUB/2008/aa_492_265.pdf).
The nine CHIANTI `fe_8.drparams` pairs match its Table 3 **experimental Fe
VIII** column, including the low-energy resonance contribution; the six
`fe_8.rrparams` values match Table 4's theoretical Fe VIII RR fit. All 18 DR
and six RR numeric coefficients were compared directly with text extracted
using `pdftotext -layout` from the final journal PDF, and all matched. The
copy is stored at
`/home/kiseon/MoCHII/references/Schmidt_2008_AA_492_265.pdf` with SHA-256
`002bbbaaea2420a1ad67d1142a1095f940d06cd38ea132d1a1c4bdf2820f1cf6`.
Here Fe
VIII is the recombining ion in Fe VIII to Fe VII, so it is the correct reverse
of terminal Fe VII to Fe VIII ionization. The paper's experimentally derived
plasma DR+TR rates cover electron thermal energies 0.2–1000 eV, or about
2,321–11.6 million K. Its analytic Table 3 fit is reported accurate within
1% to the derived curve over 1 meV–100 keV. Thus the 1,000 K end of this
diagnostic uses an analytic extension below the paper's stated experimental
plasma range, while 100,000 K is inside it. This distinction does not change
the need to retain the outgoing physical process.

For the other ten terminal recombining ions, local M=1 records were located
and evaluated, but their individual published fitting ranges have not been
verified here. Voronov (1997) source coefficients were read directly from
the authors' local `cfit.f`; a transition-specific error bound for each
terminal ion across this temperature interval was not established. The source
fit is therefore suitable for sensitivity analysis, not proof of a certified
truncation error. The ADS API lookup succeeded using authorized network access
and returned bibcode `2008A&A...492..265S`, DOI
`10.1051/0004-6361:200810834`, publisher EDP, and publisher PDF/HTML links.
The ADS publisher PDF link and DOI both reached the A&A publisher but returned
HTTP 403. The author-hosted PDF is the paginated final A&A article, not an
unpublished preprint; it was downloaded and checked locally as described
above. The publisher-hosted copy was not obtained.

The comparison with local MOCASSIN v2 is about computed scope, not matching
physics. In `mocassin-mocassin.2.02.73.2`, `set_input_mod.f90:82` sets
`nstages=7`. `update_mod.f90:1459` limits stage loops to
`min(elem,nstages-1)`. Its `collIon` is calculated for hydrogen at line 1378,
then explicitly set to zero for every metal at line 1485 before the
ionization-ratio equations at lines 1490–1496. It also sets charge-exchange
rates to zero outside 6,000–50,000 K at line 1466. Consequently, v2's
metal-ion results cannot validate MoCHII's metal CI treatment or the present
64-stage thermal boundary. No behavior from MOCASSIN v3 is attributed to v2.

## Decision required for a production change

Keep the current 64-stage files as an explicitly provisional photon-limited
set. A production correction needs a chosen physical operating domain and a
representation of every ion reachable by the included CI, photoionization,
and Auger processes. If the domain remains thermally open, add higher stages
with their reverse recombination, cooling, and output treatment. If a finite
domain is chosen, specify it in the user interface and enforce it at runtime;
also establish a physically justified bound on terminal charge transfer and
nonadjacent channels. Changing the public stage roster, file layout, or user
temperature limits needs prior design approval. None was changed by this
assessment.
