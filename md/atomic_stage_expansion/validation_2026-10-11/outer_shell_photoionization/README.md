# Consistent outer-shell photoionization component

A private `outer_shell_photoionization` routine is implemented in the production
species module. Below the supplied PHINNER boundary it uses the existing PH2
outer-shell fit. Above that boundary it uses the unique PH1 row assigned the
physical ionization threshold by the generator. Both the current 32-transition
files and staged 66-transition files have exactly one such row for every
transition. Compressed row position is not interpreted as a subshell index.
The heat coefficient is sigma times (1 - ionization threshold / photon energy).
This is an explicit restricted absorption model, not a complete vacancy cascade.

## Actual code and numerical checks

The production cross-section routine calls this component below PHINNER,
retaining the former arithmetic. The private component is not yet selected
above that boundary in production: current occupied-shell sums and their
known product/energy limitations remain. No new public input was introduced.
The diagnostic compilation exposes the private routine only in its own module
object, retaining the rest of the production physics and linking the actual
production object list. The full production target also rebuilt successfully.

The current 32 transitions passed 576 samples; the staged 66 transitions
passed 1,188, for 1,764 total. Independent Python PH1/PH2 expressions match
within 4.440892098500626e-16 relative. Below PHINNER, all 379 current-roster
samples and 887 expanded-roster samples retain exactly the previous cross
section and heat values. Additional absorption from other shells appears in
165 and 243 sampled cases, respectively. These counts are samples, not
physical event fractions or integrated spectral effects.

| Input condition | Linked diagnostic | Independent expression |
| --- | --- | --- |
| Fits and physical outer threshold | Existing production or staged records | Identical records |
| Photon energy | 15 fixed energies from 10 to 10,000 eV and 0.99/1/1.01 times each threshold | Same values |
| Outer PH2/PH1 boundary | Stored PHINNER | Same boundary |
| Additional opacity, radiation field, thermal equilibrium | Not solved | Not solved |
| Ion fractions and abundances | Not used in coefficient comparison | Not used |

This establishes implementation of the specified outer absorption and heating
formulas. It does not certify publication-level atomic accuracy, a continuous
energy error bound, or a spatial plasma solution. Sulfur's existing product
atom and fit are supplied as setup dependencies; their tables were not
recomputed. [Root review](root_review.json) verifies sources and outputs.
The finite controller completed both steps; it is no longer running.

## MOCASSIN v2 opacity and rate distinction

The actual v2 `ionization_mod.putOpacity` loops over subshells and accumulates
absorption cross sections. Its `update_mod` cell ionization and heat integrals
select `outShell` alone. `ph_mod.phFitEl` supplies a single selected subshell,
not a hidden sum over all occupied shells. Consequently inner-vacancy opacity
is not paired with a cell ionization/product/heating calculation for that
vacancy. Simply retaining full shell opacity while disabling all its rate and
energy consequences would preserve that inconsistency.

There is a further reference restriction: `hydro_mod.limitShell` limits the
4s table above the lower-pointer edge of 3d, or 3p when the 3d lower pointer
is below three. The 3s/3p/3d tables instead stop below the K edge, and 1s/2s/2p
extend to the configured maximum. Cell rate and heat paths use those stored
windows, with different endpoint pointer conventions. The physical PH1 data
for Fe I/II have 3d edges 14.7/21.93 eV and 3p edges 66.0/76.17 eV; Ca I/II
have empty-3d placeholder edges and 3p edges 34.43/40.9 eV. Actual window
selection depends on photon-bin input, not just the atomic edge.
[The direct table and source audit](mocassin_outer_windows.json) retains these
values and hashes. No MOCASSIN execution or matched bin grid is claimed.

The physically retained outer-shell tail in the new component is therefore
not literal equality with every v2 stored window. Dropping a valid valence
absorption tail solely to copy that window is not the acceptance criterion.
A default restriction must consistently select opacity, ionization, and heat.
Full inner-shell selection additionally requires accepted charge products,
fluorescence, and electron energy; the existing Auger limitations remain.

## Remaining integration

A reviewable input design must state that the same shell selection applies to
all three consumers. Any enabled full-shell calculation must reject or resolve
known unsupported vacancy/product regimes, rather than silently keep an
adjacent product where Auger emission changes the charge. Input choices and
remaining physical validity conditions must be settled before changing the
public parameter type. This report does not claim a completed selector or
correction of the current above-inner-edge production calculation.
