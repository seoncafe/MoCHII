# Qualified upper-temperature endpoint model: candidate review

## Physical model and acceptance state

Hold effective collision strengths at their original final source nodes for
C I, S III, and Ar VI over the proposed 1000-100000 K numerical model interval.
Only Upsilon is clamped; excitation exponentials, de-excitation factors, and
statistical equilibrium retain actual gas temperature. This follows the
original v2 upper-endpoint assumption while keeping MoCHII's corrected
positive log-PCHIP interpolation within the supplied table. It is not exact
v2 cubic-spline equality or verified high-temperature atomic accuracy.
The model has not been accepted or installed. Production inputs, atomic files,
and physical expressions remain unchanged.

The evaluation uses copied diagnostic Fortran modules with explicit model
bounds at the point of use. No new public parameter or atomic input format
was introduced. The existing source atoms are read unchanged. Above-model
values such as 100000.1 K are rejected. A representational rounding value
100000.00000000001 K is accepted by the same log-T bound used for registration;
this does not authorize a physical extension beyond the proposed model scope.

## Preserving the original interpolation

Simply appending a constant synthetic hot node changes the original PCHIP
endpoint derivative and therefore changes strengths inside the old table.
This was measured on 201 temperatures within each final original interval.
The private evaluator instead clamps only evaluations above the original
upper node and preserves the original polynomials. All 1212 in-range
cooling/population profile comparisons were exactly identical to the unchanged
native evaluator's printed values; no old source nodes were modified.

| Element and stage | Original upper node (K) | Levels | Maximum in-range strength change with appended node | Maximum registered/direct cooling difference |
| --- | ---: | ---: | ---: | ---: |
| C 1 | 20000 | 6 | 2.829553% | 0.134281% |
| S 3 | 20000 | 5 | 1.271571% | 0.087795% |
| AR 6 | 40000 | 2 | 0.100592% | 0.159297% |

The appended-node differences are collision-strength differences, not a
measured change in total nebular cooling. They justify retaining the original
interpolant rather than silently changing an already validated temperature
interval. Mg I's separately prepared constant model has no nonzero source
slope and does not have this particular interpolation issue.

## Completed actual numerical checks

The direct private evaluator passed 2412 combined samples: 201 temperatures
for each of three ions and four electron densities. The low/hot grids join
at each original upper node. Native cooling agrees with the independent
statistical-equilibrium matrix to at most 2.9976021664879227e-15; maximum
absolute population difference is 1.4016565685892601e-15. The matrix uses the
same levels, positive A values, source PCHIP strengths, actual temperature,
and ne=1e-6, 100, 1e6, and 1e10 cm^-3. No photoexcitation, trapping, spatial
flow, gravitational potential, spectrum, or boundary-condition comparison is
included. This validates the stated approximation, not atomic accuracy.

The actual registered private cooling tables passed 120 direct/apply tests,
including regular temperature nodes and intermediate temperature/density
locations. The maximum interpolation discrepancy is 0.159297%, for Ar VI,
below the 0.5% criterion. Diagnostic access measured 301 temperature nodes
and 161 density nodes, or 48461 cells, for each table. Public production
access was not added. The generic validation wrapper now reports measured
shape only when the driver provides it; older outputs remain explicitly
unmeasured and are not regenerated.

Nine further guard probes checked the nominal upper limit, one rounded
representation above it, and an actual unsupported temperature for all three
ions. Root verified the recorded sources and outputs and visually inspected
all three [full-profile pages](../upper_temperature_models/profiles/endpoint_profiles.pdf).
No clipping or overlapping labels was found. Populations below 1e-30 are
shown at that display floor; full unmodified numbers remain in the saved
profiles. See [direct model result](../upper_temperature_models/profiles/result.json),
[registered result](cooling/result.json), [guard probes](guard_probes.json),
and [root review](root_review.json).

## Preserved failure and corrected continuation

The original controller completed direct profiles, then failed registration:
its strict T comparison rejected exp(log(100000)) rounding. That controller
remains failed. Exact analysis source, completed profiles, and failed products
are preserved. The continuation verified all original sources and profiles,
then compiled a copied module using the matching log-T bound. It ran only the
unfinished registration checks and completed. The profile matrix, original
source interpolation, and accepted nominal sample temperatures did not change;
no uncertain partial output was reused.

No hot choice is inferred from the accepted cold approximation, the explicit
`metal_freefree=.true.` exception, or automatic task continuation. A future
public source-model selection must coordinate this rule with atomic source
bounds, thermal cooling, populations, charge transfer, and line emission.
