# Mg I constant effective collision strength candidate

## Physical choice awaiting acceptance

The supplied MOCASSIN v2 Mg I model has one effective-collision-strength node
at 10000 K. Its invoked one-node calculation holds that value constant while
using the actual temperature in the Maxwellian excitation/de-excitation rates.
The prepared five-level candidate applies that explicit approximation over
1000-100000 K. The upper-temperature extension is a model, not new measured
or calculated collision data. Its physical accuracy above 10000 K has not
been established. The previously accepted cold approximation alone does not
authorize this hot extension. Production remains unchanged.

Acceptance would choose the supplied v2 constant-strength model through
100000 K, with an explicit validity qualification and rejection above the
model ceiling. Otherwise Mg I installation needs suitable temperature-dependent
data; current evidence does not support silently extending a finite table.

## Numerical evidence and limitations

The finite controller completed 404 new hot samples at 101 temperatures from
10000 to 100000 K and four electron densities: 1e-6, 100, 1e6, and 1e10 cm^-3.
Native cooling agrees with the independent statistical-equilibrium matrix to
1.1690648449302898e-13; maximum population difference is
3.3306690738754696e-16. It rejected 100000.1 K. Reusing 404 verified unchanged
cold samples gives 804 unique combined samples, including a shared 10000 K
node. This numerical agreement validates the same approximation, not its
atomic accuracy. Both calculations use identical levels, radiative data,
collision strengths, temperatures, and imposed densities. There is no spatial
flow, spectrum, boundary-condition comparison, or independent nebular model.

Root verified recorded sources and artifacts. The full cooling and population
curves were visually inspected without clipping or overlap. See
[results](products/result.json), [full profiles](products/constant_excitation_profiles.pdf),
and [root review](root_review.json). No table was installed and no hot physical
choice is inferred from the user's `metal_freefree=.true.` exception.
