# Additional MAPPINGS excitation-source inspection

The local MAPPINGS V 5.2.1 tree contains invoked legacy two-level data for
Mg I, Si I, Cl I, and Cl V. This is a narrower finding than a complete
multi-level atom. Its `MULTIDAT.txt` reader was followed through
`mapinit.f::readmultilevel`; parsing all 174 atom records found none of the
eight remaining target ions. `RESONDAT2.txt` also lists neutral Cl/Ar rows,
but its read block and the old `reson` implementation are commented out.
Those particular file rows do not establish an invoked excitation source.

`mapinit.f::read2level` actually reads `INTERDAT2.txt`, stores its weights and
constant collision strength, and `cool.f` calls `inter`. That routine forms
`q12 = rka Upsilon exp(-E/kT)/(g1 sqrt(T))` and
`q21 = rka Upsilon/(g2 sqrt(T))`. The Maxwellian option and an isolated
two-level balance are used in this audit; the optional kappa correction is
excluded. The data contain no individual source references or temperature
accuracy bounds for these constant strengths.

## Confirmed Cl I statistical-weight defect

The invoked Cl I ground fine-structure row has energy `1.752743e-13 erg`,
`A21=0.0123 s^-1`, constant strength `3.120`, and weights `(g1,g2)=(2,4)`.
The inspected NIST-attributed Stout level file instead identifies the ground
`3s2.3p5 2P*3/2` with `g=4`, and its `882.351 cm^-1` upper
`2P*1/2` with `g=2`. The MAPPINGS weights are reversed for this physical
ground fine-structure transition. Its actual implementation then uses those
weights in the excitation/deexcitation coefficients. The external source is
preserved for comparison; no external file was corrected.

Holding the unknown source strength and energy fixed isolates this defect.
The source excitation coefficient is twice the value with physical weights,
while its downward coefficient is half that value. Direct evaluation of the
two-level population and cooling over 501 temperatures from 1,000 to
100,000 K and four electron densities (`1e-6,100,1e6,1e10 cm^-3`)
gives source/physical-weight cooling ratios from `2.0000000000000435` to
`2.92068442260722`. These are consequences of substituting weights in the
same conditional model, not an estimate of true Cl I electron cooling.

| Input condition | Source calculation | Physical-weight diagnostic |
| --- | --- | --- |
| Temperature grid | 501 points, 1,000--100,000 K | Identical |
| Electron density | Four values listed above | Identical |
| Spectrum and boundaries | No photons; local two-level model | Identical |
| Other colliders and charge reactions | Excluded | Identical |
| Transition energy and A | Original MAPPINGS row | Identical |
| Effective strength | Original constant 3.120 | Identical; accuracy unverified |
| Statistical weights | 2,4 | 4,2 from identified physical states |

The [full temperature profiles](cl_i_weight_profiles.pdf) and
[saved numerical arrays](cl_i_weight_profiles.npz) cover all four densities.
This is a local rate/population comparison; it has no spatial temperature,
velocity, or density solution and does not establish whole-code agreement.
The [source inventory and hashes](result.json) and
[reproducible audit](audit_mappings_sources.py) identify every inspected file.
No existing product was regenerated and no production MoCHII data changed.

These MAPPINGS constants therefore cannot be imported as verified replacement
data for the missing full-range cooling atoms. In particular, the existing
Cl I source search must distinguish this actual electron-labeled calculation
from both a physically validated electron strength and the separate
Hollenbach--McKee proton substitution.
