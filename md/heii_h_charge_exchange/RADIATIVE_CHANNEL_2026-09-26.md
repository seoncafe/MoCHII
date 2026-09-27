# He+ + H radiative charge transfer: GJ07 R26 or SLD98 (19)?

Date: 2026-09-26. Script, digitized data and figure:
`calculations/radiative_channel_comparison.py`,
`calculations/courtney2021_fig3_rct_digitized.txt`,
`calculations/radiative_channel_comparison.png`.

## The two fits

| Fit | Expression | Source it cites |
|---|---|---|
| Glover & Jappsen (2007) Table 1, R26 | 1.25e-15 (T/300)^0.25 | Zygelman et al. (1989) |
| Stancil, Lepp & Dalgarno (1998) Table 1, reaction (19) | 1.20e-15 (T/300)^0.25 | "Zygelman et al. 1989, multiplied by 0.25 to account for approach probability factor (see Stancil & Zygelman 1996)" |

Both are power laws through the same data, Zygelman et al. (1989) Table I
"Direct" column times the singlet approach factor 1/4, which is tabulated at
1-1000 K only (5.36, 4.21, 4.41, 4.50, 4.83, 5.99e-15 cm^3 s^-1 at 1, 10,
100, 200, 400, 1000 K, before the 1/4). They differ by a constant 4%, a
choice of normalization, and above 1000 K both are extrapolations of the
same T^0.25. MoCHII used GJ07 until 2026-09-26; EXHALE uses SLD98.

## The independent check: Courtney et al. (2021)

Courtney et al. (2021, ApJ 919, 70; read) recompute radiative charge
transfer (their reaction 9) from 1 K to 1e6 K, with the approach probability
of their eqs. 4-6, and state that the smooth background agrees with
Zygelman et al. (1989). Their rate coefficient is given only as Figure 3
(LTE and NLTE zero-density limit, with and without stimulated processes);
the tabulated data are on the UGA Molecular Opacity Project website, which
could not be reached on 2026-09-26. The NLTE zero-density curve without
stimulated processes (red solid; stimulated emission matters only above
~5e4 K in a strong field) was digitized from the published figure, ~3% per
point.

| T [K] | Courtney (digitized) | GJ07 / Courtney | SLD98 / Courtney |
|---|---|---|---|
| 1 000 | 1.62e-15 | 1.040 | 0.999 |
| 3 000 | 2.13e-15 | 1.043 | 1.002 |
| 5 000 | 2.46e-15 | 1.028 | 0.987 |
| 8 000 | 2.72e-15 | 1.046 | 1.004 |
| 10 000 | 2.82e-15 | 1.064 | 1.021 |
| 15 000 | 3.00e-15 | 1.110 | 1.065 |
| 20 000 | 3.15e-15 | 1.132 | 1.087 |
| 30 000 | 3.35e-15 | 1.180 | 1.133 |

Below 1000 K both fits fall below Courtney, to 0.53-0.55x at 100 K and
0.36-0.37x at 10 K.

## Assessment

- **SLD98 is the better of the two.** From 1000 K to 1e4 K it matches the
  modern calculation to within the digitization accuracy (-1.3% to +2.1%);
  GJ07 is 3-6% high there. The two fits share their source, so this is only
  a normalization choice, and SLD98's is the one the later calculation
  supports.
- **Neither is right outside 1e3-1e4 K.** Above 1.5e4 K both overshoot
  (SLD98 +7% at 1.5e4 K to +13% at 3e4 K), because the T^0.25 extrapolation
  keeps rising while the calculation flattens. Below 1000 K both are low by
  up to 2-3x.
- **Better than either:** a fit to Courtney et al.'s tabulated NLTE
  zero-density RCT rate over 1e2-1e5 K. That needs the published data file,
  which cannot be obtained (2026-09-27): the UGA Molecular Opacity Project
  site no longer exists; ExoMol ExoPhoto (exomol.com/exophoto) carries only
  the UGAMOP HeH+ photodissociation cross sections (Miyake, Gay & Stancil
  2011); LAMDA carries excitation data only; UMIST RATE22 row 502 is the
  SLD98 fit itself (1.20e-15 (T/300)^0.25, 10-41000 K, DOI 10.1086/306473).

## Decision (2026-09-26)

Both codes use SLD98 (19), 1.20e-15 (T/300)^0.25. The EXHALE session
reached the same conclusion independently from the published pages
(recorded at its charge_exchange.f90, note B2), and MoCHII switched
(`ion_balance_mod::cx_HeII_H`; gate `check_helium_charge_exchange.f90`).
Both codes state the caveats: Zygelman's table stops at 1000 K, and the
T^0.25 above it overshoots the Courtney calculation by up to 13% at 3e4 K. Nothing in MoCHII's non-radiative channel or detailed
balance changes: the radiative channel has no collisional reverse.
