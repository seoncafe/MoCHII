# Hydrogen charge-transfer boundaries and adopted approximation

On October 10, the user adopted the reduced-mass Langevin upper-limit model for
the ten q >= 5 target boundaries. The generator and staged native consumers
are validated; see [adoption and measured checks](langevin_adoption_validation.md).
The species-specific physical transfer probability remains unverified.

The lower-charge physical decision is to use the published Kingdon and Ferland (KF96)
total recombination fits for N 4, Cl 3–4, Ca 3–4, and Fe 4 where selective
staging needs them. Their corresponding Kingdon and Ferland (KF99) Table 1
mean kinetic heat values are available. Fe 4 is a qualified estimate: KF96
marks its Landau–Zener result with the warning that radiative charge transfer
may be faster. None of these total rates establishes a state-resolved reverse
reaction. Charges five and above have a published Langevin capture-rate upper
limit, but no measured species-specific electron-transfer rate in the inspected
sources.

The reaction notation `X q` means `X^(q+) + H(1s) -> X^((q-1)+) + H+`.
For the rows below, `k(T) = 10^-9 a T4^b [1+c exp(d T4)] cm^3 s^-1`,
`T4=T/(10,000 K)`. Values at 10,000 K are direct evaluations of those
coefficients, rather than quoted paper values. Positive mean heat adds kinetic
energy to gas in the recombination direction, assuming ground-state reactants
and low-density radiative relaxation. The KF99 calculation excludes local
absorption of cascade photons and collisional redistribution of products.

| Boundary | KF96 source | `(a,b,c,d)` | Supported temperature evidence | `k(10,000 K)` (cm^3 s^-1) | KF99 mean heat (eV) | Qualification |
|---|---|---|---|---:|---:|---|
| N 4 | Table 1, ref. 11 | `(3.28, 0.52, -0.52, -0.19)` | 1,000–300,000 K printed fit range | 1.8695385e-9 | +11.00 | Cloudy has a changed, unattributed fit. |
| Cl 3 | Table 2, new calculation | `(1.88, 0.32, 1.77, -5.70)` | Calculated at 1,000, 3,162, 10,000, and 31,623 K | 1.89113403e-9 | +8.00 | Landau–Zener total. |
| Cl 4 | Table 2, new calculation | `(7.27, 0.29, 1.04, -10.14)` | Same four calculated temperatures | 7.27029842e-9 | +9.00 | Landau–Zener total. |
| Ca 3 | Table 1, ref. 17 | `(0.0317, 2.12, 12.06, -0.40)` | 1,000–30,000 K printed fit range | 2.87964694e-10 | +6.60 | Fit to earlier calculation. |
| Ca 4 | Table 1, ref. 17 | `(2.68, 0.69, -0.68, -4.47)` | 1,000–30,000 K printed fit range | 2.65913841e-9 | +9.90 | Fit to earlier calculation. |
| Fe 4 | Table 2, new calculation | `(14.60, 0.0357, -0.92, -0.37)` | Same four calculated temperatures | 5.32205647e-9 | +10.00 | KF96 footnote `c`: radiative charge transfer may be faster. |

KF96 Table 2 gives no individual temperature-range column. Section 3.2 says
that its new Landau–Zener rates were calculated at `log10(T/K) = 3, 3.5, 4,
4.5`; this sampled span must not be relabeled as a published fit interval.
Cloudy `ctrecombdata.dat` independently uses 1,000–30,000 K for Cl 3–4 and
Fe 4. A numerical table on 1,000–30,000 K stays within the calculated span,
but endpoint holding beyond that span is only a computational continuation.
KF96 itself cautions against distant extrapolation of fits.

For N 4, both local MOCASSIN versions reproduce KF96's coefficients. Cloudy
c25.00 instead records `(2.95, 0.55, -0.39, -1.07)` over 1,000–1,000,000 K,
giving `2.5553687e-9 cm^3 s^-1` at 10,000 K versus the KF96 value
`1.8695385e-9`. KF96 cites Stancil et al. (1995) as the original source.
No published primary source for Cloudy's revision was established. Superseded
on 2026-10-11 (author's decision): Ferland et al. (2013, RMxAA 49, 137,
Sec. 2.1) state that Cloudy's charge-exchange rates are "an updated version of
the Kingdon & Ferland (1996) database", so N4+ (and Si2+, the only other row
where Cloudy differs) now take the Cloudy coefficients and intervals.

KF96 and KF99 also support Cl 3–4 and Ca 3–4, which were
omitted from the earlier nine-candidate audit. Cl 3–4 and Fe 4 are in KF96
Table 2, whereas Ca 3–4 are in Table 1; their heat values are in KF99
Table 1. Local MOCASSIN v2 and v3 contain
Ca coefficients but assign them to `chex(18,3:4)`, overwriting Ar instead of
placing them at atomic number 20; those rows cannot verify MOCASSIN's active
Ca rate. Both versions' `update_mod.f90` evaluate only their explicitly
populated `chex` rows and contain no rows above charge four in the inspected
table. Cloudy's `HCTRecom` reads the first four rows of its data file and
uses `HCTAlex*(ion+1)` for higher charges; the default `HCTAlex` is
`1.92e-9 cm^3 s^-1`. `ChargTranSumHeat` uses the file's heat values for the
first four charges and a charge-indexed `2.86*ion` eV estimate above them.
MAPPINGS IV `charex.f` evaluates tabulated rate coefficients, while `cheat.f`
multiplies those coefficients by separate stored energy defects. Neither
routine determines the missing high-charge electron-transfer probability.

KF99 Section 2 states that it uses the Langevin rates of Ferland et al.
(1997, *ApJ* 481, L115–L117; DOI `10.1086/310655`) for `q>4`, plus a mean
energy-defect estimate based on an avoided crossing at 9.5 bohr. The final
Ferland et al. journal pages were read from the author-affiliated Harvard
DASH repository. Its Section 2 assumes electron capture on *every* collision
that passes the centrifugal barrier. Equation (2) gives the Langevin rate
`3.85e-8 z sqrt(alpha*m/mu) cm^3 s^-1`, where `z` is excess nuclear charge,
`alpha` the neutral's polarizability, and `mu` the ion–atom reduced mass.
Equation (3) specializes this to Fe ions and atomic H:
`1.92e-9 z cm^3 s^-1`. For an ion `X^(q+)`, `z=q` in the article's charge
convention. The authors explicitly call equation (3) an *upper limit* and
expect it to become more accurate at higher charge. They give no temperature
validity interval; the model rate is independent of temperature. Section 3
tests Fe ions above charge four with the equation (3) rate and with zero
charge transfer to show sensitivity, not to establish a measured rate or
individual product channels. The selected implementation uses this formula as an explicitly
labeled unit-probability upper-limit approximation. It must not present `1.92e-9 q` as a
verified physical rate for all species or label 1,000–100,000 K as a
published validity range. The electron-transfer probability, rate by product
level, and radiative energy disposition remain unknown for N 5, S 5–6,
Ar 5–6, Fe 5–6, Cl 5–6, and Ca 5. KF99 mean heating and the capture bound
alone do not settle the physical rate for these boundaries.

The reverse `X^((q-1)+) + H+` process cannot be inferred by applying a
single Boltzmann factor to an excited-product total. KF96 Section 3.2 notes
that excited recombined levels generally have small populations; its Table 3
lists only selected ground-state reverse reactions. The mean heat in KF99 is
not a reverse threshold. Selective implementation should therefore use the
forward recombination rate and positive `CXHEAT` value, with `CXI=0`, pending
state-specific products, level populations, and detailed-balance data.

## Verification record

KF96 is the final *ApJS* 106, 205–211 PDF at
`../heii_h_charge_exchange/references/Kingdon_Ferland_1996_ApJS_106_205.pdf`;
Section 2, Section 3.2, Tables 1–3, and the Table 2 footnotes were read with
`pdftotext -layout`. KF99 is the final *ApJ* 516, L107–L109 publisher PDF at
`references/Kingdon_Ferland_1999_ApJ_516_L107.pdf`; Section 2 and Table 1
were read. The ADS API was attempted first but DNS resolution failed in the
sandbox. The final Ferland et al. article was subsequently obtained through
Harvard DASH's public item, bundle, and bitstream API at
`https://dash.harvard.edu/server/api/core/bitstreams/b6078519-b47b-4d7b-9980-576915b6b6b4/content`.
The PDF is `/nfs/mocafe/kiseon/MoCHII/references/Ferland_Korista_Verner_Dalgarno_1997_ApJ_481_L115.pdf`
(SHA-256 `612bf229d498dd20e4019dd8440282373af410237c2c7a0d4a6554ae68cd8af0`),
the intended permanent reference location; its source copy was first read
from `references/Ferland_Korista_Verner_Dalgarno_1997_ApJ_481_L115.pdf`;
it contains a Harvard cover page and the printed *ApJ* pages L115–L117.
The publisher URL `https://iopscience.iop.org/article/10.1086/310655/pdf`
served a CAPTCHA, while the ADS scan URL
`https://articles.adsabs.harvard.edu/pdf/1997ApJ...481L.115F` returned
HTTP 403. The Harvard copy was inspected directly with `pdftotext -layout`,
including its Section 2 derivation, Figure 1 caption, and Section 3
application. No preprint was substituted for published-method verification.
Code paths
inspected: Cloudy c25.00 `source/atmdat_char_tran.cpp`, `data/ctrecombdata.dat`,
`source/atmdat.cpp`; MOCASSIN v2/v3 `source/update_mod.f90`; MAPPINGS IV
`src/mastercode/charex.f` and `src/mastercode/cheat.f`. No production atomic
files or external code were modified, and no simulation suite was run because
this was a source and code-path investigation.
