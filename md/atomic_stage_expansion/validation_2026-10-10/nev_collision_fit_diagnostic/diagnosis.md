# Ne V 6–15 collision-strength fit

**Numerical-representation result only.** The source cubic spline produces
negative strength between nonnegative nodes; preserving its clipped zero is
not a physically acceptable interpolation correction. See
`positivity_interpolation_assessment.md` for the subsequent source audit.

The 0.2% native fit criterion is achievable with the existing 16-coefficient,
type-2 storage format. The source effective collision strength remains the
CHIANTI cubic-spline evaluation. Only the invertible coordinate used to store
its approximation changes.

The source is `~/RT_Codes/CHIANTI/dbase/ne/ne_5/ne_5.scups`, SHA-256
`387bdb0d55dec9cce35edf07a9d99faa31b3e4814f24495ed5b4398b14415584`.
Pair 6–15 has source type 2, `de = 2.144 Ry`, and source scaling constant
`C = 0.9492`. All 15 tabulated scaled strengths are nonnegative. The first
two are 0 and 0.000147. The cubic spline nevertheless reaches
`-8.3216173394e-5` near 1729.12 K. The current `upsilon()` implementation
clips this negative interpolation to zero. Its exact positive crossing is
3622.576354426 K. This is an interpolation artifact; it does not establish
negative collision strengths in the underlying calculation.

The existing degree-15 fit has a maximum relative error of
0.00283478527513 on both the 201-point production grid and a 1601-point
logarithmic grid from 1000 to 100000 K. The worst 1601-point error occurs at
3736.801 K. Its fitted strength is zero wherever the source evaluation is
zero on that grid.

For type 2, the storage coordinate is `s = (kT/de)/(kT/de + C)`, while the
descaled strength is the fitted polynomial itself. Thus a new positive `C`
changes only the coordinate of the approximation, provided the target values
are still evaluated using the source `C = 0.9492`. The recommended native
coordinate uses `C = 0.04`. A linear minimax fit of a degree-15 Chebyshev
polynomial uses relative-error constraints at the 1601 positive target values,
nonpositive polynomial constraints in the zero region, and exact source
zero-crossing value and slope constraints. The printed native fields are
rounded to their actual file precision; the constant coefficient is then
adjusted at that precision to minimize the root residual.

The rounded row is in `result.json` at `selected_native_storage.native_row`.
Evaluated with those rounded fields, its maximum relative error is
0.0005519743 on 201 points, 0.0005519897 on 1601 points, and 0.0005521539
on 160001 points. The fitted result is zero at every sampled source-zero
temperature on all three grids. The 160001-point worst temperature is
97000.729 K. The row uses all 16 available coefficients, `logf = 0`, and
requires no reader or data-layout change. It must be generated from the source
by fitting code rather than hand-editing a produced table.

Relative error is ill-conditioned arbitrarily close to the source's artificial
zero crossing. Even after rounding-aware adjustment, the polynomial residual
at the exact crossing is `1.6687e-12`; the relative error is about 0.20% at
crossing +0.01 K and 0.0199% at crossing +0.1 K. The dense-grid result is
therefore a sampled-domain result, not a uniform mathematical bound at the
zero. The source cubic spline itself may deserve a separate interpolation
review. Replacing it with a shape-preserving interpolator would change the
physical data evaluation and must be assessed separately from this
representation repair.

`diagnose_nev_6_15.py` reproduces the source audit, the coordinate scan, the
root-constrained fit, and the rounded native row. It changes no shared source
or production output.
