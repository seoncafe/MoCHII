# data/atomic

Fitted coefficient files produced by `tools/fitting`. Every file carries a
provenance header: source and version of the atomic data, fit form, fit
range, and maximum fit error.

## Collisional line-cooling fits

The cooling used in the thermal balance is stored as closed-form fits, one
file for each ion (`cooling_<element>_<stage>.txt`, e.g. `cooling_o_3.txt`
for O III). `cooling_fit_parameters.txt` collects all of them into a single
listing, generated from those files by
`tools/fitting/make_cooling_table.py`.

The fitted quantity is the rate coefficient

```
Lambda_ion(T) = T^(-1/2) * sum_{i=1}^{N} A_i * exp(-T_i / T)
```

| symbol | meaning | unit |
|---|---|---|
| `Lambda_ion` | line-cooling rate coefficient of the ion | erg cm^3 s^-1 |
| `A_i` | amplitude of term *i* | erg cm^3 s^-1 K^(1/2) |
| `T_i` | excitation temperature of term *i* | K |
| `T` | electron temperature | K |
| `N` | number of terms (1 to 8, one line for each) | — |

The volumetric cooling rate of the ion is `n_e * n_ion * Lambda_ion(T)`
in erg cm^-3 s^-1; MoCHII sums it over the ions of the active registry.

File layout of `cooling_<ion>.txt`: comment header (lines starting with
`#`), then `N`, then `N` lines of `A_i  T_i`.

The fits are built from CHIANTI v11.0.2 level energies and effective
collision strengths (Burgess–Tully descaling) in the optically thin
n_e -> 0 limit of statistical equilibrium, over 10^3–10^5 K: all ions sit
in their lowest level, collisional excitation happens out of it only, and
each excitation radiates the full excitation energy through the cascade, so
`Lambda_ion` is the excitation power out of the ground level.  Transition
energies are the observed level differences, the same energies the n-level
solve uses in both the Boltzmann factor and the emitted photon.  The fits
therefore carry no density suppression above the critical density, and
overestimate the cooling of the low-n_crit lines ([C II] 158 um, the
[O III] 52/88 um pair, the [O II] and [S II] optical doublets) in dense
gas; the density-dependent line emissivities used for the diagnostics come
from the separate n-level solve (`nlevel_<ion>.txt`).

Maximum fit errors are recorded in each header: at or below ~0.7% across
the default C/N/O/Ne/S set except C III (5.5% full range, 0.50% where the
ion cools) and S II (2.2%), with larger errors (up to ~28% for Si III, 7–9%
for the Cl ions) for the heavily depleted Si, Cl, and Ca ions, whose
contribution to the thermal balance is ~10^-4 at their abundances.

## Other files

| file | content |
|---|---|
| `nlevel_<ion>.txt` | n-level model: levels, A-values, Chebyshev fits of the effective collision strengths |
| `element_<el>.txt` | photoionization, recombination, collisional-ionization, and charge-exchange data of one element |
| `badnell_rr.dat`, `badnell_dr.dat` | Badnell radiative and dielectronic recombination fit parameters |
| `sh95_*.txt`, `hei_porter_caseB.txt` | recombination-line emissivity tables |
| `hei_metastable.txt` | He I 2^3S metastable rate coefficients and the 10830 A line constants |
| `gauntff_vh14.dat`, `gamma*.dat` | free-free Gaunt factors and free-bound emission tables |
