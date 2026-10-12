# Lexington benchmarks with Cloudy c25.00

Reference runs of the two Lexington H II region benchmarks and the PN150
planetary nebula with a current version of Cloudy, kept here so the
comparison in the paper can be reproduced without a Cloudy installation.

The published GF column of the benchmark tables is the 1990s Cloudy, so a
current run is what tells us how much of a difference between MoCHII and the
published values is code and how much is thirty years of atomic data.

## Inputs

| file | benchmark | source |
|---|---|---|
| `hii40_c25.in` | HII40 (40 kK blackbody) | physics block of `tsuite/auto/hii_paris.in`, unmodified |
| `hii20_c25.in` | HII20 (20 kK blackbody) | physics block of `tsuite/auto/hii_coolstar.in`, unmodified |
| `hii40_cool.in` | HII40 | `hii40_c25.in` physics block plus the zone-by-zone cooling (`save cooling`, `save cooling each`), heating and element saves read by `fixed_state_cooling.py` |
| `pn150_c25.in` | PN150 (150 kK blackbody) | physics block of `tsuite/auto/pn_paris.in`, unmodified (see the PN150 section) |
| `wood_lines.dat` | line list read by `save line list` | written for this comparison |
| `pn150_lines.dat` | PN150 line list read by `save line list` | written for this comparison |
| `orig_hii_paris.in`, `orig_hii_coolstar.in`, `orig_pn_paris.in` | untouched copies of the distributed decks | kept so the physics block can be diffed against ours |

Only output commands were appended to the distributed decks; the physics
(source, density, radius, abundances, stopping criterion) is untouched, so
the runs are the Cloudy test-suite benchmarks as shipped.

## Reproducing

```
cd tests/cloudy_c25
~/CLOUDY/c25.00/source/cloudy.exe -r hii40_c25
~/CLOUDY/c25.00/source/cloudy.exe -r hii20_c25
```

Each takes about 90 seconds.

## Outputs used by the paper

| file | content | used for |
|---|---|---|
| `*_abs.lin`, `*_rel.lin` | line list, absolute and relative to H$\beta$ | the C25 column of the benchmark tables |
| `*.ovr` | overview: T_e and the main ion fractions against depth | C25 temperature; the volume-averaged `Te(NeNp)` is read from the `Averaged Quantities` block of `*.out` |
| `*.ele_{h,he,c,n,o,ne}` | fraction of each element in every ionization stage against depth | the Cloudy curves overplotted on the radial-structure figures |
| `*.rad` | radius against depth | radius axis of those figures |
| `*.out` | full log | volume-averaged quantities, convergence, provenance |
| `*.dr` | zone thicknesses | not used; kept as the run record |

The blend composition assumed when comparing a Cloudy line with a benchmark
row (which components are summed into `blnd 2326`, `blnd 3727`, and so on)
is documented in `wood_lines.dat` and in
`MoCHII_paper/tables/make_wood_tables.py`; it was read from c25.00's own
`data/blends.ini`.

## Consumers

The paper sources live beside this repository, in `../MoCHII_paper/`
relative to the repository root; its scripts read these files by relative path.

- `MoCHII_paper/tables/make_wood_tables.py` — the C25 column of the benchmark tables
- `MoCHII_paper/python/make_wood_benchmark.py` — the Cloudy curves on the radial figures
- `fixed_state_cooling.py` (run with `run_fixed_state.sh`, which builds the n-level
  driver `fixed_state_nlevel.f90` against `src/*.o`) — MoCHII's cooling terms
  evaluated on the `hii40_cool.*` state and integrated channel by channel
  against Cloudy's, plus the temperature-scaling estimate of the paper's
  fixed-state cooling section. Every metal file the driver reads (fit, n-level
  atom, hydrogen-impact rates) comes from one excitation model,
  `data/atomic/excitation/<model>/`, resolved as a run resolves it; the code
  default unless `--excitation-model chianti|mocassin_v2` is given. An ion the
  model gives no fit takes the n-level table, as `species_mod::metal_cooling`
  does.

## PN150 (Lexington/Meudon "Paris" planetary nebula)

`pn150_c25.in` is the physics block of `tsuite/auto/pn_paris.in` (lines 1
and 3-23, kept as distributed in `orig_pn_paris.in`) with output commands
appended. Its line 2, `set save prefix`, is left out so the save files carry
the names given; its own `normalize to "Ca b" 4861` is an output command and
is not carried, so `pn150_c25_rel.lin` is relative to the H I 4861 total.
The parameters match Jin, Kewley & Sutherland (2022, ApJ 927, 37) Table 1
entry for entry: 150 kK blackbody of radius 1e10 cm, n_H = 3000 cm^-3,
R_in = 1e17 cm, sphere, He/C/N/O/Ne/Mg/Si/S as published, Cl/Ar/Fe off. The
deck has no stopping command, so the run ends at Cloudy's default 4000 K
floor, at x(H0) = 0.967. No CMB, no grains, one iteration.

`pn150_lines.dat` carries the Jin et al. Table 4 PN150 lines, one entry for
each component of a published blend, plus H I 4861, He II 4686 and their
case-B values. Its header lists the composition of every Cloudy blend used,
read from c25.00 `data/blends.ini`. Jin's row "N II] 1749+" is the N III]
1750 multiplet (Cloudy `blnd 1750`); N II has no line there.

```
cd tests/cloudy_c25
~/CLOUDY/c25.00/source/cloudy.exe -r pn150_c25
```

The run takes about 75 s (187 zones). Measured on 2026-10-11:

| quantity | Cloudy c25.00 |
|---|---|
| L(H I 4861) | 2.8966e35 erg/s (Jin median 2.68e35, ratio 1.081) |
| L(Ca B 4861) / L(H I 4861) | 0.99297 (to renormalize a Ca B-normalized list, x1.0071) |
| outer radius (4000 K) | 4.1706e17 cm = 0.13516 pc |
| H front, x(H0) = 0.1 / 0.5 / 0.9 | 0.1294 / 0.1327 / 0.1345 pc |
| He III front, x(He2+) = 0.9 / 0.5 / 0.1 | 0.0805 / 0.0899 / 0.0959 pc |
| Te(NeNp), volume | 1.23e4 K |
| T at the illuminated face | 1.880e4 K |

Ratios to the Jin et al. Table 4 median, relative to H I 4861 (the
`Ca B`-normalized values are 1.0071 times these): [S II] 6716+6731 1.717,
[S II] 4069+4076 1.827, [S III] 9532+9069 0.983, [S III] 18.7 um 0.961,
[S III] 33.5 um 0.894, He II 4686 0.816, [Ne III] 3869+3968 0.723,
[O III] 5007+4959 1.081. The complete comparison, with MoCHII, is
`results/continuous_energy/m13_restoration/compare_pn150.py`.
