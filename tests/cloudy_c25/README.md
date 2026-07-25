# Lexington benchmarks with Cloudy c25.00

Reference runs of the two Lexington H II region benchmarks with a current
version of Cloudy, kept here so the comparison in the paper can be
reproduced without a Cloudy installation.

The published GF column of the benchmark tables is the 1990s Cloudy, so a
current run is what tells us how much of a difference between MoCHII and the
published values is code and how much is thirty years of atomic data.

## Inputs

| file | benchmark | source |
|---|---|---|
| `hii40_c25.in` | HII40 (40 kK blackbody) | physics block of `tsuite/auto/hii_paris.in`, unmodified |
| `hii20_c25.in` | HII20 (20 kK blackbody) | physics block of `tsuite/auto/hii_coolstar.in`, unmodified |
| `wood_lines.dat` | line list read by `save line list` | written for this comparison |
| `orig_hii_paris.in`, `orig_hii_coolstar.in` | untouched copies of the distributed decks | kept so the physics block can be diffed against ours |

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
`paper/tables/make_wood_tables.py`; it was read from c25.00's own
`data/blends.ini`.

## Consumers in this repository

- `paper/tables/make_wood_tables.py` — the C25 column of the benchmark tables
- `paper/python/make_wood_benchmark.py` — the Cloudy curves on the radial figures
