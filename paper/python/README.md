# Figure scripts for the MoCHII code paper

Every figure in `mochii.tex` is produced by one script in this directory.
The scripts read the existing run outputs under `../../tests/` (no MoCHII
run is needed) and write the figure into `../figures/` as a PDF, so they
can be run from anywhere:

```
cd paper/python
python3 make_g1_stromgren.py
```

| script | figure file | paper figure |
|---|---|---|
| `make_g0_gamma.py` | `g0_gamma_check.pdf` | photoionization rate vs the analytic solution |
| `make_g1_stromgren.py` | `g1_stromgren_check.pdf` | Stromgren H and He ionization structure |
| `make_g2a_thermal.py` | `g2a_thermal_check.pdf` | electron temperature with thermal balance |
| `make_g4_refine.py` | `g4_refine_check.pdf` | ionization-front re-refinement |
| `make_g4_tng.py` | `g4_tng_maps.pdf` | Illustris-TNG cutout maps |
| `make_dust_stromgren.py` | `dust_stromgren.pdf` | H II region size versus dust treatment |
| `make_dust_maps.py` | `pahlive_maps.pdf`, `dustemis_maps.pdf` | PAH survival and dust-band maps |
| `make_pdr_profiles.py` | `pdr_profiles.pdf` | photodissociation-region radial profiles |
| `make_wood_benchmark.py` | `wood_hii40.pdf`, `wood_hii20.pdf` | Lexington HII40/HII20 radial structure, in the panel layout of Wood et al. (2004) |

Conventions used by all of them:

- **PDF output** (vector), matching the rest of the paper.
- **`text.usetex=True`**, so every label, title, and legend entry is a
  LaTeX string and carries no raw Unicode.
- **Axis, tick, and title fonts are scaled 1.5x** over the matplotlib
  defaults for readability at the printed size.
- **Monte Carlo point clouds are randomly sub-sampled for plotting only**
  (`THIN` at the top of each script; a fixed seed keeps a figure
  reproducible). Every number quoted in a script's console output or in
  the paper is computed from the full cell set, never from the sample.

The atomic-data tables of the appendix are generated separately by
`../tables/make_atomic_tables.py`.
