# examples

Run templates, following the MoCafe examples/ pattern (one directory per case
with an input namelist + `run.sh`). Each `run.sh` launches the case with
`mpirun -np 8 ../../MoCHII.x <name>.in`; outputs are written next to the input
with the `out_file` basename.

Most cases build a uniform-density sphere or shell directly from the namelist
(`grid_type='car'` with `nx/ny/nz`, `xmax/ymax/zmax`, `nH_const`, and `rmax` or
`rmin`+`rmax`), so no grid file is needed. The cases that read a grid or density
file keep that file in their own directory.

| Example | What it demonstrates |
|---|---|
| `stromgren_sphere` | Classic Stromgren sphere: H + He photoionization equilibrium around a central hot star, at fixed electron temperature. |
| `stellar_spectrum` | Photoionization driven by a young stellar-population SED (`par%ion_spectrum`) instead of a blackbody; `make_spectrum.py` converts a 10 Myr FSPS population, which is soft (its hot stars have died) and barely ionizes helium. |
| `physical_spectrum` | Stellar spectrum in PHYSICAL units (`par%spectrum_type='per_ang'`, L_lambda vs lambda): the ionizing luminosity is derived from the file, not set with `par%luminosity`, and printed in the log; `make_physical_spectrum.py` writes a 4e4 K blackbody normalized to the `stromgren_sphere` luminosity, so the ionization structure matches that example. |
| `hii_region` | H II region with metal cooling, a self-consistent electron temperature, and emission-line + nebular-continuum output. |
| `dusty_nebula` | Dusty photoionized sphere: grains absorb and forward-scatter in the ionizing band, and grain heating sets an equilibrium dust temperature with an infrared spectrum. |
| `pdr` | Photodissociation region: the FUV band drives grain photoelectric heating while metal electrons and metal-line cooling set the temperature beyond the ionization front. |
| `cosmological_cutout` | Illustris-TNG density cutout ionized by a source at the density peak, with the octree re-refined on the ionization front. |
| `adaptive_refinement` | Solution-driven octree re-refinement of the ionization front, resolving the front on a coarse base grid. |
| `peeloff_imaging` | Peel-off imaging: direct and dust-scattered light in the ionizing and FUV bands projected onto an image plane toward an observer. |
| `density_cube` | A 3D FITS density cube read directly onto a Cartesian grid (`make_cube.py` builds a clumpy cube). |
| `external_field` | A cloud in an isotropic external ionizing field (`par%ext_intensity`, no internal star): the field ionizes a skin around a neutral core, and the peel-off image shows the core in silhouette against the background (`direc0/direct` = the line-of-sight `e^tau` map). |
| `isrf_cloud` | A dense neutral cloud bathed in the standard interstellar radiation field (an ISRF preset, `par%ext_spectrum='draine'`, FUV-only, requires `add_fuv`): the FUV field carries no ionizing photons but photoionizes low-threshold carbon and magnesium in a dust-attenuated skin, leaving a C II skin over a C I-dominated core (`par%metal_ne` supplies the electrons that let carbon recombine). |
| `multiple_sources` | Two stars of different temperature and luminosity in one uniform medium (`par%nsource=2`, `src_lum`/`src_tstar`): overlapping H II regions of different size, with an He III core and a larger He II zone around the hotter star; `make_two_spectra.py` writes the tabulated multi-column spectrum alternative (`par%src_spectrum_file`). |
