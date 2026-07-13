# MoCHII v1.00 <img src="docs/mochii_icon.png" width="64" alt="MoCHII icon" align="top">

MoCHII (**Mo**nte **C**arlo for **H II** regions) is a Monte-Carlo
photoionization and radiative-transfer code (Fortran 90 + MPI) for **dusty
photoionized nebulae** on adaptive octree (AMR) and uniform Cartesian grids.
One self-consistent radiation field drives both the gas and the dust: the
code transports ionizing/FUV photon packets from the source (and from the
recombining gas), solves the H/He ionization and thermal balance of every
cell together with a trace-metal ionization cascade, and produces the
observables of the nebula — emission-line luminosities and emissivity maps,
recombination lines, the nebular continuum, dust temperatures, and the
infrared dust/PAH spectrum.

## Features

- **Radiation field**: log-spaced ionizing bins (13.6–100 eV, extendable
  into the FUV down to ~6 eV as a separate band segment), Planck or
  tabulated source spectra, an analytic zero-variance estimator for the
  direct field, explicit diffuse recombination packets (case A) or
  on-the-spot (case B), dust absorption and Henyey–Greenstein scattering
  in the band.
- **Gas physics**: H/He photoionization equilibrium and thermal balance;
  metals (C, N, O, Ne, S, Ar, Mg, Fe, Si, Cl, Ca) as trace species with Verner et al.
  (1996) cross sections, Badnell radiative + dielectronic recombination,
  Voronov (1997) collisional ionization, and charge exchange; all atomic
  rates are fitted offline from CHIANTI and evaluated as closed forms at
  run time.  Adding an ion is a data operation, not a code operation.
- **Diagnostics**: collisional line luminosities from n-level statistical
  equilibrium, Storey & Hummer (1995) H I recombination lines, the nebular
  continuum (free–bound, free–free, two-photon), and leaf-by-leaf line
  emissivity output for map making.
- **Dust and PAHs**: grain absorption/scattering competing with the gas for
  ionizing photons, ionization-dependent dust survival with a PAH split,
  equilibrium dust temperatures, and stochastic dust/PAH emission spectra
  via the SEDust library (astrodust, DL07, Zubko grain models).
- **Imaging**: peel-off images of the direct and dust-scattered EUV/FUV
  field toward arbitrary observers, and dust-band emissivities for infrared
  maps.
- **PDR zone (optional)**: FUV photoionization of low-threshold metals
  beyond the ionization front, metal electrons in the charge balance,
  metal photoheating, and grain photoelectric heating (Bakes & Tielens
  1994).
- **Grids**: adaptive octree read from a generic AMR file (RAMSES /
  Illustris-TNG converters available), with optional solution-driven
  re-refinement at the ionization front; or a uniform Cartesian grid
  (raster storage, integer-arithmetic traversal, no tree in memory).
- **Parallelism and I/O**: MPI with MPI-3 shared memory (one grid copy per
  node); HDF5 or FITS output through a format-agnostic interface; Python
  readers and a 2D map maker under `tools/python/`.

MoCHII shares its transport engine with the author's dust radiative-transfer
code [MoCafe](https://github.com/seoncafe/MoCafe) and is validated against
analytic Strömgren solutions, the Lexington/MOCASSIN H II-region benchmarks,
and PyNeb emissivities.

## Build and run

```
make                              # -> MoCHII.x  (MPI Fortran + HDF5)
mpirun -np 8 ./MoCHII.x input.in
```

See `docs/MoCHII_UserGuide.pdf` for the input-parameter reference, output
formats, and worked examples; `docs/MoCHII_physics.pdf` for the atomic
data, algorithms, and validation results; and `docs/MoCHII_fitting.pdf` for
the CHIANTI fitting pipeline.

## Author

Kwang-Il Seon (KASI / UST)

---

Last updated: 2026-07-13 12:32 KST
