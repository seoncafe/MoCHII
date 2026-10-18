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

- **Radiation field**: continuous stellar photon energies (13.6–100 eV) are
  sampled from a Planck or tabulated-spectrum CDF, with cross sections
  evaluated at each sampled energy.  Ionizing bins are diagnostic spectral
  tallies only; their count and edges do not enter the transport.  The band is
  extended into the FUV as a separate segment with its own packets
  (`add_fuv`), from 4.78 eV (the He I 2^3S edge) whenever something in the
  run absorbs below the H I edge: the He I 2^3S population (with helium),
  a metal whose neutral ionization potential lies below 13.598 eV, or dust.
  An input that sets `add_fuv` or `efuv_min` keeps its own values. The band
  has an analytic zero-variance estimator
  for the direct field, explicit diffuse recombination packets (case A,
  `recombination_continuum='transported'`, the default: the H/He ground
  continua, the He I and He III excited cascades, the metal ground
  free-bound continua and the free-free emission inside the band) or
  on-the-spot (case B), dust absorption and Henyey–Greenstein scattering
  in the band.  The three ground recombination continua are sampled from
  their detailed-balance Milne free-bound spectra
  (`diffuse_energy_model='milne'`), using the same cross sections transport
  absorbs with, so a diffuse packet is absorbed with the cross section it
  was emitted with.  The optional He I excited cascade splits into its
  `2^3S`, `2^1S` and `2^1P` branches with temperature-dependent fractions
  from one internally consistent set of effective recombination coefficients
  (Benjamin, Skillman & Smits 1999), and samples the `2^1S` two-photon
  continuum from the Drake, Victor & Dalgarno (1969) distribution that the
  nebular continuum already used.  The trapped He I 584 A resonance photon is
  not assumed to ionize hydrogen: it is split locally between hydrogen
  absorption and conversion to `2^1S` (the 2.06 micron decay of `2^1P`),
  which then decays as
  the two-photon continuum, with the escape-probability closure of Wood, Mathis
  & Ercolano (2004), and the nebular continuum uses the same split, so the
  transported packets and the written continuum describe one nebula.
  Quasi-random launch (`launch_sequence='sobol'`, the default;
  Owen-scrambled Sobol) for the stellar and diffuse packets — the
  launch set is independent of the MPI task count and the direct-field
  error scales as 1/N (effective photon gains of 29–222x measured on the
  analytic gate); see `docs/MoCHII_QMC.pdf`.
- **Sources**: any number of internal point sources, each with its own
  luminosity and spectrum (Planck temperatures or one multi-column
  spectrum file), plus — independently or in combination — an isotropic
  external radiation field of given mean intensity entering through the
  box faces or a bounding sphere (cosine-weighted entry), so a cloud
  embedded in an interstellar radiation field is modeled directly.
- **Gas physics**: H/He photoionization equilibrium and thermal balance;
  metals (C, N, O, Ne, S, Ar, Mg, Fe, Si, Cl, Ca) as trace species through
  ionization stages I-VII (77 stages, 66 transitions). Metal photoabsorption
  is Verner's `phfit2` (the Verner et al. 1996 outer-shell fit below the
  first inner-shell edge, the Verner & Yakovlev 1995 subshell sum above it);
  Ne III uses the Nahar (2019) R-matrix ground-term cross section below
  178 eV. A band above the lowest terminal ionization potential of the
  elements present (Cl VII, 114.2 eV) is refused. Badnell radiative +
  dielectronic recombination, Voronov (1997) collisional ionization
  (`metal_collisional_ionization`, on by default), and charge exchange with
  hydrogen in both balances (45 reactions with all eleven elements;
  `langevin_charge_exchange` adds the q >= 5 Langevin upper limits, off by
  default) with its kinetic heat (`charge_exchange_heating`, on by default).
  Metal excitation, cooling and lines come from one atom per stage, selected
  by `metal_excitation_model`: `'chianti'` (the default) takes full CHIANTI
  11.0.2 atoms for 67 stages and Stout atoms for Cl VII, Ar II and Ca IV;
  `'mocassin_v2'` takes the MOCASSIN v2 atoms, for parity runs. Every stage
  with an atom cools and emits; atomic-hydrogen impact on the C I, C II and
  O I fine structure (`neutral_hydrogen_excitation`) is on by default. The
  He I 2^3S population is coupled into the ionization and thermal balance
  (`hei_metastable_coupling`, on by default). See
  [atomic data](data/atomic/README.md) for coverage and the
  [current results](md/ATOMIC_PHYSICS_RESTORATION_RESULTS_2026-10-11.md).
  Removed and renamed inputs (`metal_inner_shell_absorption`,
  `additional_line_cooling`, `high_charge_exchange`, `fe_levels_full`,
  `pe_scale`, `metal_excitation_model='extended'`) are refused with their
  replacement named.
  Optional Shull & van Steenberg (1985) secondary ionization by fast
  photoelectrons, scaled in energy after Dalgarno, Yan & Liu (1999)
  (`use_sec_ion`, for hard fields and partially neutral gas),
  and net-charge-weighted metal free-free cooling (`metal_freefree`, on by
  default).  The thermal-balance line cooling is suppressed to the local
  electron density from the same n-level populations that emit the diagnostic
  lines (`cooling_model='local_ne'`, the default, matching the modern
  reference codes).
- **Diagnostics**: collisional line luminosities from n-level statistical
  equilibrium, Storey & Hummer (1995) H I / He II and Porter et al. He I
  recombination lines, the nebular
  continuum (free–bound, free–free, two-photon), the He I 2^3S metastable
  population with the 10833 A absorption opacity and column-density maps
  (`hei_metastable`), and leaf-by-leaf line emissivity output for map making.
- **Dust and PAHs**: grain absorption/scattering competing with the gas for
  ionizing photons, ionization-dependent dust survival with a PAH split,
  equilibrium dust temperatures, and stochastic dust/PAH emission spectra
  via the SEDust library (astrodust, DL07, MRN, Zubko grain models; MRN is
  the classic graphite + silicate power law and carries no PAHs at all, so a
  run on it has no aromatic features).  One grain model supplies both the
  transport cross sections and the emission, so the energy a cell absorbs is
  set by the same grains that radiate it back out: the transport reads the
  model's precomputed size-integrated extinction curve (`/kext` of
  `SEDust/data/<model>/sedust_<model>.h5`, the same product the cross
  sections come from) on the model's own wavelength grid, and a gate
  (`tests/grain_kext`) checks that curve against the size integral over the
  model as built, so a table belonging to other grains cannot pass
  unnoticed.  The astrodust, DL07 and MRN grids reach the whole ionizing band
  on their own: the T-matrix optics table spans 1.0e-4 to 3.981e4 um (1762
  wavelengths), so below the Lyman limit the grains are the b/a = 1.4 oblate
  spheroid solved on the DH21 dielectric function rather than a sphere
  approximation, and no extension is built at any photon energy the code
  transports.  MRN shares the DL07 axis, both being Mie on the D03
  dielectric functions.  Zubko is on its own DustEM grid, 0.001-10000 um.
- **Dust emission transport**: the reradiated infrared is itself launched and
  transported through the same grain opacity, and the reabsorbed part returns
  to the grain heating until it converges (`dust_emis_transport`, on by
  default).  The escaping infrared spectrum is written alongside the emitted
  one.  H II regions reabsorb well under a percent of their own infrared, so
  this converges in a couple of passes.
- **Imaging**: peel-off images of the direct and dust-scattered EUV/FUV
  field toward arbitrary observers (optionally the unattenuated direct image
  too, so its ratio to the direct image is the line-of-sight optical depth),
  bin-resolved image cubes with a self-describing band grid, and dust-band
  emissivities for infrared maps.
- **PDR zone**: FUV photoionization of low-threshold metals beyond the
  ionization front, metal electrons in the charge balance and metal
  photoheating (on with the metals), and, on by default with dust, grain
  photoelectric heating and accretion cooling (Weingartner & Draine 2001,
  for the grains of `dust_model`, with PAHs weighted by their survival
  x_HI + `f_ion_pah` x_HII) and neutral gas-grain energy exchange
  (Hollenbach & McKee 1979).
- **Grids**: adaptive octree read from a generic AMR file (RAMSES /
  Illustris-TNG converters available), with optional solution-driven
  re-refinement at the ionization front; or a uniform Cartesian grid
  (raster storage, integer-arithmetic traversal, no tree in memory) built
  from the namelist, from a leaf list, or from a 3D density cube
  (`density_file`, FITS/HDF5) read directly onto the grid.
- **Plane-parallel slab** (`xy_periodic`): a horizontally infinite slab —
  the medium varies only in z (a 1D column) or is an arbitrary 3D tile
  repeated periodically in x and y — illuminated through the top and/or
  bottom face by a collimated beam at a chosen incidence angle or an
  isotropic field, each face with its own strength
  (`source_geometry='slab'`).  The incident flux is the sole normalization
  (results scale with it and are independent of the tile area); the
  emergent observable is the boundary intensity I(mu) with the
  reflected/transmitted/absorbed energy budget.  Internal point sources
  also work under the same periodic boundary.
- **Parallelism and I/O**: MPI with MPI-3 shared memory (one grid copy per
  node); the photons of a transport pass are split statically over the
  ranks, or served in batches on demand by rank 0 with `use_master_slave`
  (worth ~20% of the ionizing-band pass at 70 ranks, and it cuts the
  run-to-run spread of that pass from 28% to 1.5%); HDF5 or FITS output
  through a format-agnostic interface; Python readers, a 2D map maker, and
  slice/cutout viewers for the grid and the outputs under `tools/python/`,
  plus an AMR octree grid builder (`tools/python/AMR_grid/`).

MoCHII shares its transport engine with the author's dust radiative-transfer
code [MoCafe](https://github.com/seoncafe/MoCafe) and is validated against
analytic Strömgren solutions, the Lexington H II-region benchmarks (against
the published multi-code results and against current Cloudy and MOCASSIN
runs), and PyNeb emissivities.  On the Lexington HII40 and HII20 cases the
mean electron temperature `<T[NpNe]>` is 8176 K and 6942 K, within 0.41% and
0.46% of Cloudy c25.00 (8210 K and 6910 K), and L(Hbeta) is 2.025e37 and
4.85e36 erg/s, inside the published ranges (converged runs of 2026-10-11,
`results/continuous_energy/m13_restoration/`).

## Build and run

The SEDust dust-emission library is self-contained under `SEDust/`. On a
fresh checkout, build it once before the first `make`:

```
cd SEDust/sed && ./build_lib.sh   # -> SEDust/sed/lib/libsedust.a
```

Then, from the repository root:

```
make                              # -> MoCHII.x  (MPI Fortran + HDF5)
mpirun -np 8 ./MoCHII.x input.in
```

`par%sed_data_dir` is resolved automatically relative to the executable at run
time, so dust emission works from any working directory without setting a
path. Everything the dust model is made of resolves inside that one directory.
`par%cext_dust`, the reference extinction every leaf's gray opacity is built
on, is likewise resolved from the grains the run carries unless a run states
one outright, so naming `par%dust_model` is the whole dust physics.

The optics of the named `par%dust_model` come from
`SEDust/data/<model>/sedust_<model>.h5`, which carries the wavelength axis, the
cross-section tables and the size-integrated extinction curve together. The
default astrodust axis runs `1.0e-4`--`3.981e4` um (1762 wavelengths), spanning
the ionizing band MoCHII transports.

See `docs/MoCHII_UserGuide.pdf` for the input-parameter reference, output
formats, and worked examples; `docs/MoCHII_physics.pdf` for the atomic
data, algorithms, and validation results; `docs/MoCHII_fitting.pdf` for
the CHIANTI fitting procedure; `docs/MoCHII_cooling_analysis.pdf` for the
line-cooling verification against other photoionization codes; and
`docs/MoCHII_code_comparison.pdf` for the algorithmic and physical comparison
with M^3, MAPPINGS V, ionize, CMacIonize, Torus, and MOCASSIN.

## Author

Kwang-Il Seon (KASI / UST)

---

Last updated: 2026-10-18 16:25 KST
