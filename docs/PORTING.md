# MoCHII port map

What each source module comes from and what changes on the way in.
Rule: copy, never edit the origin trees. Keep the origin file's header
comment and add one line: `! MoCHII: copied from <tree>/<path> (YYYY-MM-DD)`.

## From MoCafe v2.00 (`~/MoCafe/new/MoCafe_v2.00/src/`)

Take as-is (infrastructure, no physics changes expected):

| Module | Role |
|---|---|
| `octree_mod.f90` | AMR tree build, leaf find, cell exit, neighbor walk |
| `read_generic_amr.f90` | snapshot reader (generic AMR format) |
| `memory_mod_mpi.f90` | MPI shared-memory arrays |
| `random_mt.f90` | Mersenne-Twister RNG |
| `mathlib.f90`, `utility.f90` | numerics + helpers |
| `hdf5io_mod.f90`, `fitsio_mod.f90`, `iofile_mod.f90` | I/O |
| `peelingoff_mod.f90` | (not copied — ADAPTED into `ion_peel_mod`, see below) |
| `scattering_car.f90`, `sources_mod.f90`, `mrw_mod.f90` | (not taken — the band transport implements its own HG scattering inside `raytrace_amr`; MRW is a dust-band device) |

Take and extend:

| Module | Extension |
|---|---|
| `define.f90` | add gas leaf state (`nH`, `x_HI`, `x_HeI`, `x_HeII`, `T_e`, `n_e`, registry ion-fraction block); ionizing-band parameters |
| `grid_mod_amr.f90` | persist `nH`/`xHI` (today deallocated after the dust-density step); allocate gas state |
| `jtally_mod.f90` | ionizing-band J tally (nbin_ion x nleaf) alongside the SED tally |
| `raytrace_amr.f90` | band-dependent gas+dust opacity in the ionizing bins (no grey `s_ext` rescale there) |
| `lucy_mod.f90` | (not copied — the gas iteration is `main.f90`'s loop over `ionizing_field_and_rates` (`ionizing_field_mod`): the tally is rebuilt from zero each iteration since the opacity changes) |
| `run_simulation_mod.f90` | (not copied — ADAPTED into `photon_schedule_mod`, see below) |
| `sed_mod.f90` | copied (band/SED grids) |
| `dustemis_mod.f90` (+ SEDust) | ADAPTED into `sedust_mod` (SEDust `dust_lib` built from source under `SEDust/sed`): leaf spectra shaped by the local J, normalized to Heat_dust V; live PAH channel weighting; dust-band emissivities |
| `observer_mod.f90` + `peelingoff_mod.f90` | ADAPTED into `ion_peel_mod`: MoCafe observer geometry + direct/scattered peel for the ionizing/FUV band, tau by a tally-free walk, hook into the transport |
| `setup.f90`, `main.f90` | trimmed drivers (new files, MoCafe conventions) |

## From MOCASSIN (`~/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/`)

| Item | Notes |
|---|---|
| `source/gaunt.f90` | free-free Gaunt factors (Hummer 1988). Free-standing, zero dependencies — verbatim copy. DONE: the van Hoof et al. (2014) table is available as `par%gaunt_vh14` (`gaunt_vh14_mod`, `data/gauntff_vh14.dat`); the two agree to <= 0.5% in T_e |
| `emission_mod.f90::fb_ff` block + `ph_mod.f90::initLogGamma*` | nebular continuum: free-bound (gamma tables) + free-free + two-photon. Ported into `nebcont_mod`; the tables ship in the ES06 compact paired-edge format and were numerically verified fb+ff combined (no swap needed) |
| `data/gammaHI.dat`, `gammaHeI.dat`, `gammaHeII.dat`, `HeI2phot.dat` | continuum data files -> `data/` |
| (reference only) `tests/` suite | 3D validation targets: Stromgren/HII20/HII40 Te and line fluxes |

## From EXHALE (`~/RT_Codes/Exoplanetary_Atmospheres/ATES/EXHALE/`)

| Item | Notes |
|---|---|
| `cooling_data/chianti_cooling.py` | dependency-light CHIANTI v11 reader + Burgess-Tully descaling -> seed of `tools/fitting/` |
| `cooling_data/*` notebooks + `cno_formula_coefficients.txt` | fit pattern Lambda(T) = T^(-1/2) sum A_i exp(-T_i/T); C/N/O already fitted (coronal limit, 1e3-1e5 K) |
| rate treatments (docs + code) | Verner+96 photoionization cross sections, Badnell RR+DR, Voronov collisional ionization, Kingdon & Ferland charge exchange (`docs/charge_exchange_table4.md`) |

## New modules (no origin)

| Module | Role |
|---|---|
| `gas_state_mod` | leaf gas state + species registry |
| `gas_rates_mod` | Gamma / heating integrals from the ionizing-band J tally |
| `ion_balance_mod` | H/He (+ element cascade) equilibrium solve |
| `thermal_mod` | heating - cooling balance for T_e (Tier-1 fitted Lambda) |
| `nlevel_mod` | generic n-level atom solve (Tier-2, output-time diagnostics) |
| `gas_opacity_mod` | ionizing-band absorption coefficients from the state |
| `photo_xsec` | VFKY96 + VY95 analytic photoionization cross sections |
| `ion_band_mod` | band grids, source sampling, FUV extension |
| `recomb_mod` | Hui & Gnedin H/He recombination + Voronov CI |
| `cooling_mod` | thermal-loop cooling (rec, ff, CI, Tier-1 lines) |
| `species_mod` | 11-element registry: cascade, cooling/heating, opacities, n_e closure |
| `diffuse_mod` | explicit diffuse field (3 ground channels + He I excited channel) |
| `amr_refine_mod` | solution-driven I-front re-refinement + window recycling |
| `dust_temp_mod` | equilibrium mixture T_dust + IR spectrum |
| `sh95_mod` | Storey & Hummer H I / He II case-B line tables |
| `lines_mod` | line luminosities + leaf emissivity output |
| `nebcont_mod` | nebular continuum (fb+ff tables, Milne, two-photon) |
| `gaunt_vh14_mod` | van Hoof et al. (2014) free-free Gaunt table |
| `ion_peel_mod` | peel-off imaging (adapted from MoCafe, see above) |
| `sedust_mod` | SEDust coupling (adapted from MoCafe, see above) |
| `ionizing_field_mod` | one ionizing-band transport pass + the rate integrals it feeds: stellar and diffuse launch, transport, reduction, `gas_rates_compute`. `ionizing_field_and_rates(label, with_peel)` serves the iterating passes, the final consistency pass and the peel-off imaging pass, so the field the written rates carry is built exactly like the field that drove the iteration. Lifted out of `main.f90`, where it existed twice (once for the iteration, once copied for imaging) |
| `photon_schedule_mod` | how one pass's global photon indices are split over the ranks, behind `photon_schedule_begin`/`_next`/`_end` so the loop body is written once. ADAPTED from MoCafe's `run_simulation_mod.f90`, which keeps the two schedules as two whole routines with the transport loop inlined in each (`run_master_slave`, `run_equal_number`); separating the schedule from the loop body is what lets the stellar loop, the diffuse LyC loop and the peel-off pass share one implementation. `par%use_master_slave = .false.` (MoCHII default; MoCafe defaults `.true.`) reproduces the previous static cyclic sequence index for index; `.true.` makes rank 0 serve batches of `par%num_send_at_once` on demand and transport nothing |
