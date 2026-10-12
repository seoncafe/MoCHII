# Conditional full-bound O V numerical validation

This private sequence tests the current complete O V bound atom: 166 levels,
2,995 radiative rows, 13,695 collision pairs, and 205,425 original scaled
collision nodes. Its results are conditional numerical diagnostics. Original
photon-frequency and radiative-strength conventions remain under physical
review, so every result records `physical_model_accepted=false` and
`source_accuracy_verified=false`. No product is installed automatically.

The distributed Stout O V collision table gives 13 temperature nodes from
5,000 to 50,000,000 K. Inversion of the rounded CHIANTI BT abscissae gives
4,980.51 to 5,021.38 K for the first finite node of all 13,695 pairs. The
low-40-level source pairs agree with Stout at 5,000 K within 0.3174% after
rounding and interpolation. The final 2014 Be-like paper's Methods text
instead states 16,000 to 160,000,000 K for its calculated effective collision
strengths. The relation between that sentence and the distributed O V grid
remains unresolved. This sequence uses the distributed 5,000 K boundary and
the approved constant-Upsilon continuation below it. Maxwellian rates and
detailed balance retain actual gas temperature. No accuracy is claimed below
the distributed boundary or over the 1,000 to 100,000 K numerical interval.

The 166-level energy selection combines observed lower levels with identified
2014 ICFT target energies. Original `wgfa` rates mix the lowest-20-level MCHF
calculation with the 2014 calculation. Raw wavelength-column photon gaps and
inferred photon gaps from `A/gf` differ from selected level gaps for some
transitions; these differences alone do not identify an incorrect rate. The
source-specific oscillator-strength and frequency conventions require further
review before physical adoption. [The source audit](source_audit.json) records
the original hashes, all first-node temperatures, and the measured gap
differences. The final article was read from the preserved published PDF.

`prepare_plan.py` verifies the isolated Intel source and GNU coefficient-reader
builds, then lists the original CHIANTI files, Stout energy and collision
tables, paper, selection and fitting code, driver objects, executables,
controller scripts, and source audit as frozen inputs. After root review,
`launch_controller.sh` starts one persistent `tmux -L mochii-workflows` session
when `OV_CONDITIONAL_DIAGNOSTIC_APPROVED=yes` is set. It refuses to overwrite
an existing controller or product. The five stages are:

1. Check bound identities, pair uniqueness, directed radiation, and the
   5,000 K source floor.
2. Generate the full source-node atom. Compare compiled source strengths,
   cooling, and all 166 stationary populations on 401 temperatures and four
   densities with independent Python values. Require finite differences below
   `1e-10`; verify exact serialization of all 205,425 original BT nodes.
   Save full arrays, cooling PDFs, and population heatmaps at each density.
3. Verify local and global stationarity and collisional-energy input against
   radiative output from saved native populations. Require gross-flow and
   gross-power-normalized errors below `1e-8`.
4. Fit full-source low-density cooling on the existing 201-point grid. Validate
   its serialized coefficients on 801 temperatures to below 0.3%, including
   the actual GNU coefficient reader, and save a vector PDF plus full arrays.
5. Check the actual Intel registered density table: 64,561 cells, 40 sampled
   comparisons, and interpolation/application errors below 0.3%.

| Input condition | Independent source calculation | Compiled atom calculation |
| --- | --- | --- |
| Gas temperatures | 401 log-spaced values, 1,000 to 100,000 K | Identical values |
| Electron densities | 1e-6, 100, 1e6, 1e10 cm^-3 | Identical values |
| Levels and collisions | Current complete 166-level bound source, 13,695 pairs | Generated 166-level atom, original source nodes |
| Radiation | Current 2,995 selected-source rates | Same serialized rates |
| Geometry, velocity, and spatial profiles | None; statistical equilibrium only | None; statistical equilibrium only |

All stages use one BLAS thread and `I_MPI_FABRICS=shm`, an 8 GiB address-space
limit, and finite CPU/wall limits. The controller records source hashes,
commands, exit status, and checked products. A future physically corrected O V
model needs a separate preserved validation sequence.
