# Conditional full-bound O IV numerical validation

This directory prepares a private, source-hash-frozen numerical sequence for
the current complete O IV bound atom: 147 levels, 2,319 radiative rows, 10,731
collision pairs, and 160,965 original scaled collision nodes. The previous
35-level atom is not an accepted complete-source installation candidate. The
current 147-level radiative-strength provenance remains under independent
physical review. Every result in this sequence has
`physical_model_accepted=false`, and no product is installed automatically.

The 1,000–100,000 K test range is a numerical comparison interval, **not** a
published atomic-accuracy interval. The present numerical model has no accepted
full radiative-strength provenance. Source-level energy selection combines
observed valence states with original higher theoretical target energies.
The published collision calculation covers 3,200–32,000,000 K. Below 3,200 K,
the adopted approximation holds the collision strength at its first valid
temperature while retaining actual gas temperature in Maxwellian rates and
detailed balance. The 1,000–3,200 K numerical tests therefore check the
implemented approximation, not published collision accuracy.
Unshifted scattering resonances and the consistency of all radiative line
strengths with those selected energies remain physical uncertainties. Native
agreement, stationary balance, fit accuracy, and density-table accuracy do not
resolve those source questions.

`prepare_plan.py` freezes the current CHIANTI source files, ionization
potentials, selection and fitting code, controller scripts, source and GNU
coefficient-reader executables, Intel n-level registry driver, linked objects,
and build manifests. It verifies the Intel source/objects against the isolated
build record. The GNU coefficient check uses `cooling_mod` and its driver;
the current Intel registry independently checks the n-level module. All source
file hashes, commands, exit statuses, and step transitions are retained by
`~/.codex/scripts/calculation_workflow.py` under `controller/`.

After root review, the persistent `launch_controller.sh` requires
`OIV_CONDITIONAL_DIAGNOSTIC_APPROVED=yes`. It emits one plan and starts tmux
session `mochii_oiv_full_20261010`. It refuses to overwrite an existing plan,
state directory, or product. The prepared stages are:

1. Check complete bound identities, pair uniqueness, and radiation direction.
2. Generate the source-node atom; compare compiled source strengths and
   401-temperature full populations/cooling at four electron densities with
   independent Python source calculations. Require finite errors below
   `1e-10` for strengths, cooling, and populations. Verify all 160,965 original
   scaled collision nodes round-trip exactly. Preserve full profiles, source
   arrays, PDFs, an all-147-level population heatmap at every density, and atom
   hashes.
3. Form an independent level-stationarity and collisional-energy versus
   radiative-energy ledger from the saved native populations and current
   original collision nodes. Require gross-flow-normalized stationarity and
   energy residuals below `1e-8`.
4. Fit the full low-density O IV cooling source on the existing 201-point
   fitting grid, validate the serialized fit on an independent 801-point
   grid with maximum relative difference below 0.3%, and check the actual GNU
   coefficient reader. Preserve the full source/fit/difference arrays and PDF.
5. Read the generated atom and fit with the actual optimized Intel registered
   cooling driver. Require 64,561 table cells, 40 comparison samples, and
   relative interpolation/application differences below 0.3%.

The sequence uses one BLAS thread, an 8 GiB address-space limit, and finite
stage CPU/wall limits of 180, 2,400, 900, 3,600, and 1,800 seconds. A later
physical source correction requires a separate preserved validation sequence,
not an overwrite of this conditional diagnostic.

| Input condition | Independent source calculation | Compiled atom calculation |
| --- | --- | --- |
| Gas temperatures | 401 log-spaced values, 1,000–100,000 K | Identical values |
| Electron densities | 1e-6, 100, 1e6, 1e10 cm^-3 | Identical values |
| Levels and collisions | Current complete 147-level bound source, 10,731 pairs | Generated 147-level atom, original source nodes |
| Radiation | Current 2,319 selected-source rates | Same serialized rates |
| Geometry, velocity, and spatial profiles | None; statistical equilibrium only | None; statistical equilibrium only |
