# Mg I retained in the extended excitation selection

## Confirmed defect and correction

The selected production atom path is constructed by
`metal_excitation_data_file` in `src/nlevel_mod.f90`. Its thermal consumers
in `species_setup`/`nlevel_cooling_add` and its line consumer in `lines_write`
read that same directory. Before this correction, the approved Mg I atom
existed in `mocassin_v2` but was absent from `extended`. Selecting `extended`
therefore lost Mg I excitation cooling and emitted lines; the existence of
the reference file did not make those selected consumers read it.

The identical approved five-level Mg I atom was added to `extended`, with
its manifest entry and absolute-cooling-without-fit record. The installation
generator now includes it in both selections. Its constant collision strengths
retain the existing 1000-100000 K model interval and explicitly unverified
atomic accuracy. No new physical approximation was adopted. The previous
manifest is preserved in `previous_manifest.json`.

## Validation

The existing namelist/model-selection diagnostic was compiled against the
current production objects and executed with explicit `extended` input.
The invoked selector and production reader loaded 62 atoms with 1951 levels,
including Mg I with 5 levels, 10 radiation rows, and 10 collision rows.
The two Mg I files are identical and their installation hashes are recorded.
`reader_result.json` records the executed source, linked objects, executable,
input, and scope. `installation.json` records the data correction.

The new read-only selected-stage inventory confirms that Ca VII is now the
only supplied model missing in either selection. There are no selected
radiation rows with nonpositive photon gaps. Counts do not establish state
identity or atomic accuracy. No cooling table, fit, or numerical profile was
regenerated, since this change adds the same accepted data at the selected
path. No full regression, thermal-equilibrium calculation, or spatial
comparison was performed. The overall 77-stage goal remains incomplete.
