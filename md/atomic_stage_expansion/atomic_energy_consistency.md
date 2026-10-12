# CHIANTI level-energy consistency audit

The CHIANTI level index is an identity, not an energy rank. A level-energy
choice must define one set of physical wavenumbers for Boltzmann factors,
detailed balance, spontaneous emission, and emitted photon energy. The present
observed-when-available choice is inconsistent with some directed `.wgfa` and
`.scups` rows. A nonpositive physical gap cannot be made into an excitation
merely because the file lists the second level index after the first.

`tools/fitting/atomic_level_energies.py` reads the original CHIANTI files and
audits three level-energy choices: observed only, theoretical only, and observed
with theoretical values for missing observations. It keeps reversed and equal
energy pairs visible. It does not change external files or generate cooling
products. Run `python tools/fitting/atomic_level_energies.py n_4 --output
/tmp/n_4_energy_audit.json` for a detailed ion report.

## Source behavior checked

- CHIANTI v11 `pro/read/read_elvlc.pro`, lines 290-294, constructs its `energy`
  field from observed values, replacing a missing observed value with the
  theoretical value. This is a source policy, not a proof that all connected
  transition energies remain physically ordered.
- Cloudy c25.00 `source/atmdat_chianti.cpp`, lines 930-950 and 1115-1210,
  provides experimental, theoretical, and mixed selections. It sorts selected
  energies and remaps original transition indices. Its treatment shows why an
  index-order assumption cannot establish a physical direction. Cloudy's
  theoretical selection falls back to experimental energy when theoretical
  energy is missing; all theoretical energies in these 20 ions are present, so
  this distinction does not affect the measurements below.
- Cloudy Stout `data/stout/n/n_4/n_4.nrg` is a separate energy dataset. Its
  first 29 levels are ordered and include measured values, but its level labels
  and later content cannot be treated as a drop-in replacement for CHIANTI
  `.scups` and `.wgfa` without a verified level correspondence.

## Direct measurements on 20 added ions

The following values come from executing the audit on the raw CHIANTI v11
`.elvlc`, `.wgfa`, and `.scups` files on 2026-10-09. A radiative count includes
only rows with positive Einstein A. Collision counts include all parsed
`.scups` rows. An equal gap is listed separately from a negative gap.

| Energy choice | Missing level energies | Ions with reversed radiative rows | Reversed radiative rows | Equal radiative rows | Ions with nonpositive collision gaps | Reversed collision rows | Equal collision rows |
|---|---:|---:|---:|---:|---:|---:|---:|
| Observed only | 3230 | 2 | 28 | 0 | 14 | 532 | 46 |
| Theoretical only | 0 | 1 | 16 | 3 | 8 | 338 | 12 |
| Observed, then theoretical | 0 | 9 | 134 | 0 | 14 | 2063 | 48 |

For N IV, the observed-then-theoretical choice produces 14 reversed radiative
rows, 81 reversed collision rows, and 17 equal collision rows. The theoretical
choice removes its radiative reversals but retains two equal collision gaps.
For example, N IV levels 107 and 116 have observed energies 595296 and 593690
cm^-1, respectively; their theoretical energies are 588134 and 593223 cm^-1.
The `.wgfa` row 107 to 116 has A = 74600 s^-1 and wavelength 62266.5 Angstrom.
Its observed gap is -1606 cm^-1 and its theoretical gap is +5089 cm^-1. The
positive wavelength corresponds approximately to the magnitude of the observed
gap, but the row's declared downward direction conflicts with that observed
ordering. Swapping its indices would assign A to a different upper state,
which is not justified by the source data.

Choosing theoretical energies everywhere is not a general repair. Fe VII still
has 16 reversed radiative and 327 reversed collision rows; Fe V gains seven
reversed collision rows; Ca VI gains two equal radiative and two equal collision
rows. Equal energy pairs also need explicit treatment because a zero photon
energy and a positive Einstein A do not describe an ordinary radiative line.

## Consequence for generation

Do not silently omit reversed or equal collision transitions and do not change
the direction of an Einstein A row solely from the selected energy ordering.
The 20 complete models require an ion-specific physical decision, supported by
level identification and source wavelengths, before their cooling tables are
accepted. A consistent implementation must then apply the chosen level set to
all transition gaps and to the rate matrix. This audit is diagnostic only; it
does not certify any of the three policies for every ion.

Validation was limited to the new diagnostic and its three focused tests,
including the N IV source row. Cooling products and unrelated regression tests
were not regenerated because this change does not modify their generation path.
