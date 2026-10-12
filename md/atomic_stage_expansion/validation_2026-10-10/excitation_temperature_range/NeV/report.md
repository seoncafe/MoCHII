# NeV excitation temperature-range diagnostic

| Input | Full source and selected atom |
|---|---|
| Temperature | 401 logarithmic samples from 100 to 1,000,000 K |
| Electron density | 1e-06, 100, 1e+06, 1e+10 cm^-3 |
| Levels | Full 304; selected first 25 source indices |
| Radiative transitions | Full 7530; selected 100 |
| Electron collision transitions | Full 46056; selected 300 |
| Potential and boundaries | Identical source atom; no gravitational potential or spatial domain in this statistical equilibrium calculation |
| Spectral input | None; no radiation field |
| Species and reactions | Same ion; spontaneous emission and electron excitation/de-excitation only |

Diagnostic completed: True. Full-domain source cutoff accepted at 0.3%: False.
Numerical floor: 1.42404726944460889e-306 erg cm^3 s^-1. Counts below the floor are recorded in `result.json`; they provide no relative-error evidence.
The full cooling, absolute cutoff, and density-ratio curves are in the three PDF files. All values and numerical masks are in `full_domain_profiles.npz`.
