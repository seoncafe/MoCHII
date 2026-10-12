# Cl V ground-term electron excitation: provisional source audit

The final Saraph and Storey article is *Astronomy and Astrophysics Supplement
Series* **134**, 369–375 (1999), DOI
[10.1051/aas:1999437](https://doi.org/10.1051/aas:1999437), ADS bibcode
`1999A&AS..134..369S`. ADS resolves its publisher PDF to
`https://aas.aanda.org/articles/aas/pdf/1999/02/ds1586.pdf`. That URL and the
publisher landing route returned HTTP 403; the public ADS PDF route also
returned HTTP 403, and the ADS `/full/1999A%26AS..134..369S` route returned
an empty HTTP 202 scan response. The CiteSeerX PDF route returned HTTP 404,
and the ResearchGate full PDF was inaccessible. No matching article PDF was
found in `~/MoCHII/references`. The final journal-formatted text indexed by
[ResearchGate](https://www.researchgate.net/publication/253614097_Atomic_data_from_the_IRON_Project_XXX_Collision_data_for_the_Po_1over_2_-_Po_3over_2_fine-structure_lines_of_P_iii_S_iv_and_Cl_v)
and [CiteSeerX](https://citeseerx.ist.psu.edu/document?doi=15e037d5725168fa42aab91856e8d82d4ab83d20&repid=rep1&type=pdf)
exposes the methods and Table 5. The numbers transcribed here are **pending
verification in the original final PDF**. This audit cannot certify a
published-data import.

The indexed Section 3 describes a 21-state close-coupling calculation for Cl V
and says that scattering energies through 1.5 Ryd are sufficient for rates
through 30,000 K. Indexed Table 5 prints 26 Cl V effective collision strengths
from 1,000 through 40,000 K and stops there. Thus even its final printed row
extends above the temperature the methods explicitly justify. A 1,000–100,000 K
Cl V cooling fit is not supported by this article. A later published calculation
covering 100,000 K was not identified in focused ADS searches for Cl V collision
strengths, aluminum-like chlorine excitation, and in the APAP source index;
this is a search result, not proof that none exists.

The actually invoked MOCASSIN v2 file is
`/home/kiseon/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/data/clv.dat`
(SHA-256 `819052a750275325438e384bed46649b69627b2fc04d86125161d81ded415a52`).
`data/fileNames.dat` lists it, `source/hydro_mod.f90::makeElements` loads it,
and `source/emission_mod.f90::forLines` calls `equilibrium`, which interpolates
its collision strengths versus log temperature and uses them in collisional
rates and line emission. The file contains one two-level collision transition,
15 nodes at 2,000–200,000 K, `g=(2,4)`, a 1490.8 cm^-1 gap, and
`A_21=2.98e-2 s^-1`. Eleven of its nodes coincide with printed Table 5
temperatures between 2,000 and 40,000 K; their strengths agree within 0.005,
the journal's two-decimal printing precision. The largest discrepancy between
the MOCASSIN natural cubic curve and the provisional Table 5 log-log PCHIP
curve over 2,000–40,000 K is 0.511%. The MOCASSIN 50,000, 100,000, and
200,000 K nodes lie beyond the printed table and are **not verified as published
results**. Its `in prep.` header predates the final article.

Cloudy c25 Stout `cl_5.nrg` (SHA-256
`b67cfcc09ca10a35469b0e2893dd90e10520f014cfb089b23388dd5d138bb05a`)
gives the ground fine-structure gap 1490.5 cm^-1 and weights 2 and 4;
`cl_5.tp` (SHA-256
`9348285cd5864d6ae694c19b954016a340ddf154e137d106d66149efe7ae76e9`)
gives `A_21=2.98e-2 s^-1`. The provisional atom uses these energy and
radiative data, so it differs from the MOCASSIN gap by 0.3 cm^-1. The original
MOCASSIN A attribution was unpublished; the independent Stout value supports
the candidate, but the final journal article does not establish that A.

| Input or method | Provisional candidate | MOCASSIN v2 |
| --- | --- | --- |
| Collision strengths | Indexed final Table 5, 26 printed values | `clv.dat`, 15 values |
| Temperature nodes | 1,000–40,000 K | 2,000–200,000 K |
| Interpolation | PCHIP of log Upsilon versus log T | Natural cubic Upsilon versus log10 T |
| Level gap | 1490.5 cm^-1, Cloudy Stout/NIST | 1490.8 cm^-1 |
| Radiative rate | 2.98e-2 s^-1, Cloudy Stout | 2.98e-2 s^-1 |
| Levels and collider | Two ground-term levels; electrons | Same |
| Density normalization | Cooling per electron and Cl V ion | Same two-level rates |

`check_candidate.py` writes `clv_published_table_candidate.txt` in source-node
version 2 format and compares it with the existing optimized native source
driver. It tests 801 temperatures from 1,000 to 40,000 K and electron
densities `1e-6`, `1e2`, `1e6`, and `1e10 cm^-3`. The independent analytic
two-level equation includes Maxwellian excitation and de-excitation, spontaneous
decay, and the actual source temperature. Measured maximum relative native
Upsilon error is `8.88e-16`; maximum relative cooling error is `1.33e-15`.
The [full cooling plot](candidate_cooling_profiles.pdf) was visually inspected;
the dashed native curves overlay the analytic curves throughout the grid.
`result.json` records source and executable hashes, node comparisons, and
numerical checks. Run with:

```sh
MPLCONFIGDIR=/tmp/mplconfig_clv python \
  md/atomic_stage_expansion/validation_2026-10-10/clv_ground_excitation/check_candidate.py
```

This two-level result omits cooling through higher excited Cl V states. It
does not provide a full-ion cooling convergence test. The original final PDF
must be read before claiming the provisional Table 5 values and methods are
verified. A physically supported 100,000 K product additionally needs a
published collision source that covers that temperature; the MOCASSIN
high-temperature nodes alone do not establish one. No production atom,
cooling fit, or registry entry was changed.
