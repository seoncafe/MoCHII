# PN150 L(Hbeta) residual against Cloudy c25.00 (2026-10-14)

MoCHII `../pn150_p27_nodrift` gives L(Hbeta) 0.893 of Cloudy's. The
measurements here localize it to metal recombination radiation; the
findings and the plan are in
`md/METAL_RECOMBINATION_DIFFUSE_PLAN_2026-10-14.md`.

| script | what |
|---|---|
| `budget.py` | MoCHII vs Cloudy: ionized volume, EM, R_B(H), Hbeta per recombination, He recombinations |
| `budget_runs.py <base>...` | R_B per Q(H) for MoCHII runs |
| `metal_absorption.py`, `absorption_runs.py <rates.h5>` | photoionizations by absorber (MoCHII) |
| `cloudy_absorption.py` | the same for Cloudy, from `gam.gam` (save gammas) |
| `euv_lines.py` | MoCHII metal lines shortward of 911.75 A |
| `cloudy_variants.py <base>...` | R_B per Q, L(Hbeta), printed EUV metal lines for Cloudy decks |

Cloudy decks (`*.in`, run with `~/CLOUDY/c25.00/source/cloudy.exe -r <name>`):
`base`/`faint` (the reference deck), `iter` (iterate to convergence),
`nometal`/`nometal_f` (C N O Ne Mg Si S off), `nocx`, `noind`, `nolev2`,
`noots`, `nooutl`, `nooutc`, `nometal_noots`, `nometal_nooutc`, `gam`.
MoCHII L6 runs: `l6_metals.in`, `l6_nometals.in` (`run_l6.sh`, grid
`pn150_L6.fits`, binary `MoCHII_l6.x`, the tree of 2026-10-13).
