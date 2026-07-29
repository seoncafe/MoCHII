# H/He recombination for MoCHII: the Mao & Kaastra (2016) + Badnell/TAMOC method

Analysis date: 2026-07-14. Status: analysis only (no code changed yet).

## 1. Motivation

MoCHII currently computes the H/He recombination coefficients that enter the
ionization balance from the **Hui & Gnedin (1997)** analytic fits
(`src/recomb_mod.f90`: `alphaA_HII`/`alphaB_HII`, `alphaA_HeII`/`alphaB_HeII`,
`alphaA_HeIII`/`alphaB_HeIII`; He III is the Z=2 hydrogenic scaling of the H
fit). Two things prompt a revisit:

1. **Consistency.** The metals already draw their RR+DR from Badnell 2023
   (`data/atomic/badnell_{rr,dr}.dat`). H/He are the last species still on a
   separate (Hui & Gnedin) source.
2. **The HII40 He over-ionization residual.** MoCHII keeps too little He I
   (He I ~0.16 vs converged MOCASSIN 0.21); a comparison of the H/He
   recombination across codes (see section 5) showed MoCHII's He II -> He I
   coefficient sits ~13% below MOCASSIN's Verner & Ferland value, in the
   direction that would explain part of the residual.

The exact hydrogenic photoionization cross section for H I and He II was just
adopted (`photo_xsec.f90::sigma_hydrogenic`); recombination is the natural
companion, since the ground-state recombination is the detailed-balance (Milne)
partner of exactly that cross section.

## 2. The two resources

### 2a. Mao & Kaastra (2016), A&A 587, A84: level-resolved RR

Path: `~/RT_Codes/CHIANTI/RR_data/J_A+A_587_A84/` (`rr_spex3.dat`, `ReadMe`,
`plot_rr.ipynb`). VizieR J/A+A/587/A84.

Level-resolved **radiative** recombination rate coefficients for the H-like
through Na-like isoelectronic sequences, Z = 1-30, 39696 records. For the
**H-like sequence the coefficients are computed from the exact quantum
hydrogenic photoionization cross sections under detailed balance** (the same
physics as `sigma_hydrogenic`); He-like through Na-like use the ADAS archival
data. Each record is one nlj level (columns `n`, `sp` = 2S+1, `L`, `J`, `cfg`).

The parameterization (T in eV, `T_eV = T[K] * k_B`, `k_B = 8.617333e-5 eV/K`):

```
alpha_level(T) = a0 * 1e-10 * T_eV^(-b0 - c0*ln T_eV)
                          * (1 + a2*T_eV^(-b2)) / (1 + a1*T_eV^(-b1))   [cm^3/s]
```

with seven fitted parameters `a0, b0, c0, a1, b1, a2, b2` per level (fit
accuracy < 5% for ~99% of levels). Because it is **level-resolved**, the case
decomposition is a simple sum:

- `alpha_1` (ground, direct capture to n=1) = the single n=1 level,
- `alpha_A` (case A, total) = sum over all levels,
- `alpha_B` (case B) = `alpha_A - alpha_1` = sum over n >= 2.

The recombined ion labels the sequence: H II -> H I is the H-like sequence
(s=1, z=1); He III -> He II is H-like (s=1, z=2); He II -> He I is He-like
(s=2, z=2). The n=1 ground level in the data is `1s` for the H-like ions and
`1s2` (1S_0) for the He-like ion.

### 2b. Badnell / TAMOC: resolved RR + DR and the cascade method

Path: `~/RT_Codes/CHIANTI/MgII_Badnell/` (a Mg II recombination-line study).
Source: the Strathclyde AMDPP TAMOC database,
`https://amdpp.phys.strath.ac.uk/tamoc/DATA/{RR,DR,PI}/` (also on OPEN-ADAS).

The `nrb##...dat` files are Badnell's resolved recombination data (adf48-style),
providing RR and DR at three resolutions: nlj-resolved (n=3-8), nl-resolved
(n=1-10), and n-resolved (n=1-999 on a non-uniform grid), plus the matching
photoionization. The MgII work (`implementation_plan.md.resolved`,
`run_cascade.py`) computes an **effective recombination coefficient** for a
diagnostic line, e.g. alpha_eff(3p -> 3s) for Mg II 2796/2803, by cascading the
level-resolved RR+DR down through the bound levels:

```
alpha_eff(level, T) = sum_{n,l} [alpha_RR(n,l,T) + alpha_DR(n,l,T)] * P(n,l -> level)
```

using CHIANTI A-values for the low levels and hydrogenic A-values for the high
levels. This is the machinery for **recombination-line** effective coefficients
(the diagnostic side), not the total recombination that closes the ionization
balance; it is the same data ecosystem (Badnell RR + Badnell DR) and is the
right tool if/when MoCHII wants recombination-line emissivities for metals
computed self-consistently rather than from the SH95/Porter tables.

The **total** Badnell RR and DR fits MoCHII already carries
(`data/atomic/badnell_{rr,dr}.dat`, verified md5-identical to the Cloudy c23.01
copies of the TAMOC `clist_K` files) are the converged sums of these resolved
rates:

```
alpha_RR(T) = A / [ sqrt(T/T0) (1+sqrt(T/T0))^(1-B') (1+sqrt(T/T1))^(1+B') ],  B' = B + C*exp(-T2/T)
alpha_DR(T) = T^(-3/2) * sum_i c_i exp(-E_i/T)
```

## 3. Method for the MoCHII ionization balance

The ionization balance needs `alpha_A` (with the diffuse field) or `alpha_B`
(on-the-spot, the MoCHII default). Badnell gives only the **total**
(`alpha_A`); Mao & Kaastra give the level resolution needed to split off the
ground term. The clean combination:

- **`alpha_A` = Badnell RR (+ Badnell DR for He II)** — the converged total,
  including the full high-n tail. (See the caveat below on why not to sum the
  Mao & Kaastra levels for the total.)
- **`alpha_1` = Mao & Kaastra n=1 level** — the ground-state direct capture,
  exact for the H-like ions (H II, He III) and ADAS-based for He II -> He I.
- **`alpha_B` = `alpha_A` - `alpha_1`.**

This uses both resources, is on the modern Badnell/detailed-balance footing
throughout, and reproduces the textbook H values to < 1% (section 4).

**High-n tail caveat.** Summing the tabulated Mao & Kaastra levels alone
undershoots the total: the H-like table stops at n=16, and at 10^4 K the sum
over n <= 16 gives alpha_A = 3.96e-13 vs the converged Badnell 4.19e-13
(~5.6% low). The missing n > 16 levels are a real ~5% of the total. So take
alpha_A from the Badnell converged fit and use Mao & Kaastra only for the
ground level alpha_1 (a single level, no tail issue) and hence the A/B split.

Badnell DR coefficients for H and He: **H II and He III have none** (bare
nuclei, no core electron to excite). **He II -> He I** has a Badnell DR fit but
it is proportional to exp(-46300/T) and is < 1e-29 at 10^4 K — negligible below
~5e4 K, matching MOCASSIN, which omits H/He DR entirely.

## 4. Verification (T = 10^4 K)

Recommended values (alpha_A = Badnell, alpha_1 = Mao & Kaastra n=1,
alpha_B = alpha_A - alpha_1) against the current Hui & Gnedin case B and the
textbook (Osterbrock) values:

| Reaction | alpha_A (Badnell) | alpha_1 (Mao n=1) | alpha_B (new) | alpha_B (H&G, current) | new / H&G | textbook alpha_B |
|---|---|---|---|---|---|---|
| H II -> H I | 4.193e-13 | 1.587e-13 | 2.606e-13 | 2.592e-13 | 1.005 | 2.59e-13 |
| He III -> He II | 2.190e-12 | 6.537e-13 | 1.536e-12 | 1.545e-12 | 0.994 | ~1.55e-12 |
| He II -> He I | 4.372e-13 | 1.565e-13 | 2.807e-13 | 2.616e-13 | 1.073 | ~2.73e-13 |

Cross-checks:
- alpha_1(H II) = 1.587e-13 matches the textbook ground-state H recombination
  1.58e-13 to 0.4%; alpha_B(H II) = 2.606e-13 matches 2.59e-13 to 0.6%.
- **A from-scratch Milne integral written during this analysis gave
  alpha_1(H) = 7.9e-14 — a factor-of-two low (a statistical-weight error).**
  The published Mao & Kaastra parameterization avoids exactly this class of
  mistake; use it rather than re-deriving the Milne integral.

## 5. Cross-code comparison (why this matters)

Recombination in the ionization balance, at 10^4 K, cm^3/s:

| Reaction | MoCHII/EXHALE (Hui & Gnedin 97) | MOCASSIN (Verner & Ferland 96, case A) | Cloudy (Badnell 06/23) | This method |
|---|---|---|---|---|
| H II -> H I | aA 4.30e-13 / aB 2.59e-13 | 4.19e-13 | 4.19e-13 | aA 4.19 / aB 2.61e-13 |
| He II -> He I | aA 4.23e-13 / aB 2.62e-13 | 4.84e-13 | 4.37e-13 | aA 4.37 / aB 2.81e-13 |
| He III -> He II | aA 2.23e-12 / aB 1.55e-12 | 2.19e-12 | 2.19e-12 | aA 2.19 / aB 1.54e-12 |

Code sources: MoCHII `recomb_mod.f90`; EXHALE `Cool_coeff.f90` (identical Hui &
Gnedin coefficients, case B only); MOCASSIN `update_mod.f90` -> `radrec.dat`
(Verner & Ferland 1996 total, case A + explicit diffuse field, no H/He DR);
Cloudy H-like/He-like iso-sequence solver (total consistent with Badnell RR +
Badnell DR, case A + diffuse).

> **Note (2026-07-17): EXHALE default moved to Badnell + Mao & Kaastra.**
> The "MoCHII/EXHALE (Hui & Gnedin 97)" column and the "EXHALE `Cool_coeff.f90`
> (identical Hui & Gnedin coefficients)" description above are the historical
> shared baseline. On 2026-07-17 EXHALE switched its default H/He case-B
> recombination to the same Badnell (2023) RR (+ He II DR) minus Mao & Kaastra
> (2016) alpha_1 construction as MoCHII (the "This method" column: aB 2.61e-13
> H II, 2.81e-13 He II, 1.54e-12 He III at 10^4 K), so the two codes now share
> the default coefficients. The Hui & Gnedin fits remain reproducible in EXHALE
> via `legacy_hhe_rates`.

Findings:
- **H II and He III (hydrogenic): the new method equals the current Hui &
  Gnedin to within 0.5%.** Hui & Gnedin's H fit was already excellent;
  switching changes essentially nothing for H and He III.
- **He II -> He I is the one real change: +7.3%** (alpha_B 2.62e-13 ->
  2.81e-13). Higher He recombination -> more He I -> moves MoCHII toward
  MOCASSIN's higher He I fraction. This is the correct direction for the HII40
  He over-ionization residual, though only part of it (MOCASSIN's Verner &
  Ferland He value, 4.84e-13, is itself the high outlier; the modern Badnell
  4.37e-13 sits between Hui & Gnedin and Verner & Ferland).
- The new He II -> He I alpha_B = 2.81e-13 is close to the Hummer & Storey
  textbook ~2.73e-13.

## 6. Proposed implementation for MoCHII

Follow the established atomic-data convention (fit offline, evaluate closed
forms online):

1. **Ground term alpha_1.** Add a data file, e.g.
   `data/atomic/rr_ground_HHe.txt`, holding the Mao & Kaastra n=1 parameters
   for the three ions (H-like z=1, H-like z=2, He-like z=2), with a provenance
   header pointing at J/A+A/587/A84 `rr_spex3.dat`. `recomb_mod` evaluates the
   Mao & Kaastra closed form for alpha_1(T).
2. **Total alpha_A.** Evaluate the Badnell RR (+ He II DR) closed form in
   `recomb_mod` from the three coefficient sets already verified here (H II
   Z1N0, He III Z2N0, He II Z2N1 of `badnell_{rr,dr}.dat`). These can be
   hardcoded constants with a provenance comment (as the VFKY96 rows are in
   `photo_xsec`), or read from the existing Badnell files.
3. **alpha_B = alpha_A - alpha_1.** The public interface
   (`alphaA_*`/`alphaB_*`, elemental in T) is unchanged, so `ion_balance_mod`
   and `diffuse_mod` (which forms `alpha_1 = alpha_A - alpha_B` for the diffuse
   packets) need no change — the diffuse ground-recombination rate becomes
   exactly the Mao & Kaastra alpha_1, which is more correct than the
   Hui & Gnedin difference.
4. **Keep Hui & Gnedin available** behind a switch (e.g.
   `par%recomb_model = 'badnell_mao' | 'hui_gnedin'`, default the new one) so
   the recorded gates can be reproduced and the He II swap can be A/B tested in
   place, exactly as `par%ci_model` did for collisional ionization.

Gate impact expected: H and He III are unchanged to < 0.5%, so the Stromgren
and HII20 gates should move negligibly; HII40 should show a small reduction in
He (and metal) over-ionization from the +7.3% He II -> He I recombination.
Re-baseline together with the pending `ion_align_edges` rebaseline.

## 7. The Badnell cascade method (future, separate)

The MgII_Badnell cascade path (section 2b) is the tool for **effective
recombination-line coefficients** computed self-consistently from Badnell
resolved RR+DR, rather than the total needed for the balance. It is the natural
route if MoCHII later wants metal recombination lines beyond the current
SH95 (H I), He II SH95, and Porter He I tables — for instance recombination
contributions to [O III]/[N II]/C II lines, or the Mg II 2796/2803 effective
coefficient. Out of scope for the ionization-balance change above; recorded
here because it is the same data source and was the method the analysis drew
from.

## 8. Open items / caveats

- The He II -> He I alpha_1 uses the He-like ADAS data (not exact hydrogenic);
  it is the ground 1s2 1S_0 direct capture. Case B for He is more intricate
  than for H (singlet/triplet, the 2^3S metastable); alpha_B = alpha_A -
  alpha_1(1s2) is the on-the-spot approximation MoCHII already makes and is
  adequate at the current level.
- alpha_A from Badnell and alpha_1 from Mao & Kaastra are two sources; they are
  mutually consistent (both converge to the textbook totals), and alpha_B is
  positive and smooth over 5e3-2e4 K, but this is a two-source construction by
  design (exactly the Cloudy decomposition: total minus resolved ground).
- Temperatures checked: 5000, 10000, 20000 K. The Mao & Kaastra fits are valid
  1e3-1e6 K (10^-1 to 10^2 eV region well covered); the Badnell fits span the
  full nebular range.

## Data and script paths

- Mao & Kaastra: `~/RT_Codes/CHIANTI/RR_data/J_A+A_587_A84/rr_spex3.dat`
- Badnell RR/DR (in MoCHII): `data/atomic/badnell_{rr,dr}.dat`; source
  `https://amdpp.phys.strath.ac.uk/tamoc/DATA/{RR,DR}/`
- Badnell resolved + cascade example: `~/RT_Codes/CHIANTI/MgII_Badnell/`
- Analysis scripts (session scratchpad, not committed): `mao.py`, `final.py`,
  `rec_cmp.py`, `milne.py`
