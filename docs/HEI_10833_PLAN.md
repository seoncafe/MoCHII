# He I 10830 Å metastable (2$^3$S) diagnostic — plan

Date: 2026-07-24 (rev 2; status updates same day).  2026-07-26: the two open
gates H2b and H4 spelled out in Section 6 (Cloudy version corrected to
c25.00 with the decks and the save command, what H2b tests and how to
attribute a disagreement; the lifetime argument, the EXHALE entry points, and
the remaining work for H4) and their order in Section 8.
Status: **Stages 1 AND 2 IMPLEMENTED AND GATED.**  Stage 1:
`src/hei_metastable_mod.f90`, `par%hei_metastable`,
`data/atomic/hei_metastable.txt` via `tools/fitting/make_hei_metastable.py`.
Stage 2 (pure Python + LINE10833 data rows): `HeiMetaData` in
`tools/python/mochii_output.py` with `column_map` / `tau10833_map` /
`kappa0_10833` (NIST f/A values adopted; CHIANTI A=1.15e7 rejected).
Gates in `tests/hei_meta/` — see the results notes in Section 6.
H2b (Cloudy anchor) and H4 (EXHALE advection reference) remain open.
Rev 1 superseded the original draft, which was built on a double-counting
error (adding a collisional 10830 row on top of the Porter table; Section 1;
see `docs/MOCHII_PLANE_PARALLEL_HEI10833_REVIEW.md` Section 3 for the review
that caught it).  Rev 2 makes the plan concrete: the primary data source is
the **EXHALE triplet implementation** (same author, already verified, newer
than OH18 in places), the OH18 rate labels are corrected, and the gates gain
analytic limits and an EXHALE advection reference.

Goal: a metastable 2$^3$S (triplet) diagnostic for He I whose primary
observables are the 2$^3$S number density $n_3$, the metastable column
density, the 2$^3$S photoionization rate $\Phi_3$, and the 10830 Å
absorption opacity / optical depth.  These are the quantities that matter
where the Porter case-B emissivity assumptions break down (exoplanet transit
spectroscopy, chromospheres, PDRs) — not extra optically-thin emission.

## 1. What MoCHII already outputs, and the error this plan corrects

MoCHII writes the He I 10830 Å line as `HeI10830` from the Porter et al.
(2012/2013) case-B table `data/atomic/hei_porter_caseB.txt` (parsed from
Cloudy c23.01's `he1_case_b.dat` by `tools/fitting/make_hei_lines.py`).

**That table is a case-B collisional-radiative TOTAL emissivity, not a pure
recombination cascade.**  Verified directly: the tabulated
$4\pi j / (n_e\,n_{\rm He^+})$ for 10830 at $T_e = 10^4$ K rises by a factor
of **6.6** from $n_e = 10$ cm$^{-3}$ ($2.81\times10^{-25}$) to
$n_e = 10^4$ cm$^{-3}$ ($1.85\times10^{-24}$).  A pure recombination cascade
coefficient is nearly density-independent; a 6.6x rise with $n_e$ is the
signature of collisional excitation feeding the upper level.  So the Porter
case-B table ALREADY includes collisional excitation out of the metastable
2$^3$S reservoir into 2$^3$P.  (The steepest part of that rise sits between
$n_e = 10^{2}$ and $10^{4}$ cm$^{-3}$ — exactly where the analytic critical
density $n_{\rm crit} = A_{31}/(q_{31a}+q_{31b})$ of Section 3 falls, an
internal consistency check on the interpretation.)

**Consequence:** adding a separate collisional 10830 row and SUMMING it with
the Porter 10830 — the original draft's central idea, by analogy with the
`hc` collisional H I Balmer rows — would double-count and is wrong.  The H
analogy fails for a specific physical reason: SH95 H I is genuinely
recombination-only (the H ground state is not metastable), so a collisional
H I component weighted by $n_{\rm HI}$ is a distinct channel that
legitimately adds.  Porter He I already carries the metastable-collisional
term; under case-B conditions there is nothing separate left to add.

## 2. The actual gap: a metastable 2$^3$S diagnostic

The Porter case-B emissivity is correct in the dense, optically-thin,
recombination-dominated H II regime it was built for.  It does NOT give:

- the 2$^3$S population $n_3$ and its column density along a line of sight —
  the primary observable for exoplanet transit and chromospheric 10830 work
  (10830 is seen there in ABSORPTION against the stellar disk);
- the 10830 line-center absorption opacity and optical depth, which depend
  on $n_3$ and the line profile, not on any emissivity;
- the 2$^3$S photoionization rate $\Phi_3$, which depletes the metastable
  reservoir wherever a soft-UV field is present.

The goal is therefore a metastable-population module, not an emission row.

## 3. Physics: the 2$^3$S balance (OH18, labels corrected)

Oklopcic & Hirata (2018, ApJ 855, L11 — OH18) track three He states: the
neutral singlet ground $f_1$ (1$^1$S), the metastable triplet $f_3$
(2$^3$S), and the ion.  The triplet balance is

```
  source:  alpha_3 * n_e * n(He+)        (recombination into the triplet system)
         + q_13   * n_e * n(He 1^1S)     (e-impact excitation 1^1S -> 2^3S)

  sink:    n_3 * [ A_31                  (radiative 2^3S -> 1^1S, A_31 = 1.272e-4 s^-1)
                 + Phi_3                 (2^3S photoionization, threshold 4.78 eV)
                 + n_e * (q_31a + q_31b) (e-impact 2^3S -> 2^1S and 2^3S -> 2^1P)
                 + n_HI * Q_31 ]         (Penning/associative ionization with H0)
```

**Label corrections relative to rev 1** (verified against the EXHALE
implementation, `EXHALE/src/modules/radiation/util_ion_eq.f90::HeITR_coeffs`):

- $q_{31a}$ is 2$^3$S $\to$ 2$^1$S and $q_{31b}$ is 2$^3$S $\to$ **2$^1$P**
  (spin-changing e-impact transfers into the singlet system, which then
  decays to 1$^1$S).  Rev 1 mislabeled them as "$\to 1^1$S, $\to 2^1$S";
  the direct e-impact 2$^3$S $\to$ 1$^1$S channel is not the OH18 sink.
- The H-collision term uses the NEUTRAL hydrogen density $n_{\rm HI}$, and
  $Q_{31}$ is the sum of Penning and associative ionization
  (He 2$^3$S + H$^0$ $\to$ products), not a de-excitation.
- $A_{31} = 1.272\times10^{-4}$ s$^{-1}$ (2$^3$S lifetime ~7870 s).
  **Pitfall found during this planning pass**: CHIANTI `he_1.wgfa` lists
  $A(1\!-\!2) = 1.73\times10^{-4}$ — ~36% above the accepted value.  Do NOT
  take $A_{31}$ from the CHIANTI file; use 1.272e-4 (as EXHALE does).

Because $f_3 \ll f_1$ always ($n_3/n_{\rm He} \sim 10^{-7}$–$10^{-6}$),
$n(1^1{\rm S}) \simeq n_{\rm HeI}$ and the balance is LINEAR in $n_3$ — a
closed-form quotient in each leaf, no iteration:

```
  n_3 = [ alpha_3 n_e n_HeII + q_13 n_e n_HeI ]
        / [ A_31 + Phi_3 + n_e (q_31a + q_31b) + n_HI Q_31 ]
```

Two analytic limits anchor the gates (Section 6):

- low-density radiative limit ($A_{31}$ dominates):
  $n_3/n_{\rm He^+} = \alpha_3 n_e / A_{31}$
  ($= 1.65\times10^{-7}$ at $n_e = 100$, $T_e = 10^4$ K);
- high-density collisional plateau: $n_3/n_{\rm He^+} =
  \alpha_3/(q_{31a}+q_{31b})$, density-independent; the crossover
  $n_{\rm crit} = A_{31}/(q_{31a}+q_{31b})$ lands at a few
  $\times 10^3$ cm$^{-3}$, matching where the Porter 10830 table steepens.

**Advection caveat (unchanged from rev 1).**  OH18 solve
$v\,\mathrm{d}f_3/\mathrm{d}r = \rm source - sink$ along an outflowing
hot-Jupiter profile.  A leaf-by-leaf `source = sink` solve is a static /
local-equilibrium approximation: valid in dense H II regions (2$^3$S
destruction time short against any flow time), NOT a faithful reproduction
of the OH18 benchmark.  MoCHII's generic AMR reader
(`src/read_generic_amr.f90`) discards velocity input, so the exoplanet gate
is rescoped to local equilibrium (Section 6, H4) unless an advection solver
is added.  EXHALE, which evolves the triplet as a species in its 1D flow,
is the natural in-house reference for quantifying the advection error.

### 3.1 The 2$^3$S photoionization rate $\Phi_3$

```
  Phi_3 = 4 pi * integral J_nu * sigma_3(nu) / (h nu) dnu
```

with threshold 4.78 eV (2593 Å).  Two concrete points beyond rev 1:

- **Evaluate over the FULL band, not just the FUV segment.**
  $\sigma_3(E)$ remains significant far above 13.6 eV (Cooper minimum near
  34.7 eV, secondary bump near 45.6 eV, slow high-E tail — see the EXHALE
  fit below), so the integral runs over all `nnu_band` bins.  The reduced
  tally `jt_ion(inu, il)` already provides $J_\nu$ in each leaf and bin; the
  computation is the exact pattern of `gas_rates_mod`'s
  $\Gamma_i = \sum_\nu 4\pi J_\nu \sigma_i \Delta\nu / (h\nu)$ — reuse it.
- **Band configuration**: the soft part needs `add_fuv` with
  `efuv_min <= 4.78` eV.  The FUV segment is log-uniform (edge alignment
  currently applies to the ionizing segment only), so simply set
  `efuv_min = 4.78` to pin the threshold at the segment edge; extending
  `ion_align_edges` to FUV-segment thresholds is an optional refinement.

**Correction found at implementation time (2026-07-24)**: the original
"$\Phi_3 = 0$ without `add_fuv`" framing was wrong — the ionizing band
itself (13.6–100 eV) photoionizes 2$^3$S ($\sigma_3 \sim 1$–$2\times
10^{-18}$ cm$^2$ across the band), so $\Phi_3 > 0$ wherever any ionizing
field exists.  The module therefore ALWAYS integrates over the full band;
without `add_fuv` it only misses the soft-UV part $[4.78, {\rm floor})$
(rank-0 note).  In a dense H II core the EUV-only $\Phi_3$ is small
($\Phi_3/A_{31} \lesssim 10^{-2}$ in the gate run); the soft-UV segment is
what makes $\Phi_3$ dominant in the exoplanet / PDR regime.

## 4. Data: port the EXHALE triplet coefficients

EXHALE (`~/RT_Codes/ExoAtmosphere/EXHALE`) already implements the complete
OH18 triplet chemistry, verified in production, with two upgrades over the
raw OH18 choices.  MoCHII ports these expressions (offline into
`data/atomic/` closed forms with provenance headers, per convention —
Fortran only evaluates):

| quantity | EXHALE routine (src/modules/) | expression / source |
|---|---|---|
| $\alpha_3(T)$ | `radiation/Cool_coeff.f90::rec_HeII_23S` | $2.10\times10^{-13}(T/10^4)^{-0.778}$ (OH18) |
| $\alpha_1(T)$ (1$^1$S channel) | `Cool_coeff.f90::rec_HeII_11S` | singlet-channel counterpart (for photon bookkeeping) |
| $q_{13}(T)$ | `Cool_coeff.f90::coex_HeI_1S_23S` | $2.10\times10^{-8}\sqrt{13.6/kT}\,e^{-19.81/kT}\,\Upsilon(T)$ (Bray et al. 2000 via OH18) |
| $q_{31a}, q_{31b}(T)$ | `coex_HeI_23S_21S`, `coex_HeI_23S_21P` | Bray et al. 2000 via OH18 |
| $A_{31}$ | `util_ion_eq.f90::HeITR_coeffs` | $1.272\times10^{-4}$ s$^{-1}$ |
| $Q_{31}(T)$ | `Cool_coeff.f90::penning_HeI_23S` | **Taylor et al. (2025, ApJ 989, 68) Table 2** two-branch fit — see the verification note below |
| $\sigma_3(E)$ | `functions/cross_sec.f90::sigma_HeI23S` | **two VFKY96 wings + log-linear bridge**, ~4% fit to Norcross (1971) + TOPbase, threshold 4.78 eV through the high-E tail |
| $q(2^3{\rm S}\to2^3{\rm P})$ | `Cool_coeff.f90::coex_rate_HeI23S_10833` | $1.16\times10^{-20}\sqrt{T}\,e^{-13179/T}$ (for the optional emission/absorption source function) |

Notes:

- **$Q_{31}$ verified against the papers directly (2026-07-24).**  Taylor
  et al. (2025, ApJ 989, 68; DOI 10.3847/1538-4357/ade3c9) is real and its
  Table 2 fit matches the EXHALE implementation exactly:
  $1.9\times10^{-9}(300/T)^{0.07}$ for $T \le 4000$ K,
  $9.1\times10^{-9}(300/T)^{0.50}$ above — a Maxwell-Boltzmann average of
  the Morgner & Niehaus (1979) + Cohen & Lane (1971) cross sections,
  superseding the OH18/Lampon constant $5\times10^{-10}$ (Roberge &
  Dalgarno 1982) by a factor ~3.  BUT it is not the only 2025 revision:
  **Garcia Munoz (2025, A&A 698, A199)** independently refit the same
  channel two months earlier from the NEWER Movre & Meyer (1997) cross
  sections ($k = 10^{-9}\exp(c/T + d_1\ln T + d_2\ln^2 T + d_3\ln^3 T)$,
  <2% over 200–10000 K); neither paper cites the other.  Numerically
  Taylor is 1.2–2.2x ABOVE Garcia Munoz over 300–10$^4$ K (both agree
  RD82 is ~3x low; the cross-section sets themselves are mutually
  consistent per GM25 Fig. 5, so the residual spread is fit/averaging
  choice).  Taylor's two-branch form is also DISCONTINUOUS at 4000 K
  (jump ~1.6e-9 $\to$ 2.5e-9).  Decision: port the Taylor fit (EXHALE
  consistency), record the GM25 fit as the alternative in the data-file
  header, and treat the ~1.2–2x spread as the documented systematic
  uncertainty of the $n_{\rm HI} Q_{31}$ sink (it matters only in the
  neutral/PDR zone; an H2 gate variant swapping the two fits bounds the
  effect).  ADS search (2025–2026) finds nothing newer for this channel.
- $\sigma_3$ being in VFKY96 form is a gift: MoCHII's `photo_xsec` already
  evaluates `sigma_VFKY96`, so the data file carries two standard 8-parameter
  lines; only the two-wing bridge (a log-linear join between 34.7 and
  45.6 eV) is new logic, a few lines.
- **Consistency with the existing code**: `diffuse_mod` already uses
  `HEI_F3S = 0.75` for the He I diffuse triplet branch.  At $10^4$ K,
  $\alpha_3 = 2.10\times10^{-13}$ vs
  $0.75\,\alpha_B({\rm He\,II}) = 0.75\times2.807\times10^{-13} =
  2.11\times10^{-13}$ (badnell_milne) — consistent to 0.2%.  Document this;
  no re-derivation needed.
- **Draine (2026, ApJ 999, 3) cross-check (2026-07-24)** — resonant
  scattering of the 1.0833 um triplet in H II regions, directly this
  diagnostic's regime.  Population physics AGREES: same balance (his
  Eq.~1; $\Phi_3$/Penning noted subdominant, matching our numbers),
  $\alpha_{\rm trip}$ within 6% (Del Zanna & Storey 2022 vs OH18),
  $k_d$ within ~3% (both Bray 2000), $n_{\rm crit}$ 3855 vs our 3945,
  and his Eq.~4 column reproduces our low-density example center-to-edge
  $N(2^3{\rm S})$ to ~25%.  His $\tau$ formula, however, is exactly
  $3\times$ ours ($\tau_{\rm tot} = 262$ vs our 87.5 for each
  $10^{14}$ cm$^{-2}$ at $b=10$ km/s): his Table 1 lists $g_u = 3, 9, 15$
  (sum 27) for 2$^3$P$^o_J$ where the true degeneracies are $2J+1 = 1, 3,
  5$ (sum 9), inflating Eq.~8 by 3.  The NIST $A = 1.0216\times10^7$ and
  $f = 0.0599/0.1797/0.2996$ are mutually consistent ONLY with
  $g_u = 2J+1$ (verified numerically), and our $\sigma_0$ matches the
  p-winds/OH18 standard validated against transit observations — so the
  factor 3 is on the paper's side (likely erratum material), and MoCHII's
  $\tau_0$ stands.  His qualitative conclusions (broad multipeaked
  profiles, blueshift, dust-absorption enhancement at $\tau \gg 1$)
  survive since H II regions are deep in the thick regime either way.
- **Cross-checks**: (i) p-winds (in-tree at `~/RT_Codes/ExoAtmosphere/p-winds`,
  the public OH18 implementation — `p_winds/microphysics.py::
  helium_triplet_cross_section`, `p_winds/helium.py`) for $\sigma_3$ and the
  rate set; (ii) CHIANTI `he_1` (its `scups` IS Bray et al. 2000, with the
  needed transitions 1–2, 2–3, 2–4/5/6 all present) through the Tier-2
  pipeline for the $\Upsilon(T)$ fits; (iii) NIST ASD for the 10830
  fine-structure $A$/$f$-values — CHIANTI `he_1.wgfa` lists
  $A(10830) = 1.15\times10^7$ s$^{-1}$ where NIST gives
  $\simeq 1.02\times10^7$; resolve before Stage 2 uses them.

## 5. Design (two stages)

- **Stage 1 — metastable population** (`hei_metastable_mod`): evaluate the
  Section 3 closed form in each leaf on the converged state.  Inputs all
  available: $n_e$, $T_e$, $n_{\rm HeII} = x_{\rm HeII} n_{\rm He}$,
  $n_{\rm HeI}$, $n_{\rm HI}$; $\Phi_3$ from the full band tally
  (Section 3.1) — with a rank-0 note when the band floor sits above
  4.78 eV, i.e. the soft-UV part is missed.  Runs once at
  output time (a leaf loop like the other diagnostics), no hot-loop cost,
  no feedback into the ionization/thermal balance ($n_3$ is a
  $\sim10^{-7}$ trace of He; its energy budget is negligible).
- **Stage 2 — 10830 absorption diagnostic**: from $n_3$, the line-center
  opacity of each fine-structure component $i$,
  $\kappa_{0,i} = \frac{\sqrt{\pi} e^2}{m_e c}\,f_i\,\frac{\lambda_i}{b}\,n_3$,
  $b = \sqrt{2kT_e/m_{\rm He} + v_{\rm turb}^2}$ (components at vacuum
  10832.06 / 10833.22 / 10833.31 Å; $f_i$ from NIST after the
  cross-check above).  Outputs: an `n_2s3` leaf block (emis-file section
  under `emis_output`) and a line-of-sight column map — the Python side
  (`tools/python/mochii_output.py` + `leaf_field.py`) gains a small
  `column_map` (integral of a leaf field along the LOS, cm$^{-2}$), the
  flux-conserving deposit reused with $n\,dl$ weights.  $\tau_{10830}$
  follows from the column and $b$; warn when line-center $\tau > 1$ (the
  optically-thin assumption fails; line transfer is out of scope).

If explicit metastable 10830 EMISSION is ever wanted, it must REPLACE the
Porter 10830 row, never be summed with it; the recombination part then needs
a collisionless recombination-into-2$^3$P coefficient DISTINCT from the
Porter case-B total (the two describe the same photons).  The
`coex_rate_HeI23S_10833` rate above provides the collisional part.

## 6. Validation gates

**Stage-1 results (2026-07-24, `tests/hei_meta/check_hei_meta.py` — all
PASS, independently re-run and reproduced bit-identically with a fresh
build):** analytic limits from the data file — low-density ratio
$1.6509\times10^{-7}$ (spec $1.65\times10^{-7}$), plateau
$6.513\times10^{-6}$, $n_{\rm crit} = 3945$ cm$^{-3}$ ($n_e$-scan
half-plateau crossing at 3926); module $n_{2s3}$ vs the independent Python
closed form on each run's own state: worst rel $4.1\times10^{-15}$
(machine precision, all three runs); low-density run core
$n_3/n_{\rm HeII}$ within median 2.8% / max 3.9% of
$\alpha_3 n_e/A_{31}$; high-density run ($n_e = 1.1\times10^4 >
n_{\rm crit}$) at 0.676 of the plateau; $\Phi_3$ (add_fuv,
`efuv_min=4.78`) vs an independent $J_\nu$ integration: worst rel
$1.5\times10^{-9}$, max $\Phi_3 = 9.07\times10^{-5}$ s$^{-1}$;
$\alpha_3(10^4\,{\rm K}) = 0.998\times0.75\,\alpha_B({\rm He\,II})$
(HEI_F3S consistency); $\sigma_3$ vs p-winds (Norcross) within 1–7%.

**Stage-2 results (2026-07-24, same checker — all PASS):** LINE10833
constants carry $f$ = 0.0599/0.1797/0.2996 (1:3:5) and $A =
1.0216\times10^{7}$ s$^{-1}$ each (NIST ASD / Drake 2007, confirmed by
p-winds; the CHIANTI `he_1.wgfa` $A \simeq 1.15\times10^{7}$, ~13% high,
is NOT used — the Section-4 cross-check item is resolved); synthetic
uniform cube: column conservation $\sum N A_{\rm pix} = \sum n V$ and
interior column $= n_3 \ell$ exact, line-center $\tau_0$ (each component
and the blended $J'=1,2$ doublet) $=$ hand integral to $<2\times10^{-16}$;
thin case silent, $\times10^5$ scaled case fires the $\tau_0>1$ warning;
high-density run: column conservation $9\times10^{-15}$,
doublet/$J'{=}0$ $\tau_0$ ratio 8.003 $=$ the $f$-ratio 8.002 (that run is
genuinely optically thick at line center, $\tau_0 \sim 10^4$ — flagged as
designed).

- **H1 — Porter interpolation fidelity**: the existing 10830 interpolation
  reproduces the source grid across the Porter $(T_e, n_e)$ range (checks
  the current emissivity path; no new physics).
- **H2 — analytic limits**: the module reproduces (i) the low-density
  radiative limit $n_3/n_{\rm He^+} = \alpha_3 n_e/A_{31}$
  ($1.65\times10^{-7}$ at $n_e=100$, $T_e=10^4$ K) and (ii) the
  high-density plateau $\alpha_3/(q_{31a}+q_{31b})$, with the crossover at
  $n_{\rm crit} = A_{31}/(q_{31a}+q_{31b})$; each sink term
  ($\Phi_3$, e-impact, $n_{\rm HI} Q_{31}$) exercised on/off against a
  one-zone Python reference using the same coefficients.
- **H2b — external anchor**: $n_3/n_{\rm He^+}$ against Cloudy's He I
  model-atom level populations for a standard H II slab (Cloudy solves the
  full atom including 2$^3$S).  This is the FIRST external check of the
  coefficient set: H1, H2 and H5 all passed as internal self-consistency
  (the module reproduces the analytic limits and the one-zone Python form
  built from the same coefficients to machine precision), so what remains
  untested is whether the coefficients themselves are right as a group.
  * **What it actually tests.** The largest approximation in the closed form
    is $\alpha_3 = 0.75\,\alpha_B({\rm He\,II})$, i.e. that every triplet
    recombination cascades down to 2$^3$S (spin statistics 3/4).  Cloudy's
    100-level He I atom resolves that cascade, so part of its triplet
    recombination reaches 2$^3$P and bypasses 2$^3$S.  H2b also exposes
    $q_{13}$ (Bray 2000), the two e-impact sinks, and the Penning term
    $n_{\rm HI} Q_{31}$ as a group.
  * **Version, corrected 2026-07-26.** Use **c25.00**, not c23.01: the
    built binary on this machine is `~/CLOUDY/c25.00/source/cloudy.exe`,
    while `~/CLOUDY/c23.01/` holds data only (no executable), and the paper
    benchmarks already use c25.00, so this keeps one reference version.
    c23.01 remains the provenance of `he1_case_b.dat` (Section 1) and that
    does not change.
  * **How to run it.** Append output commands to the existing decks
    `tests/cloudy_c25/hii{40,20}_c25.in` and run the MoCHII side of the same
    two models with `par%hei_metastable` on.  The two codes already share
    density, composition, and source there, which is why these decks are the
    cheap route.  The c25.00 command is `save species densities` (renamed
    from `save species populations`; `source/parse_save.cpp:2330`).
  * **Criterion** (not set in rev 2; proposed here).  Tens of percent is the
    size of the closed-form approximation and should be reported as such; a
    factor-level disagreement is a coefficient error.  A disagreement can be
    attributed by its radial shape, because the two source terms dominate in
    different zones: $\alpha_3 n_e n_{\rm He^+}$ inside the ionized region
    (where $n_{\rm HeI}$ is small) and $q_{13} n_e n_{\rm HeI}$ beyond the
    front.
- **H3 — Porter-equivalent limit**: under case-B conditions (dense,
  optically thin, $\Phi_3$ negligible) an explicit total 10830 built from
  $n_3$ (collisional part) plus a collisionless recombination part matches
  the Porter table.  Compare the explicit total to the Porter total — never
  "Porter + explicit collisional" to anything (the rev-0 double count).
  This gate is OPTIONAL (only if explicit emission is implemented).
- **H4 — advection error, EXHALE reference**: run the local-equilibrium
  module on an EXHALE hot-Jupiter profile ($n$, $T$, $n_e$, $x_{\rm HeII}$
  from the EXHALE output) and compare $n_3(r)$ and the metastable column
  against EXHALE's own advected triplet solution.  This quantifies the
  local-equilibrium error directly against the in-house code rather than
  gating on OH18 paper figures; MoCHII must NOT claim the OH18 benchmark
  without an advection term.
  * **Why the exoplanet profile is the stress case, and an H II region is
    not.**  The 2$^3$S radiative lifetime is $1/A_{31} \simeq 7.9\times10^3$ s
    (~2.2 h).  In an H II region every flow time is far longer, so local
    equilibrium is safe (this is why the H II gates say nothing about the
    assumption).  In an escaping planetary atmosphere, $v \sim 10$ km s$^{-1}$
    over $r \sim 10^{10}$ cm gives a flow time of the same order as the 2$^3$S
    lifetime, so the advection term enters the balance at leading order.
  * **What exists on the EXHALE side** (checked 2026-07-26,
    `~/RT_Codes/ExoAtmosphere/EXHALE`): the merged He-triplet + metals system
    `System_HeH_TR_metals` in `src/modules/radiation/ionization_equilibrium.f90`
    (the triplet occupies `sys_x(4)`; the triplet kinetics keep the MINPACK
    solve, no analytic Jacobian), and a finished HD189733b run whose
    `output/` already holds `Ion_species.txt` and the advected
    `Ion_species_adv.txt` (plus `Hydro_ioniz{,_adv}.txt` for $n$, $v$, $T$).
    Nothing needs to be re-run on that side.
  * **Remaining work**: identify which column of `Ion_species_adv.txt` is
    $n(2^3{\rm S})$ (read the writer; the files carry no header), and write a
    thin driver that feeds the 1D radial profile to the local-equilibrium
    form.  MoCHII itself is not touched: this is a quantification gate, not a
    fix, because `src/read_generic_amr.f90` discards velocity input and the
    code cannot advect today.
- **H5 — absorption**: analytic column-to-$\tau$ check for a uniform slab
  (Stage 2 formula vs hand integral); the $\tau > 1$ warning fires where
  expected.

## 7. Integration and parameters

- `par%hei_metastable` (default off): compute the 2$^3$S diagnostic on the
  converged state.
- Output: $n_3$ (`n_2s3` block) and its column map as primary products; a
  `tau10833` map from Stage 2.  The Porter `HeI10830` emission row is
  UNCHANGED and nothing is summed onto it.
- Depends on the converged state like the other diagnostics; $\Phi_3$ is
  always integrated over the full band (the EUV bins contribute), and
  `add_fuv` with `efuv_min = 4.78` adds the soft-UV part (a rank-0 note
  fires when the band floor is above 4.78 eV).
- New data file: `data/atomic/hei_metastable.txt` (the Section 4 table:
  rate fits + the two VFKY96 $\sigma_3$ wings), generated by a
  `tools/fitting/make_hei_metastable.py` that transcribes the EXHALE
  expressions and cross-checks them against p-winds and CHIANTI at a grid
  of $(T, n_e)$ — the golden-table pattern the CX audit recommended.

## 8. Scope, order, and caveats

- **Order**: Stage 1 with $\Phi_3 = 0$ + H2/H2b first (self-contained,
  gives $n_3$ in the dense H II regime); then $\Phi_3$ from the band
  (H2's $\Phi_3$ limb + the exoplanet/PDR regime); Stage 2 absorption (H5)
  rides along; H4 when an EXHALE profile is prepared.
- **Order of the two open gates (2026-07-26)**: H2b first.  It is
  self-contained, the c25.00 binary and the two benchmark decks are already
  in place, and a failure points straight at the coefficient set.  H4 needs
  the EXHALE column identified and a 1D driver written, and its product is a
  number for the local-equilibrium caveat rather than a pass/fail.  H3 stays
  optional; if it is ever done, the explicit 10833 must REPLACE the Porter
  row, never be added to it (the Porter table is a collisional-radiative
  total).
- **Caveats to document at the code site**: (i) local equilibrium — the
  advection term is absent and H4 quantifies the error; (ii) 10830 can be
  optically thick in dense metastable regions (H5 warning; line transfer
  out of scope); (iii) the $n_{\rm HI} Q_{31}$ channel matters in the
  neutral/PDR zone — include it from the start (the PDR C II / O I
  H-impact precedent in `species_mod`); (iv) `hei_diffuse` currently
  assumes every triplet-channel recombination emits the 19.82 eV photon
  ($A_{31}$ route); at $n_e > n_{\rm crit}$ the collisional transfer to the
  singlets diverts that branch — a later consistency refinement between
  `hei_diffuse` and this module, negligible for the current gates.
- **Relation to the existing line**: the Porter case-B `HeI10830`
  emissivity is a collisional-radiative total and stays exactly as is.
  The metastable module adds a NEW diagnostic ($n_3$, column,
  $\tau_{10830}$); it does not add or modify any emission row.
