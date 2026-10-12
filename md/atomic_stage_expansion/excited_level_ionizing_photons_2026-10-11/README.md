# Metal ionizing photons from excited levels: survey, measured size, Part B design

Date: 2026-10-11 (KST). Plan row 4.6 of
[ATOMIC_PHYSICS_RESTORATION_PLAN_2026-10-11.md](../../ATOMIC_PHYSICS_RESTORATION_PLAN_2026-10-11.md).
Scope confirmed by the author: **Part A only (survey, measurement, design)**.
MOCASSIN v2 computes none of this physics. Its `emission_mod.f90::fb_ff` is a
ground-state, outer-shell Milne continuum, and `photon_mod.f90:925-940` lets
every line photon escape. So this row lies outside the v2 scope, and nothing in
`src/` or `data/` was changed.

The missing photons are those at or above 13.598 eV from metal ions, in three
groups:

1. **Free-bound continua of excited levels.** Radiative recombination into
   excited levels whose binding energy (IP - E_level) plus the electron energy
   exceeds 13.598 eV. Today, channel 6 carries only captures into the
   **ground term**.
2. **Recombination lines.** The radiative cascade lines that follow recombination
   (RR and DR), and the dielectronic-stabilization photons, at 13.598 eV and above.
3. **Collisionally excited lines** at 13.598 eV and above (for example
   O III 2s2p3 -> 2p2, O IV 2s2p2 -> 2s2 2p, N IV, Ne V, S V).

Each number below states how it was obtained. The measurement code and its
outputs are in [`measurement/`](measurement/) (see "Reproduction").

---

## 1. Summary

### Measured size (leaf-integrated, erg/s)

Abundances are the HII20 set (C, N, O, Ne, S). The excited-level sum is
(1) + (2) + (3); "ground term" is the Badnell partial rate into the ground term,
the physics of channel 6.

| state | (1) fb excited | (2) RR cascade | (2) DR stabilization | (2) DR cascade | (3) collisional | excited-level sum | / total diffuse | / channel 6 | / ground term (Badnell) |
|---|---|---|---|---|---|---|---|---|---|
| HII40, 2 iterations (Phase 4 run) | 3.6e34 | 2.9e35 | 6.8e35 | 1.5e35 | 4.7e32 | **1.14e36** | **1.9e-3** | 2.5 | 1.8 |
| 1e5 K blackbody, 2 iterations | 7.3e35 | 1.8e36 | 2.9e36 | 2.0e36 | 1.4e34 | **7.46e36** | **6.3e-3** | 3.7 | 5.0 |
| HII40 converged (m11, 43-stage roster) | 8.8e35 | 2.5e36 | 2.9e36 | 1.1e36 | 7.6e32 | **7.29e36** | **1.0e-2** | -- | 3.4 |

- **Excited levels dominate the metal recombination photons.**
  - Their ionizing output is 1.8-5 times the ground-term continuum.
  - Counting the ground term too, 57-72% of all metal recombination radiation,
    n_e n_+ (alpha IP + beta), lies at or above 13.598 eV.
  - Channel 6 transports only the ground-term part.
- **Main emitters.**
  - Converged HII40: C2+ -> C II and O2+ -> O II.
  - 1e5 K blackbody: those two plus C3+, N3+, O3+, O4+.
  - The C2+ DR stabilization photons alone are 2.0e36 erg/s in the converged
    HII40 state.
- **The He-ionizing share (>= 24.587 eV)** of the excited-level photons is 0.12
  (HII40, 2 iterations), 0.40 (HII40 converged) and 0.63 (1e5 K blackbody).
- **(3), the collisionally excited lines, is negligible in H II regions.**
  - It is 1e-4 to 2e-3 of the recombination part, and at most 1.2e-5 of the
    diffuse field.
  - Nearly all of it lies at 13.6-24.6 eV.
  - The leading emitters are C II 904 A and O II 2s2p4 -> 2p3 at
    T ~ 7-10 kK, and O IV / N IV at 1e5 K.
  - At T >~ 3e4 K the EUV share of the line cooling grows to 15-50% per ion
    (O III 0.15, O IV 0.32, S V 0.51 at 10^4.5 K; table in
    `measurement/euv/*.dat`). This matters for hot gas, not for these states.
- **Sensitivity to the cascade approximation for high-n captures** (Section 4.3):
  without the direct high-n -> low-n photons the totals change by <= 0.3%.
- **For comparison, Cloudy's heuristic** of one photon at 0.75 IP plus one at
  0.25 IP per RR+DR recombination, on the same states, gives an ionizing part of:

  | state | Cloudy heuristic | this measurement |
  |---|---|---|
  | HII40, 2 iterations | 1.5e36 | 1.1e36 |
  | 1e5 K blackbody | 8.9e36 | 7.5e36 |
  | HII40 converged | 9.3e36 | 7.3e36 |

  Cloudy emits this heuristic **in addition to** its ground-state Milne
  continuum (Section 2.4).

**Condition of the states.**
- The two Phase 4 runs are not converged. After 2 iterations, HII40 is still
  97% O+ with a median T_e of 6.5 kK.
- The converged HII40 state is the m11 run on the 43-stage roster. Its O III is
  the top stage, so O3+ recombination is absent there, and the run had no
  channel 6.
- The magnitudes are therefore indicative, not final. All three states agree
  on the order of magnitude, 2e-3 to 1e-2 of the diffuse field.

### Data available (66 recombining ions X^q+ -> X^(q-1)+, q = 1-6, 11 elements)

| source | ions covered | content |
|---|---|---|
| Badnell partial RR, ADAS adf48 (IC) | **53 / 66** | final-state-resolved rates. IC-resolved low levels for H- to Al-like targets; nl to n = 10, and to n = 100 for the jk12/eab12 Si-/P-like files; n-shells to 999; 19 temperatures z^2 (10-10^7) K |
| Badnell partial DR, ADAS adf09 (IC) | **52 / 66** | final-state-resolved rates per core transition (Delta n = 0, 1). Same layout. Initial ground and metastable parents |
| Zhang et al. 2026 archive | Fe3+-Fe6+ DR (adf09, level-resolved, per the existing source decision `fe_recombination_source_decision.md`) | **Not read here** |
| CHIANTI 11.0.2 `.rrlvl` / `.reclvl` | 6 / 66 (C IV, C V, N V, N VI, O VI from Badnell 2006; C VI, H-like, from Mewe et al. 1982) | RR into resolved levels. CHIANTI 11.0.2 has **no `.drlvl`**. `.auto` files (two-ion DR) exist for few ions. `.fblvl` lists the levels for CHIANTI's free-bound continuum (Verner ground + Karzas-Latter excited) |
| Cloudy c25.00, MOCASSIN, MAPPINGS V, CMacIonize, COLT | -- | Section 2.4. None transports level-resolved metal recombination photons. Cloudy uses a 0.75/0.25 IP pair; MAPPINGS uses Mewe-1986 hydrogenic excited continua |

**No final-state data in any source read:**
- Fe+ -> Fe I, Fe2+ -> Fe II and Fe3+ -> Fe III. Fe4+-Fe6+ DR may be in the
  Zhang 2026 archive.
- Ca+ to Ca4+.
- Ar+ and Ar2+.
- Cl+.

Badnell's partial RR stops at P-like targets and partial DR at P-like. The
missing ions are K- to S-like, Cl-like and Mn- to Ca-like targets.

### Part B design in one paragraph

- **New diffuse channel 8, metal excited-level recombination.** The rates come
  from tables, generated offline per recombining ion, of:
  - the ionizing free-bound level groups: threshold, rate share, and a Milne
    (Kramers) shape normalized to the Badnell partial rate;
  - the line yields per recombination at discrete energies: radiative cascade
    through the MoCHII CHIANTI atom, and DR stabilization.

  Both are on a temperature grid and normalized so that each leaf emits exactly
  n_e n_+ (alpha IP + beta) - (channel 6) in total, split between the
  transported band and the escaping remainder.
- **Energies are sampled continuously:**
  - for a line, its exact energy;
  - for the free-bound continuum, the threshold plus an exponential electron
    energy;
  - for DR stabilization, the resonance energy plus the binding energy.
- **Channel 9, collisionally excited EUV lines.** It reuses the
  `nlevel_cooling_mod` (T, n_e) solve with one extra coefficient table plus a
  coarse table of line fractions, about 15 MB in all. It is cheap but small; whether it
  is on by default is an open decision.
- **Cooling stays as it is:**
  - Recombination photons draw their IP part from the ionization budget and
    their kinetic part from the existing recombination cooling beta.
  - Collisional EUV lines are already in the line cooling, and absorbed photons
    return their energy through the photoionization heating tallies.
- **No recombination feeding of the n-level solve is needed for (2):** the EUV
  cascade is radiative at nebular density (A >~ 1e6 s^-1 >> n_e q).
- **Cost:** offline generation takes about 30 s per ion (Python). Runtime
  memory is under 20 MB for channel 8 and about 15 MB for channel 9. Runtime
  work is one extra coefficient per leaf and channel, plus sampling.

Details are in Section 5, gates in Section 6, open decisions in Section 8.

---

## 2. Data survey

### 2.1 What MoCHII carries now

- **Totals (all 66 transitions).** The element files `data/atomic/element_*.txt`
  carry the total RR (Badnell 2006 fits, ApJS 167, 334) and DR (Badnell et al.
  2003 project, A&A 406, 1151) rates in `RR`/`DR` rows, plus Zhang et al. 2026
  for Fe. `species_mod::alpha_rec` uses them.
- **Channel 6.** Milne on the outer-shell VFKY96/VY95 cross section with
  ground-term weights (`milne_recomb_spectrum_mod::milne_metal_*`). It covers
  ground-term captures only.
- **Excitation atoms.** `data/atomic/excitation/chianti/nlevel_*.txt` holds
  70 atoms, 6629 levels and 129,842 radiative rows, of which 73,512 are at or
  above 13.598 eV (counted from the files).
  - Most are full CHIANTI 11.0.2 bound-level atoms.
  - Some are small: ar_2 has 2 levels; cl_2, cl_3, cl_4 have 5; ca_4, ca_5,
    ca_6 have 2-5; o_1 has 7; s_1 has 5.
  - ne_3 and ne_4 contain n = 2 configurations only.

### 2.2 Badnell final-state-resolved partials (ADAS adf48 RR / adf09 DR)

Source: `https://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/` and `.../DR/`. The
index pages list the sequence directories.
- **RR** runs from Ur-like (bare) to P-like targets.
- **DR** runs from H-like to P-like targets, with no Ur-like.
- Every file holds the partial rates from the initial **ground and metastable**
  parents (PRTI) to:
  - IC-resolved final levels of the recombined ion (H- to Mg/Al-like targets);
  - `ILREP` nl-shells and `INREP` representative n-shells up to n = 999, per
    final parent (PRTF).
- **Si- and P-like files** (jk12, eab12) give **configuration-averaged
  parents and no IC level list**, only nl up to n = 100.
- **Ion coverage**, from the directory listings fetched 2026-10-11 (key:
  recombining ion; RR file prefix | DR core-transition codes):

| element | q = 1 | q = 2 | q = 3 | q = 4 | q = 5 | q = 6 |
|---|---|---|---|---|---|---|
| C | B-like nrb05 \| 22,23 | Be nrb05 \| 12,22,23 | Li nrb05 \| 12,22,23 | He nrb05 \| 12,13,22,23 | H nrb05 \| 12 | bare nrb05 \| -- |
| N | C nrb05 \| 22,23 | B \| 22,23 | Be \| 12,22,23 | Li \| 12,22,23 | He \| 12,22,23 | H \| 12 |
| O | N \| 22,23 | C \| 22,23 | B \| 22,23 | Be \| 12,22,23 | Li \| 12,22,23 | He \| 12,22,23 |
| Ne | F \| 22,23 | O \| 22,23 | N \| 22,23 | C \| 22,23 | B \| 22,23 | Be \| 12,22,23 |
| Mg | Na \| m23,m334 | Ne \| 23 | F \| 22,23 | O \| 22,23 | N \| 22,23 | C \| 22,23 |
| Si | Al saan05 \| 33,34 | Mg \| 23,334 | Na \| m23,m33,m34 | Ne \| 23 | F \| 22,23 | O \| 22,23 |
| S | P eab12 \| 33,34 | Si jk12 \| 33,34 | Al saan05 \| 33,34 | Mg \| 23,334 | Na \| m23,m33,m34 | Ne \| 23 |
| Cl | **none** (S-like) | P \| 33,34 | Si \| 33,34 | Al \| 33,34 | Mg \| 23,334 | Na \| m23,m33,m34 |
| Ar | **none** (Cl-like) | **none** (S-like) | P \| 33,34 | Si \| 33,34 | Al \| 33,34 | Mg \| 23,334 |
| Ca | **none** (K-like) | **none** (Ar) | **none** (Cl) | **none** (S) | P \| 33,34 | Si \| 33,34 |
| Fe | **none** (Mn-like) | **none** (Cr) | **none** (V) | **none** (Ti) | **none** (Sc) | **none** (Ca) |

**Totals.** RR partials: 53 ions. DR partials: 52 ions (no DR for bare C6+).

**Files with alternative calculations.** B-like C+, C-like N+ and N-like O+
have extra `*ic222` / `*ic223` files, and C+ has both nrb00 and za00 `ic22`.
The measurement used one file per core-transition code (za00 preferred, the
222/223 alternatives skipped). Which is the recommended set is an open question
for Badnell.

**Checks made on the files read** (`measurement/checktot.py`; files were read
with PRTI = 1 only):
- **n-sum.** The summed n-resolved partials, interpolated across the
  representative n in log-log, reproduce each file's own total table to 0.5%
  over 1e3-1e5 K:
  - RR, C+ / O3+ / S+ / C6+: 0.995-1.006.
  - DR, S+: 0.995-1.001.
  - DR, O3+ `ic22`: 1.003 at 1e4 K, but 1.12-1.25 above 4.5e4 K. That excess
    is from the high-n interpolation, outside the nebular range.
- **Against the MoCHII fits.** The partial totals agree with the MoCHII Badnell
  fits at 1e4 K:
  - RR: 0.96-1.04 for all 30 ions measured.
  - DR: 0.83-1.28 for most ions.
- **Exceptions in DR:**
  - N+ (0.15) and O+ (0.07). The DR fits for these include low-temperature
    resonance contributions that the partial files do not resolve the same way.
    In the measurement the partial shape is scaled to the fit rate.
  - Targets whose DR is negligible at 1e4 K (He-like, H-like) give ~0, as
    expected.

### 2.3 CHIANTI 11.0.2 (`~/RT_Codes/CHIANTI/dbase`)

- **Level-resolved recombination.**
  - `.rrlvl` (Badnell 2006 RR partials): 34 files in the whole database. Only
    c_4, c_5, n_5, n_6, o_6 and o_7 are recombination products in the MoCHII
    roster; ne_8 and higher stages are outside it.
  - `.reclvl` (Mewe et al. 1982): c_6, n_7 and o_8 (H-like).
  - **No `.drlvl` file exists** in 11.0.2.
  - The `.auto` autoionization files (two-ion DR model) exist for c_4, c_5,
    n_5, n_6, o_6, o_7, mg_2, si_4, s_6, ar_8 and a few others.
  - So CHIANTI is not a source for this row beyond what Badnell already
    provides.
- **`.fblvl` and the free-bound continuum.** CHIANTI's free-bound continuum
  (ChiantiPy `Continuum.freeBound`) uses Verner & Yakovlev for the ground level
  and Karzas & Latter (1961) Gaunt factors for the excited levels listed in
  `.fblvl` (statistical weights, n, l, energies). It is a hydrogenic
  excited-level Milne continuum with no level-resolved rate data.
- **Header error.** The `c_5.rrlvl` header says "recombination onto C V to form
  C IV". The file is C VI -> C V. This is external data; noted only.

### 2.4 Other codes (read-only survey of the source; file:line)

**Cloudy c25.00**
- **(1) Free-bound.**
  - H- and He-like iso-sequences: level-resolved RRC for every element; the
    destroyed part goes to OTS (`rt_ots.cpp:283-288`), the escaping part
    outward (`rt_diffuse.cpp:539,586`).
  - Other ions: valence-shell **ground** Milne only, on the stored opacity
    (`rt_diffuse.cpp:247-261`), transported outward.
- **(2) Recombination lines.**
  - Every non-iso recombination (RR + DR) emits one photon at **0.75 IP** and
    one at **0.25 IP** (`ion_recomb.cpp:91-95`, checked here;
    `cont_createpointers.cpp:1062-1066`, checked here). They go to OTS and
    outward (`rt_ots.cpp:361-406`).
  - Ground captures thus emit both the Milne continuum and the pair, so the IP
    is counted twice for them. This is inferred from the structure, not
    computed.
  - DR satellite lines exist only for the iso-sequences
    (`rt_diffuse.cpp:395-406`).
- **(3) Collisional lines.** Every level2 and Stout/CHIANTI line at any energy:
  the locally destroyed part `Pop*A*Pdest` goes to OTS (`rt_ots.cpp:410-437`),
  the escaping part outward. The default atoms are truncated (15 levels, Fe 25);
  for example the O III 2p3s levels (374 A) are absent by default.
- **Energy.** Heavy-element recombination cooling is `kT*alpha_RR` only
  (`cool_eval.cpp:486-514`).

**MOCASSIN 2.02.73.2**
- **(1)** Ground Milne only, Z = 3-19 (`emission_mod.f90:472-530`).
- **(2)** Metal recombination lines are output only.
- **(3)** Collisionally excited lines are computed, but every line packet
  escapes (`photon_mod.f90:925-927`). The line cooling sum runs over levels
  1-10 (`update_mod.f90:1103-1106`) while the emission runs over 17 levels.
  This is an external-code observation, not run.

**MAPPINGS V 5.2.1**
- **(1)** Ground Milne. Excited levels are tabulated for the H- and He-like
  ions; for the other ions it uses the Mewe 1986 hydrogenic n > n0 sum
  (`freebound.f:201-240`), transported by two-stream.
- **(2)** Heavy-element H-like series only. No DR stabilization lines.
- **(3)** XR3/XRL collisional cascade atoms (CHIANTI 8) are transported and
  enter J_nu (`totphot.f`). The line cooling equals the line emission
  (`multilevel.f:432-436`, `reson.f:567-570`).

**CMacIonize, COLT**
- Neither has any of (1)-(3) for metals.
- COLT's diffuse sources are H/He free-free, free-bound and two-photon only
  (`ionization/ray_trace.cc:446-560`).

The Cloudy, MOCASSIN, MAPPINGS, CMacIonize and COLT entries were read from the
source by a survey agent. The two Cloudy lines marked "checked here" were
re-read directly. No code was run.

---

## 3. Finding: channel 6 against Badnell's ground-term partial

The measurement gives a second determination of the ground-term RR rate, the
physics of channel 6: the Badnell partial into the ground term, IC-resolved
where available. Channel 6 instead uses Milne on the VFKY96 outer-shell fit.

Ratio Milne(VFKY96) / Badnell partial at 1e4 K (`measurement/milne_check.log`):

| range | ions |
|---|---|
| 0.93-1.06 | C III-VI, N II, N IV-VI, O IV-VI, Ne I, Ne IV, Ne VI, S VI |
| 0.84-1.21 | C I (1.21), C II (1.08), N I (0.88), N III (0.84), O I (0.91), O II (0.88), O III (1.11), Ne II (1.17), Ne V (1.11), S IV (1.14), S V (0.86) |
| **outliers** | **Ne III 1.92**, **S I 20.6**, **S II 3.9**, **S III 0.16** |

- **H-like check.** C VI gives 1.001, so the Milne normalization is right. The
  measurement script first used g_n/(2 g_+) with a prefactor smaller by 2,
  found the factor of 2 on H-like ions, and corrected it. MoCHII's
  `milne_prefactor` is correct.
- **Why the outliers deviate** (resolved 2026-10-11, `measurement/ground_capture/`
  and the Phase 4 entry of `../../ATOMIC_PHYSICS_RESTORATION_RESULTS_2026-10-11.md`):
  - The S I and S II Badnell files are configuration-resolved (no IC level
    list); the S III file is IC-resolved. The statistical-weight split of
    S I/S II ignored fractional parentage (shares 1 and 1/3, not 9/15 and
    4/20), and the comparison left out DR, which VFKY96 (OP smoothed over
    resonances) contains: with RR + DR into the ground term S I is 0.73.
  - S II (2.29) and S III (0.17) are disagreements between the R-matrix
    calculations (OP, Nahar's unified ground-state rates, which channel 6
    reproduces to 2% and 26%) and AUTOSTRUCTURE; not adjudicated.
  - Ne III: 1.92 against RR alone, 0.60 against RR + DR. VFKY96 misses
    near-threshold resonances present in Nahar (2019)
    and in the ALS measurement; channel 6 was 0.56x the resonant ground capture
    at 1e4 K.  Since 2026-10-11 (later the same day) Ne III takes Nahar's table
    for transport and channel 6 (`PHTABLE`; results md, Phase 4).
- **Effect on totals.** Summed over the run, channel 6 printed 4.50e35
  (HII40, 2 iterations) and 2.00e36 (1e5 K). The Badnell ground-term partials
  on the same leaves give 6.28e35 and 1.48e36 (survey method: RR only, S I/S II
  shares by statistical weight; not recomputed).
- **Decision needed.** The ground-term rate and spectrum should come from one
  source for channels 6 and 8 together (Section 8, D1).

---

## 4. Measurement method

### 4.1 States

All states use n_H = 100 cm^-3 from the input FITS, matched to the leaves by
position.

| state | T_e and n_e | stage fractions | diffuse total |
|---|---|---|---|
| `runs/work/hii40` (Phase 4 HII40, 2 iterations) | `rates.h5` | `rates.h5` | final-pass ledger of `run.log` |
| `runs/work/bb1e5` (1e5 K blackbody, 2 iterations) | `rates.h5` | `rates.h5` | final-pass ledger of `run.log` |
| `results/continuous_energy/m11_hecx_solver/hii40_converged_rates.h5` | same file | top stages C IV, N III, O III, Ne III, S IV | 7.2451e38 from its log |

Volumes are LeafSize^3.

### 4.2 Collisionally excited lines (3)

- **Program.** `measurement/euv_emissivity.f90` is built against the Phase 4
  work tree's `define.f90` and `nlevel_mod.f90`, i.e. the code path the runs
  used.
- **Solve.** `nlevel_cooling_prepT` and `nlevel_populations` solve each data
  snapshot atom (`data_snapshot/atomic/excitation/chianti/nlevel_*.txt`) at
  log T = 3.4-5.0 (0.025 dex) and n_e = 1-1e4 cm^-3.
- **Output.** Sum of pop_u A h nu / n_e over the radiative rows with
  h nu >= 13.598 eV, with the emitting level chosen by energy as in
  `emitting_level`.
- **Leaf integration.** Bilinear interpolation in (log T, log n_e) times
  n_e n_ion V.
- **Line count.** The number of lines needed to carry 1 - 1e-4 of the EUV
  emissivity at every node of a T x n_e grid (16 x 5) is 7029 in all, summed
  over the 70 atoms (`measurement/lines_needed.txt`, per atom).

### 4.3 Recombination photons (1), (2)

Per recombining ion X^q+ -> X^(q-1)+ (`measurement/recomb_photons.py`):

- **Rates.** Badnell adf48 (RR) and adf09 (DR) partials, initial ground parent
  PRTI = 1, interpolated in log-log onto log T = 3.4-4.8.
  - IC-resolved final levels are used where the file has them.
  - The remainder of each nl-shell comes from `ILREP`.
  - n > 10 (n > 100 for Si/P-like) comes from `INREP`, interpolated across the
    representative n, with the l distribution of the highest nl-resolved n.
  - Each process is then scaled to the MoCHII total rate (element-file
    `RR`/`DR` rows).
- **Final-state energies and cascades** come from the MoCHII CHIANTI atom of
  the product ion.
  - Configurations are matched after filling implicit closed shells to the
    ion's electron count and removing closed shells. CHIANTI spells them
    `3s`, `2s22p5`, `1s2.2s` and `3p2(3P)3d` for the same kind of entity.
  - An IC level is matched by configuration, (2S+1, L, J). Otherwise the
    configuration group is weighted by g.
  - The cascade ionizing energy of each atom level is
    C(i) = sum_j B_ij [h nu_ij (h nu_ij >= 13.598 eV) + C(j)], with B_ij the
    radiative branching ratios (A values only).
- **Captures to states outside the atom** (23-83% of the RR rate at 1e4 K,
  mostly high n) descend to the highest atom configuration with the same parent and
  l' <= l, and take its C.
  - A direct jump photon is added when the remaining energy exceeds 13.598 eV.
  - Removing that jump photon changes the totals by <= 0.3%
    (`leaf_nodirect.json`).
  - Core-excited parents that the atom lacks add one core photon E_parent.
- **First photon.**
  - RR (free-bound): Kramers shape normalized to the partial rate. The
    energy-weighted distribution in the electron energy E is exp(-E/kT), so the
    mean photon energy is kT / [e^b E1(b)] with b = B/kT, and the fraction at
    or above 13.598 eV is exp(-max(0, 13.598 - B)/kT).
  - DR (stabilization): h nu = E_res + B_f, with E_res weighted by the terms of
    the MoCHII DR fit, c_i exp(-E_i/T).
- **Ground term.** Ground-term captures are separated from the excited ones.
  The ground-term part is the physics of channel 6; the excited-level
  free-bound part is (1).
- **Energy check.** The total radiated energy per leaf (all photons, ionizing
  or not) reproduces n_e n_+ alpha (IP + <E_in>) and agrees with the bound
  n_e n_+ alpha (IP + kT) to 1-2% (columns `recEtot` / `bound` in
  `leaf_luminosities.log`).

### 4.4 Results by ion

Ions carrying >= 1% of the excited-level sum, in erg/s. Full tables are in
`measurement/leaf_luminosities.log`; they are produced by `tables.py`.

**Converged HII40 (m11):**

| recombination | (1) fb excited | (2) RR cascade | (2) DR stab. | (2) DR cascade | ground term | excited sum |
|---|---|---|---|---|---|---|
| C2+ -> C II | 9.1e31 | 7.40e35 | 2.02e36 | 3.18e35 | 1.13e36 | 3.07e36 |
| C3+ -> C III | 1.2e34 | 1.8e34 | 5.5e34 | 5.6e34 | 1.1e33 | 1.40e35 |
| N2+ -> N II | 7.1e34 | 1.30e35 | 5.2e34 | 2.2e34 | 1.44e35 | 2.74e35 |
| O2+ -> O II | 7.82e35 | 1.37e36 | 6.05e35 | 5.84e35 | 6.19e35 | 3.34e36 |
| Ne2+ -> Ne II | 8.0e32 | 1.68e35 | 1.8e34 | 4.9e34 | 1.02e35 | 2.35e35 |
| S2+ -> S II | 5.7e33 | 1.7e34 | 9.4e34 | 1.4e34 | 1.4e33 | 1.30e35 |

**1e5 K blackbody (2 iterations):** the leading recombinations are:

| recombination | excited sum |
|---|---|
| O2+ -> O II | 1.73e36 |
| C3+ -> C III | 1.68e36 |
| C2+ -> C II | 1.31e36 |
| O3+ -> O III | 6.9e35 |
| O4+ -> O IV | 6.2e35 |
| N3+ -> N III | 5.3e35 |

The largest (3) emitters are O IV (9.1e33), N IV (2.2e33) and O III (1.3e33).

### 4.5 Approximations of the measurement

- **Free-bound shape.** Kramers rather than level-specific photoionization
  cross sections. It integrates exactly to the partial rate and changes only
  where within each group's continuum the energy falls.
- **DR resonance energies** come from the total-fit terms, not per final state.
- **No configuration-averaged splitting for Si/P-like targets** (S+, S2+ and
  Cl, Ar, Ca analogues): statistical weight is used.
- **Cascades** are purely radiative, with no collisional redistribution. This
  is valid for the EUV-emitting levels at n_e ~ 100.
- **Small atoms** (o_1, s_1, ne_3, ne_4) make the descent rule do more of the
  work for their products.
- **Line transfer.** No resonance-line transfer and no on-the-spot fractions:
  the numbers are emission, not absorption.

---

## 5. Part B design proposal

### 5.1 Representation

**Channel 8: metal excited-level recombination photons.** It is on with
`recombination_continuum = 'transported'`, as channels 6 and 7 are.

1. **Generator.** Offline, a new `tools/fitting/make_metal_recombination_photons.py`
   writes `data/atomic/recombination_photons_<el>_<stage>.txt` per recombining
   ion, with a provenance header (Badnell file names and checksums, CHIANTI
   atom checksum, MoCHII fit rows).
   - **Temperature grid:** log T = 3.0-5.0 at 0.05 dex (41 nodes), within the
     Badnell z^2-scaled range.
   - **Free-bound groups.** Each group k has a threshold B_k (from the CHIANTI
     level energy and the ETH ionization energy), a rate share
     a_k(T) / alpha_RR(T), and the ground-term flag. Groups are merged within
     0.05 eV, about 20-200 per ion.
   - **Lines.** Each line m has an exact energy and a yield y_m(T) in photons
     per recombination, kept separately for RR and DR. Lines below 13.598 eV
     are kept as one "escaping energy" scalar per node and are not listed.
   - **DR stabilization photons.** Per final-state group, they are stored as
     (B_f, share). Their energy, E_res + B_f, is sampled from the DR fit's
     E_i weights, which is continuous in E_res.
2. **Runtime emission per leaf:**
   n_e n_+ [alpha_RR (fb + RR-line terms) + alpha_DR (stab + DR-line terms)],
   with the MoCHII alpha. One coefficient per ion per T node gives the leaf
   luminosity; the stratified leaf choice is the one channels
   1-7 use.
3. **Sampling a photon** (continuous energy, no bins):
   - pick the ion by its leaf luminosity share;
   - pick fb / line / stab by its share, at the leaf T: choose the T node with
     the interpolation weight, an exact mixture sample;
   - **fb:** pick a group by its energy share. Then E has the photon-number
     density exp(-E/kT)/(B_k + E); sample E from the exponential and accept
     with probability B_k/(B_k + E), which is exact. Then h nu = B_k + E.
   - **line:** take the exact line energy.
   - **stab:** draw E_res from the fit weights, then h nu = E_res + B_f.
   - Only h nu >= 13.598 eV (or `eion_min`) is launched. The rest is booked as
     escaping energy in the ledger.
4. **Ground term.** Channel 6 keeps the ground term, and the channel 8 tables
   exclude it, so there is no double count. Decision D1 is whether channel 6
   should take its rate from Badnell too.

**Channel 9: collisionally excited EUV lines.** Its default is decision D3.
- `nlevel_cooling_mod` already solves every atom on a 401 x 161 (T, n_e)
  table. Two arrays are added:
  - eps_EUV(T, n_e), 0.5 MB per ion, 40 MB for 77 ions, node-shared;
  - the line CDFs on a coarse grid of about 30 T x 9 n_e nodes for the pruned
    line list (1 - 1e-4 coverage). That is 7029 lines in all, about 15 MB.
- **Sampling:** the leaf luminosity from eps_EUV, then the line from the CDF of
  the nearest coarse node chosen by interpolation weight; the exact line
  energy.

### 5.2 Energy bookkeeping

**Recombination.**
- Per leaf: channel 6 + channel 8 (transported) + escaping part =
  n_e n_+ (alpha IP + beta), to rounding, for every ion.
  - beta is the MoCHII recombination cooling, kT[(3/2) alpha + d alpha/dln T].
  - The kinetic part, beta, is carried only by the first photon (free-bound
    electron energy for RR, E_res for DR).
  - The generator rescales its Kramers kinetic energy and the E_res mean so
    that the mean kinetic energy per recombination equals beta/alpha at each
    node. This is the acceptance identity of plan row 4.6.
- **The IP part is not cooling.** It came from the photon that ionized the ion
  and was booked as ionization, not heat. Nothing changes in
  `thermal_mod`/`cooling_mod`.
- **Absorbed channel 8 photons** heat by h nu - IP_absorber through the
  existing photoionization heating tallies. This is the same as channels 1-7.

**Collisional lines.**
- They are already in the n-level line cooling, as the full sum of
  pop A h nu, and stay there.
- Transporting them moves the energy from "escapes" to "absorbed elsewhere".
  The absorber's heating tally returns it, so no cooling term changes.
- The gate: the channel 9 leaf luminosity equals the EUV part of the leaf line
  cooling, to rounding.

### 5.3 Interaction with the continuum modes

| `recombination_continuum` | channel 6 now | channels 8, 9 proposed |
|---|---|---|
| `'transported'` | transported | transported |
| `'escaping'` | escape | escape (ledger only) |
| `'on_the_spot'` | escape (documented in `diffuse_mod`) | escape, the same rule. Alternative (decision D2): absorb on the spot by the opacity shares of H0/He0/He+ at the photon energy, as `hei_recombination_photons_mod` does for the He I cascade |

### 5.4 Line transfer

- **The EUV resonance lines of the emitting ion are optically thick.**
  - Example: O III 2p3s -> 2p2 near 374 A. With f ~ 0.1, n(O2+) ~ 3e-2 cm^-3
    and a thermal width at 1e4 K, tau_0 is about 2e2 across one 0.15 pc leaf.
    This is an order-of-magnitude estimate, not computed from the atom.
  - So these photons scatter locally before the H/He continuum absorbs them.
- **Proposal: transport them as continuum photons**, without resonance
  scattering. Energy and band are exact; the absorption site is approximate.
  - The cascade lines into excited lower levels and the DR satellites (shifted
    by the spectator electron) are not trapped.
  - Resonance scattering is decision D4.
- **Gate (g) below** measures the line-centre optical depth per leaf. That
  bounds the error of this simplification.

### 5.5 Does (2) need recombination feeding in the n-level solve?

**No, not for the ionizing photons.**
- The levels that emit at or above 13.598 eV have A >~ 1e6-1e10 s^-1 (O III
  2p3s: 3.6e9 s^-1). Collisional rates at n_e = 100 are n_e q <~ 1e-5 s^-1.
- So the cascade is pure radiative branching and adds linearly to the
  collisional populations. The offline cascade tables are exact in that limit.

**Recombination feeding matters for something else.**
- It changes the metastable populations and so the forbidden lines (the known
  recombination contribution to [O II] 7320/7330 and [N II] 5755).
- That is an optical-line diagnostic, plan row "metal recombination lines"
  (`metal_recombination_lines`, off). It needs a source vector in
  `nlevel_populations`: the recombination flux into each level from the same
  tables, `R_i = n_e n_+ alpha_eff,i`.
- It is not needed for 4.6.

### 5.6 Ions without final-state data

Fe+ to Fe3+ (and Fe4+-Fe6+ unless the Zhang 2026 adf09 files hold final
states), Ca+ to Ca4+, Ar+, Ar2+ and Cl+. The options are decision D5:
- (a) the hydrogenic approximation on the CHIANTI atom: captures distributed
  by the Seaton hydrogenic partial alpha_n(z, T) over n, then the same descent
  and cascade;
- (b) requesting AUTOSTRUCTURE partials from Badnell, whose page offers
  "other elements may be added on request";
- (c) leaving them out with the bound printed.

For Fe I-III and Ca I-III, the product IP minus the lowest excitation is mostly
below 13.6 eV, so their ionizing yield is small. This can be checked with the
energy bound n_e n_+ alpha (IP + kT) on any state.

### 5.7 Cost

| item | cost |
|---|---|
| Offline generation | about 25 s per ion in the Python measurement code (30 ions in about 12 min); 66 ions about 30 min. Badnell files: 97 files (17 MB gzipped) for the 30 ions measured; all 53 covered ions roughly twice that |
| Data files | about 100 lines x 41 T nodes x 66 ions, about 3 MB of text |
| Runtime memory, channel 8 | under 20 MB (node-shared) |
| Runtime memory, channel 9 | eps_EUV 40 MB + line CDFs about 15 MB (node-shared) |
| Setup, channel 9 | one extra sum over the radiative rows at each of the 64,561 table nodes already solved. Measured: 325 nodes of the 177-level O III atom take 1.16 s including the LU solve, so the sum adds a negligible share |
| Transport | about 0.2-1% more diffuse luminosity. With the packet count proportional to luminosity, about 0.2-1% more packets |

---

## 6. Gates to write (Part B)

| gate | test | tolerance |
|---|---|---|
| (a) | Energy closure per leaf and per ion: ch6 + ch8 transported + escaping = n_e n_+ (alpha IP + beta) | 1e-10 relative |
| (b) | Kinetic part: the mean of (h nu1 - B) over sampled channel 8 first photons equals beta/alpha minus the channel 6 share | Monte Carlo 3 sigma at 1e6 samples |
| (c) | Partial-rate sums reproduce the adf48/adf09 total tables and the MoCHII fits at 1e4 K (Section 2.2) | 1% / documented ratios |
| (d) | Ground-term rate, channel 6 against the Badnell partial, per ion; outliers listed | report |
| (e) | Cascade branching: sum of yields x energy + escaping = E_level for every atom level | 1e-12 |
| (f) | H-like analogue: the C VI cascade from the generator against the SH91 hydrogenic cascade of `hydrogenic_recombination_mod`, scaled to z = 6 | 1% in the 2p, 3p yields |
| (g) | Line-centre optical depth of the top 10 EUV resonance lines per leaf on the converged HII40 state | report |
| (h) | Channel 9 luminosity = EUV part of the line cooling per leaf | 1e-12 |
| (i) | `escaping` / `on_the_spot` runs bit-identical to before for channels 1-7 | bitwise with pinned collectives |
| (j) | Measured share on HII40 and 1e5 K: converged before/after, L(Hbeta), He I/He II lines, <T>, fronts | report |

---

## 7. What was checked and what was not

**Checked.**
- Badnell partial sums against the file totals (5 files) and against the
  MoCHII fits (30 ions).
- H-like Milne normalization (C VI 1.001).
- Energy totals against the analytic bound (1-2%).
- Sensitivity to the direct-jump rule (<= 0.3%).
- The n-level EUV emissivity through the same Fortran module the runs used.

**Not checked.**
- The Zhang 2026 Fe adf09 archive.
- Mg, Si, Cl, Ar, Ca and Fe: their abundances are 0 in the measured runs.
- Resonance-line optical depths: only an order-of-magnitude estimate.
- A Badnell-recommended choice among the alternative `ic222`/`ic223` files.
- The DR resonance-energy distribution per final state.
- Converged states of the 77-stage code: none exists yet.

**Outside this task, found on the way (reported, not changed).**
- The project `CLAUDE.md` header says the default `metal_excitation_model` is
  `mocassin_v2`. Both trees' `src/define.f90` set `'chianti'`
  (`src/define.f90:1221` in the main tree), and the Phase 4 runs log the
  CHIANTI atoms.
- CHIANTI `c_5.rrlvl` header text (Section 2.3).

---

## 8. Open decisions

- **D1. One source for the ground-term capture.** Options:
  - (a) Keep channel 6 on Milne/VFKY96, consistent with the photoionization
    cross section by detailed balance.
  - (b) Take both channels 6 and 8 from the Badnell partials.

  They disagree by 1.9x (Ne III) and 0.16-21x (S I-III) at 1e4 K.
- **D2. On-the-spot mode.** Channels 8 and 9 escape, as channel 6 does today,
  or are deposited on the spot into H0/He0/He+ by opacity share.
- **D3. Channel 9 (collisional EUV lines) by default.**
  - It is 1e-4 to 2e-3 of channel 8 in H II regions.
  - It becomes relevant at T >~ 3e4 K, where 15-50% of an ion's line cooling
    is EUV.
  - Under rule 2 ("physically necessary means on"), this decides whether a
    small term is required.
- **D4. Resonance scattering.** Transport the EUV resonance lines as continuum
  photons (proposed), or treat local trapping with a destruction probability
  (Cloudy's `Pdest`).
- **D5. Ions without final-state data** (Section 5.6): hydrogenic
  approximation, a request to Badnell, or omission with a printed bound.
- **D6. Which Badnell file set to use** for C+, N+ and O+ (`ic22` vs
  `ic222`/`ic223`, nrb00 vs za00).

---

## 9. Reproduction

All files are in `measurement/`.

- **Badnell files.** The 97 Badnell files used are listed in
  `badnell_files_used.txt`. Download them from
  `https://amdpp.phys.strath.ac.uk/tamoc/DATA/{RR/RR_web/adf48,DR/DR_web/adf09}/<seq>like/`
  into `$BADNELL_DIR/{rr,dr}/`, with `%23` in the URL written as `#` in the
  file name.
- **`euv_emissivity.f90`**
  - Build it against copies of the Phase 4 work tree `define.f90` and
    `nlevel_mod.f90`: `ifx -O2 -cpp ... -qmkl=sequential`.
  - Run `euv.x nlevel_<el>_<st>.txt > euv/<el>_<st>.dat`.
  - Columns: log T, log n_e, eps_EUV, all-line cooling, then the
    [13.598, 24.587), [24.587, 54.418) and >= 54.418 eV parts, all in
    erg cm^3 s^-1 per n_e n_ion.
- **`euv_lines.f90`**: the line counts in `lines_needed.txt`.
- **`recomb_photons.py`** writes `recomb_coefficients.json`.
  - `DESCENT_NO_DIRECT=1` gives the sensitivity run.
  - `RECOMB_JSON` selects the output name.
- **`leaf_luminosities.py`** writes `leaf_luminosities.json` and `.log`.
  - It reads the three states (paths in `CASES`).
  - `EUV_DIR` points to the `.dat` tables.
- **`tables.py`** writes the markdown tables.
- **`milne_check.py`** gives Section 3.
- **`cloudy_heuristic.py`** gives the Cloudy 0.75/0.25 IP comparison.
- **`checktot.py`** gives the partial-sum checks.

References (verified in NASA ADS):
- Badnell 2006, ApJS 167, 334 (2006ApJS..167..334B).
- Badnell et al. 2003, A&A 406, 1151 (2003A&A...406.1151B).
- Karzas & Latter 1961, ApJS 6, 167 (1961ApJS....6..167K).
- Storey & Hummer 1991, CPC 66, 129 (1991CoPhC..66..129S).
