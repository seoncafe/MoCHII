# Igarashi (2023) Supplementary Data: What They Can and Cannot Tell

Date: 2026-09-25.

**Source.** `references/Igarashi_2023_psacc823supp1.pdf` holds the
supplementary tables of:
- Igarashi, A. 2023, *Physica Scripta* 98, 055402, "Cross sections in
  He⁺ + H collisions for the energy region 1–100 keV u⁻¹"
  (doi 10.1088/1402-4896/acc823, ADS 2023PhyS...98e5402I);
- the `acc823` in the file name matches the DOI;
- the method, from the ADS abstract: a semiclassical atomic-orbital
  close-coupling calculation with a two-electron treatment.

**What was not read.** The article itself was not obtained. Two things are
therefore not confirmed from the text:
- the energy convention, taken here as the lab energy per nucleon of the
  He⁺ projectile on an H target;
- which final states are summed in "capture". It is taken to be total
  electron capture into all He states.

## What the tables contain

Four sets of cross sections [cm²], at 24 energies from 1 to 100 keV/u:
- elastic, capture, excitation and ionization in He⁺ + H;
- excitation of H(2s), H(2p₀), H(2p₁), H(2p), H(n = 2);
- excitation of H(3s), H(3p), H(3d) and the Balmer-α emission.

## 1. No bearing on the thermal rate

**The energies.** 1 keV/u of He⁺ on H is a CM energy of 0.805 keV, and the
tables span 0.805–80.5 keV (CM).

**The Maxwell weight there.** The Maxwell weight exp(−E/kT) at 0.805 keV is
2.7×10⁻⁴¹ at 10⁵ K and underflows to zero at 10⁴ K.

**Consequence.** These data cannot confirm or refute the thermal
non-radiative rate from Loreau et al. (2014), which is ~1×10⁻¹⁶ cm³ s⁻¹ at
6×10³–5×10⁴ K (`LOREAU_2014_RATE_2026-09-25.md`). They also cannot settle
the conflict between Zygelman et al. (1989) and Loreau (2014) at 5–100 eV,
since that range lies below them.

## 2. Total capture versus the ground-state channel

Loreau's process (1), H⁺ + He(1s²) → H(1s) + He⁺(1s), gives the forward
**ground-state** capture He⁺(1s) + H(1s) → He(1s²) + H⁺ by reciprocity, at
the same CM energies:

  σ_f = ¼ (E_r/E_f) σ_r(E_r).

| E [keV/u] | E_cm [eV] | Igarashi capture (total) | Loreau ground channel, reciprocity | total / ground |
|---|---|---|---|---|
| 1 | 805 | 3.3×10⁻¹⁷ | 5.1×10⁻²⁰ | 650 |
| 2 | 1 610 | 4.9×10⁻¹⁷ | 6.5×10⁻¹⁹ | 75 |
| 5 | 4 025 | 1.0×10⁻¹⁶ | 7.3×10⁻¹⁸ | 14 |
| 10 | 8 051 | 2.1×10⁻¹⁶ | 2.4×10⁻¹⁷ | 9 |
| 20 | 16 102 | 2.6×10⁻¹⁶ | 4.2×10⁻¹⁷ | 6 |
| 50 | 40 254 | 7.8×10⁻¹⁷ | 1.9×10⁻¹⁷ | 4 |
| 100 | 80 508 | 1.3×10⁻¹⁷ | 4.0×10⁻¹⁸ | 3 |

**Meaning.**

- **Capture at keV energies is mostly not into the ground state.** Capture
  into He(1s²) is 1/650 of the total at 0.8 keV and 1/3 at 80 keV. The rest
  goes into excited He, He(1s nl). For those final states the triplet
  entrance channel, statistical weight 3/4, is open, whereas He(1s²) can be
  reached only from the singlet (1/4). Two consequences follow:
  - a **total** capture cross section, or a fit to one, is not the
    ground-state channel that a nebular He⁺ + H rate needs;
  - extrapolating total capture down to eV energies would overstate that
    channel by the same kind of factor.

- **The n = 2 capture channels open between 10⁴ and 10⁵ K.** They are
  endothermic, with defects from the NIST level energies:

  | Final state | Defect |
  |---|---|
  | He(1s2s ³S) + H⁺ | 8.8 eV |
  | He(1s2s ¹S) + H⁺ | 9.6 eV |
  | He(1s2p ³P) + H⁺ | 10.0 eV |
  | He(1s2p ¹P) + H⁺ | 10.2 eV |

  At 10⁴ K their Boltzmann factors are ≤ 4×10⁻⁵. At 5×10⁴–10⁵ K they are
  0.1–0.35. So in hot, partly neutral gas the He⁺ + H → He + H⁺ rate may be
  carried by excited-state capture. That capture is in neither Loreau's
  process (1) nor Zygelman's two-state (A ¹Σ⁺ / X ¹Σ⁺) calculation.
  - It would populate He(2³S) directly, which bears on the 10830 Å
    diagnostic.
  - The present tables start at 0.8 keV and do not give the threshold
    behavior, so the size of this contribution is not determined.

- **A coincidence, not validation.** Zygelman's Figure 8 non-radiative
  cross section, extrapolated with its high-energy slope (1.35) from 40 eV to
  805 eV, gives 3.4×10⁻¹⁷ cm². That equals Igarashi's total capture at
  1 keV/u (3.3×10⁻¹⁷).
  - Zygelman's calculation contains only the ground-state channel, so the
    match cannot confirm his cross section.
  - If it is more than a coincidence, it would suggest that his rising
    cross section behaves like a total capture rather than like the ground
    channel alone.
  - The ground-channel value from Loreau at the same energy is 650× lower.
  - This is recorded as an observation only.

## 3. The excitation data

The H(n = 2, 3) excitation and Balmer-α cross sections concern fast
He⁺ + H collisions (keV, e.g. edge plasmas or the solar wind). They do not
enter MoCHII. The ADS abstract states that:
- the H(n = 2) excitation agrees with Geddes et al.;
- the Balmer-α emission is about half of Donnelly et al.;
- the capture is close to Lockwood et al. and Olson et al.;
- the ionization agrees with Hsu et al.

## Limits

- **Taken from the file, not the article.** The energy convention and the
  final states summed in "capture" come from the supplementary file and the
  ADS abstract. If "capture" excludes some states, or if E is in a
  different frame, the ratios in section 2 shift. With a CM energy
  convention they would change by up to a factor of ~2 at 1 keV/u, where σ
  varies steeply.
- **The Loreau fit at these energies.** It is used within its stated range
  (to 10 MeV/u, < 25%).
- **Not a validation.** No rate was changed and no simulation was run.
