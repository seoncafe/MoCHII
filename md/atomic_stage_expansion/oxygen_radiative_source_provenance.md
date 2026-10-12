# Ground-configuration radiative probabilities in the installed atoms

This audit separates a measured source difference from a proved error in a
transition probability. The local final printed article is
`/home/kiseon/MoCHII/references/Storey_Zeippen_2000_MNRAS_312_813.pdf`
(SHA-256 `0e11e93af007e3077b1b0da570ca45c7376f0fe756e2892c52bfea86f2ed852d`).
The final citation is Storey and Zeippen, *MNRAS* **312**, 813–816 (2000),
[doi:10.1046/j.1365-8711.2000.03184.x](https://doi.org/10.1046/j.1365-8711.2000.03184.x).
Section 2.3 and Tables 1–4 were read in the final printed version. The paper
applies order-alpha-squared corrections to the one- and two-electron magnetic
dipole transition operator. Its Table 3 and Table 4 totals include E2 and the
corrected M1 probability; the earlier Galavis et al. (1997) M1 probabilities
used a different operator. These corrections are separate from relativistic
terms in the atomic Hamiltonian or improvements to wavefunctions.

The indices below are the actual selected levels in the current atomic files.
In the carbon sequence, `2` is 3P1, `3` is 3P2, and `4` is 1D2. In the oxygen
sequence, `1` is 3P2, `2` is 3P1, and `4` is 1D2. Each entry is a total decay
probability in s^-1. Values were read from the `NRAD` records, not inferred
from a source-file header.

| Ion and actual file | 1D2 to 3P1: installed / Storey–Zeippen | 1D2 to 3P2: installed / Storey–Zeippen | Source status |
| --- | ---: | ---: | --- |
| C I, `data/atomic/nlevel_c_1.txt` | 6.17e-5 / 7.154e-5 | 1.81e-4 / 2.240e-4 | Current CHIANTI ground forbidden rates derive from Hibbert et al. (1993) CIV3; the two values are substantially different from the corrected Table 3 totals. |
| N II, `data/atomic/nlevel_n_2.txt` | 9.849e-4 / 9.819e-4 | 2.919e-3 / 3.015e-3 | CHIANTI identifies Tachiev and Froese Fischer (2001) for levels 1–5. The rate ratio is close to the earlier uncorrected result, but this is not an operator-level proof. |
| O III, `data/atomic/nlevel_o_3.txt` | 6.791e-3 / 6.791e-3 | 2.046e-2 / 2.046e-2 | The installed pair now equals final Table 3. The prior atom is preserved in `validation_2026-10-10/oxygen_radiative_provenance/nlevel_o_3_before_correction.txt`. |
| Ne V, staged `validation_2026-10-10/source_density_ratio/NeV/nlevel_ne_5.txt` | 1.324e-1 / 1.222e-1 | 3.702e-1 / 3.504e-1 | This staged 25-level atom is not an installed `data/atomic/nlevel_ne_5.txt`. It derives from the CHIANTI ground rates attributed to Tachiev and Froese Fischer (2001). |
| O I, `data/atomic/nlevel_o_1.txt` | 2.101e-3 / 2.151e-3 | 6.503e-3 / 6.446e-3 | CHIANTI identifies Froese Fischer and Tachiev (2004). The rate ratio equals the earlier Table 4 value to rounding. |
| Ne III, `data/atomic/nlevel_ne_3.txt` | 5.107e-2 / 5.417e-2 | 1.653e-1 / 1.711e-1 | CHIANTI identifies Landi and Bhatia (2005). The rate ratio equals the earlier Table 4 value to rounding. |

For C I, the later final Froese Fischer (2006) paper,
[doi:10.1088/0953-4075/39/9/005](https://doi.org/10.1088/0953-4075/39/9/005),
is archived as
`/home/kiseon/MoCHII/references/Froese_Fischer_2006_JPhysB_39_2159.pdf`
(SHA-256 `154d94abb23f070096ee3140e42ec1c860946f0abb8527076bfea792a2408533`,
final printed PDF from the NIST MCHF collection). It
is especially useful because Table 7 displays separate E2 and M1 contributions
from the earlier CIV3 calculation and its newer MCHF calculation. The CIV3 M1
probabilities are 6.17e-5 and 1.836e-4. The installed first branch equals
that CIV3 value; its second branch is 1.81e-4 in both the installed atom and
its CHIANTI source, 1.4% below the value quoted in the 2006 comparison.
The CHIANTI header identifies both installed branches with the earlier CIV3
source, but this numerical difference remains unexplained. The newer MCHF M1
probabilities are 7.54e-5 and 2.24e-4, but this paper does not establish in
its methods that Storey–Zeippen's corrected M1 operator was included.
Its Table 3 explicitly uses *ab initio* energies. Thus neither numerical
agreement of one newer value nor publication date alone settles the choice.
The [Wiese and Fuhr 2007 NIST critical compilation](https://www.nist.gov/document/jpcrd3620071287ppdf),
archived as `/home/kiseon/MoCHII/references/Wiese_Fuhr_2007_JPCRD_36_1287.pdf`
(SHA-256 `26139ce2dba92e0c32c7af7b1e4c5d8fd1ef75de91dd05a8d29a9e87ef80f10d`,
final printed PDF from NIST), recommends energy-adjusted C I totals in
Section 5.1.2 and Table 6 of approximately
7.295e-5 and 2.1807e-4, differing from the corrected Table 3 pair. It
attributes these to the Froese Fischer (2006) calculation. This is a genuine
source disagreement, not evidence that the Storey–Zeippen correction is absent
from every later calculation.

For N II, the same NIST compilation, Section 6.2.2 and Table 17, recommends
separate M1 and E2 contributions: 9.8e-4 plus 9.2e-7 for 6548 A, and 2.9e-3
plus 8.7e-6 for 6583 A. It attributes the calculation to Tachiev and Froese
Fischer (2001). The table assigns both lines uncertainty class D. The final
2001 article's transition-operator details have not yet been obtained here;
therefore the installed-vs-corrected difference must remain a documented source
choice pending direct methods verification. The same applies to staged Ne V.

For O I, the final Froese Fischer and Tachiev (2004) compilation,
[doi:10.1016/j.adt.2004.02.001](https://doi.org/10.1016/j.adt.2004.02.001),
is archived as
`/home/kiseon/MoCHII/references/Froese_Fischer_Tachiev_2004_ADNDT_87_1.pdf`
(SHA-256 `a5b5de7b87429b509fcd900acc907e92d155b7887490720adac0a4a1e1d2d971`,
final printed PDF from the NIST MCHF collection). Its Section 3.7 includes
M1 and E2 transitions and validates selected lifetimes.
Its discussion of Breit–Pauli corrections to wavefunctions and E1 length/velocity
forms does not state that the one- and two-electron M1 transition-operator
correction of Storey and Zeippen was used. For Ne III, the source is actually
Landi and Bhatia, *Atomic Data and Nuclear Data Tables* **89**, 195–265 (2005),
[doi:10.1016/j.adt.2005.02.003](https://doi.org/10.1016/j.adt.2005.02.003);
the CHIANTI `.wgfa` header incorrectly gives page 139, which belongs to a
different paper in the same volume. The available final abstract says the
SUPERSTRUCTURE code computed radiative rates, but its full methods and M1
operator choice remain unverified. No O I or Ne III rate replacement is
justified by the ratio comparison alone.

The original full-source O II and O III radiation had a separate energy-order
defect: 31 O II rows and 48 O III rows reversed under the selected hybrid
energies. A rates have a defined *upper* level and cannot be repaired by
merely exchanging indices. The currently installed full 196-level O II and
177-level O III atoms use source-identified opposite-parity E1 line strengths
at their selected physical frequencies; the O II 2D fine-structure M1
direction follows its observed level energies and an independent Stout
identity check. The original CHIANTI files remain unchanged. The correction,
conditional treatment of other forbidden rates, full profile validation,
and prior installed file hashes are recorded in
`oxygen_valence_convergence.md` and its validation directory.

The O III correction was tested on the actual six-level native consumer at
3,204 temperature-density points; see
`validation_2026-10-10/oxygen_radiative_provenance/native_main/result.json`.
The 5007/4959 energy-flux ratio is 2.9839683 with the corrected rates. The
largest measured cooling difference was 0.02906775% over the tested grid.
Its low-density cooling is unchanged within numerical precision because every
excitation to this radiatively connected model eventually releases its energy
to the ground state. No cooling fit regeneration follows solely from these
ground-pair A changes. A future change to other ions still requires a direct
consumer check, including density-sensitive cooling and line ratios.

No C I, N II, Ne V, O I, or Ne III production radiative rate was changed by
this audit. Before selecting a replacement, verify the final source methods
for the M1 operator, transition energies, and any experimentally constrained
lifetimes or branching ratios; then compare the published absolute A values
for both branches at the same level identities. Preserve the original source
files for a direct numerical comparison.
