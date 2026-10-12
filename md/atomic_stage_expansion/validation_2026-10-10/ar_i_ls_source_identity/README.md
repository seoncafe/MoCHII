# Neutral argon LS-source identity audit

**Judgment:** the local OPEN-ADAS ADF04 file is structurally sufficient to
identify one 35-term LS electron-impact and spontaneous-emission system, but it
is not a physically accepted cooling atom as supplied. Its independent
radiative sink at term 2 prevents use of the ordinary ground-only coronal
limit; this is a limitation of an LS model, not by itself proof of missing
physically allowed radiative transitions. Two collision pairs have
extraordinary strengths. In particular, the
8-to-15 effective collision strength at 10,000 K gives an excitation rate of
19.645 cm^3 s^-1 under the standard Maxwellian formula. This requires source
resolution before translating the file into production data. The LS terms
must not be paired by index with Stout's 503 fine-structure states.

## Inspected source and format

The read-only source is `/tmp/arlike_cpb12_ar0.dat`, downloaded from the
[OPEN-ADAS ADF04 listing](https://open.adas.ac.uk/element/detail/adf04/arlike/arlike_cpb12%5D%5Bar0.dat)
on October 9, 2026. Its SHA-256 is
`daa65074f68d9e0ad626995e75cb4ef54fe971f2fe45979b6a54a421ae8e4128`.
The file footer attributes its scattering calculation to C. P. Ballance and
D. C. Griffin, *Journal of Physics B* **41**, 065201 (2008), DOI
[10.1088/0953-4075/41/6/065201](https://doi.org/10.1088/0953-4075/41/6/065201),
and says Connor Ballance prepared this record in December 2012. Its orbital
description says `2p5nl,2s2p6nl`, whereas the actual level labels are `3p5nl`;
the mismatch is an unresolved provenance error, not a reason to relabel the
source. The file does not separately identify the origin of its A values.
The final publisher article and its methods were not read in this audit, so
published accuracy and intended temperature restrictions remain unverified.

The [reproducible inspector](audit_source.py) parses the source without copying
its table into this repository. Its [measured results](result.json) show:

| Property | Direct result |
|---|---:|
| LS terms | 35, numbered 1–35 |
| Highest term energy / file ionization boundary | 119,843 / 127,110 cm^-1 |
| Statistical weight | `g = 2 XJ + 1 = (2S+1)(2L+1)` for every term |
| Electron-impact rows | 595 distinct pairs, exactly all `35 choose 2` |
| Temperature nodes | 14, from 100 to 2,000,000 K; includes 1,000–100,000 K |
| Radiative rows | 104 positive A; 491 zero A |
| Collision-node signs | All 8,330 tabulated strengths finite and positive |
| Bethe column | One trailing optional high-energy value, not a 15th temperature |

The 35 listed term energies are strictly below the file's ionization boundary;
no explicitly autoionizing state or ionization/recombination row occurs in its
35-term electron matrix. Five collision pairs have equal *rounded* term
energies: `(7,8)`, `(15,16)`, `(21,22)`, `(27,28)`, and `(34,35)`. Their
collision strengths are finite and positive. No positive A row has a zero or
negative listed transition energy. The file represents complete collision
connectivity at every tabulated temperature, but this does not establish that
all physical spontaneous channels are included.

## Physical limitations found in the file

Term 2 has no positive A to any lower term; terms 1 and 2 are both radiative
sinks. All other terms have an A path to either sink, but 17 terms, including
term 2, have no A path to the ground term 1. Term 2 is listed as a triplet
`3p5 4s` term above the singlet ground term; a spin-forbidden electric-dipole
decay is absent in strict LS coupling. Thus two sinks can be a coherent
approximation rather than a corrupted transition list. At finite positive
electron density the complete collision matrix can couple both sinks, but the
usual ground-only low-density cooling limit is not established. Actual weak
spin-forbidden and other radiative losses, if required at the intended
density, need an identified physical source rather than an inferred zero or
an arbitrary replacement rate. This LS source must not be silently presented
as a complete fine-structure radiative model.

The raw type line begins `1.00    3`, where `3` is the ADF04 type for
Maxwellian Upsilon versus electron temperature. Each transition line contains
two indices, A, 14 Upsilon ordinates, and an optional Bethe-limit value.
The original transition 15-to-8 source line has 18 whitespace-separated
ASCII fields: two indices, A, 14 strengths, and the Bethe value. The
`1.93+09` token is the seventh Upsilon ordinate at 10,000 K; it is
neither a negative exponent nor an A value nor the trailing Bethe limit.
Two effective collision strengths exceed `1e5` at 10,000 K. The most extreme
is the 8-to-15 row: listed energy gap 7,220 cm^-1, lower LS weight 3, and
Upsilon `1.93e9` at 10,000 K. The Maxwellian excitation expression
`8.629e-6 Upsilon exp[-Delta E/(k_B T)]/(g_lower sqrt(T))` gives
`19.645 cm^3 s^-1` at 10,000 K and `1.308e-5 cm^3 s^-1` at 1,000 K using
that row's respective nodes. These values are far beyond plausible
electron–atom excitation rates. They are reproduced by direct parsing of the
file's printed exponents; they are not interpolation overshoot. The ADF04
numeric syntax and 14-node/one-Bethe-column distinction were checked against
the locally available [official ADF04 specification](https://www.adas.ac.uk/man/appxa-04.pdf).
The source of the anomalous strengths is not established, and a guessed
exponent correction would not be justified.

The file's 100–2,000,000 K grid is a numerical tabulation range, not an
independently verified physical accuracy range. The source includes no
fine-structure line-resolved rates. A consistent LS calculation could in
principle approximate term-summed cooling after the anomalous collision rows
and radiative closure are resolved, but it cannot predict individual Ar I
fine-structure line emissivities. The OPEN-ADAS footer limits the download to
personal use and prohibits redistribution or inclusion with modeling code
without express written permission. This audit retains only derived checks;
the source table remains outside the repository.

Before numerical adoption, identify the published or author-provided meaning
of the exceptional Upsilon rows, establish the provenance and applicability of
the LS A values (particularly at densities where metastable-term losses
matter), and confirm the source
accuracy range and data-use permission. A candidate using a corrected,
internally consistent source would then need a complete finite-density and
low-density radiative-balance check, exact source-node verification, 1,000–
100,000 K cooling and population profiles, native-reader comparison, and a
source-preserving registration test. None of those cooling calculations was
run here because the present physical inputs are not acceptable.
