# Provenance of charge-transfer heat for nine added boundaries

## Physical finding

The nine positive heat values in Cloudy c25.00 are published in Table 1 of
Kingdon and Ferland (1999, *ApJ Letters* **516**, L107-L109,
[doi:10.1086/312008](https://doi.org/10.1086/312008)). The
[publisher's final PDF](references/Kingdon_Ferland_1999_ApJ_516_L107.pdf)
was checked directly with `pdftotext -layout`, including Section 2 and Table 1.
They are the authors'
**adopted energy defects** for `X^q+ + H(1s) -> X^(q-1)+ + H+`. The paper defines
positive `Delta E` as kinetic energy added to the gas, with heating and cooling
measured relative to the continuum. It explicitly uses total charge-transfer
rates and mean energy defects when several product channels exist. Thus the
table supplies a physically motivated thermal prescription, but does not give
state-specific branching fractions or identify a unique product level for each
of these nine reactions. Neither the Cloudy table nor the MAPPINGS copy should
be interpreted as a ground-state ionization-potential difference.

| Reactant | Cloudy zero-based `ion` | KF99 Table 1 / Cloudy `Delta E` (eV) | MAPPINGS `e4` | MAPPINGS heat from `e4` (eV) |
|---|---:|---:|---:|---:|
| C^4+ | 3 | +11.30 | 13.113 | 11.299909 |
| N^3+ | 2 | +6.40 | 7.4269 | 6.400007 |
| O^3+ | 2 | +5.00 | 5.8023 | 5.000035 |
| O^4+ | 3 | +8.47 | 9.8290 | 8.469977 |
| Ne^3+ | 2 | +5.82 | 6.7538 | 5.819975 |
| Ne^4+ | 3 | +8.60 | 9.9799 | 8.600012 |
| S^4+ | 3 | +8.60 | 9.9799 | 8.600012 |
| Ar^4+ | 3 | +8.60 | 9.9799 | 8.600012 |
| Mg^3+ | 2 | +5.73 | 6.6494 | 5.730010 |

The MAPPINGS column was evaluated directly as
`e4 * k_B * 10,000 K / (1.602176634e-12 erg/eV)`, with its source constant
`k_B=1.380649e-16 erg/K`. The near identity after rounding demonstrates shared
tabulation, not an independent product-state calculation. In particular, the
repeated 5.73 and 8.60 eV values are consistent with the paper's generic
avoided-crossing estimate; the paper does not flag, row by row, which of these
nine entries comes from state-resolved literature and which from that estimate.
Their individual derivations therefore remain unverified.

## Published method and limits

Section 2 of Kingdon and Ferland (1999) defines the volumetric thermal term as
`n(H) n(X^q+) alpha_CT Delta E`. The authors assume ground-state reactants and
interstellar conditions. When channel-specific results exist, their mean defect
uses the individual channel defects; when they do not, they infer a mean defect
from the assertion that strong one-electron transfer occurs at avoided crossings
roughly 7-13 bohr apart, adopting 9.5 bohr. They state that defects inferred in
this manner should generally be accurate to substantially better than 50%; this
is their estimate, not a verified uncertainty for each row. Table 1 lists the
adopted `q<=4` defects. The same paper treats radiative transfer from product
cascades as outside its direct thermal term and neglects secondary ionization by
those photons and collisionally modified cascades at very high density. These
assumptions delimit use of the values as approximate staged heat records; they
do not by themselves establish the thermal treatment of product photons in
MoCHII.

The paper's prose says the mean defect for `q>4` is `2.86 q eV`, with `q`
defined earlier as the **reactant** ion charge. The actual Cloudy code uses
`2.86*ion` for its zero-based product-charge index `ion=q-1`, yielding
`2.86*(q-1) eV`. This is also the Coulomb energy of the product ions at the
adopted 9.5-bohr crossing, `(q-1) * 27.2114/9.5 eV`, and the approximate
pattern in Table 1 is 5.73 eV for `q=3` and 8.60 eV for `q=4`. Consequently,
the paper's written `2.86 q` and the implemented/physical `2.86*(q-1)` are
inconsistent. The original intention of that one symbol cannot be settled from
these sources; the high-charge rule is an estimated capture heat and is **not**
one of the nine tabulated `q<=4` values.

## Executed code and data paths

In Cloudy c25.00, `ChargTranEval` in `source/atmdat_char_tran.cpp:44-79`
populates `CharExcRecTo[H][nelem][ion]` through `HCTRecom(ion, nelem)`.
`HCTRecom` initializes `MakeHCTData` on its first call, and that reader opens
`data/ctrecombdata.dat`, checks magic number `201903042`, and reads seven
numbers for each of four rows per element into
`CTRecombData[nelem][ion][0..6]` (`:806-868`, `:919-941`). The element index
is atomic number minus one, and `ion` is product charge `q-1`. For `ion<=3`,
the first four columns determine a total rate in `cm^3 s^-1`; the fifth and
sixth bound the temperature used for the fit. The seventh is the signed heat
defect in eV; it is **not** used to evaluate the rate. For `ion>3`, `HCTRecom`
returns the configurable high-charge coefficient instead of reading a fit.

`ChargTranSumHeat` at `:569-670` adds
`alpha_CT * CTRecombData[nelem][ion][6] * n(H^0) * n(X^(ion+1)+)`
for these rows. Each positive value in the table above therefore adds gas heat.
It converts the sum from `eV cm^-3 s^-1` to `erg cm^-3 s^-1` using `EN1EV`,
then applies the `HCharHeatOn` switch. `conv_base.cpp:909-912` calls this
routine after ionization convergence and sorts a negative net result into
cooling; `heat_sum.cpp:761` consumes positive charge-transfer heat.
`ion_cx.cpp:35-43` separately consumes the rate coefficients in the ionization
balance. This traces both the rate and heat paths, rather than merely finding a
table entry.

MAPPINGS V-521 reads `lab/data/ionisation/CHXDAT.txt` in `src/mapinit.f:3135-3210`.
The reaction rows hold the same rate fits and a final dimensionless `e4` defect,
defined in the file as `Delta E/(k_B*10,000 K)`. `src/charex.f:160-201`
computes the recombination coefficient; `src/cheat.f:41-60` computes
`chxrde * k_B * 10,000 K`, multiplies by that coefficient and both reactant
densities, and adds the result to `chgain`. This source supplies no independent
state-specific product data. An older local Cloudy c07.02.02
`source/atmdat_char_tran.cpp` already contains the same `2.86*ion` high-charge
formula and embedded low-charge defect data, establishing that these conventions
predate the c25.00 external data file.

## Source and verification record

| Directly inspected source | SHA-256 |
|---|---|
| `md/atomic_stage_expansion/references/Kingdon_Ferland_1999_ApJ_516_L107.pdf` | `2c681fe32532baddee5d19a67e472cefc51f6a0135de01ab82eb9aa182e78def` |
| `/home/kiseon/CLOUDY/c25.00/source/atmdat_char_tran.cpp` | `22a1443eed27765841679a705f53120bb25fbbf8d35757942e2096c73ce605c9` |
| `/home/kiseon/CLOUDY/c25.00/data/ctrecombdata.dat` | `4ff7f4aeaf24dd5e0c6ab450694e480ef16136fab443c4c68b7acd799c5cbe45` |
| `/home/kiseon/RT_Codes/MAPPINGS/mappings_V-521/lab/data/ionisation/CHXDAT.txt` | `b239065f0d1d974471806fd0c2b22b9164df3b2aefd6bbbbe97e8ec3f18fd69b` |
| `/home/kiseon/RT_Codes/MAPPINGS/mappings_V-521/src/mapinit.f` | `dc1ba60ab10b50f3d4592c69f2f31c48a8cc7b049ec411f6fa2f6bd2b1d0d715` |
| `/home/kiseon/RT_Codes/MAPPINGS/mappings_V-521/src/charex.f` | `423353177343b07609fbf1c5caa6b13827453954e3982973d071149691ba404d` |
| `/home/kiseon/RT_Codes/MAPPINGS/mappings_V-521/src/cheat.f` | `066abbaf5690cddc0704dd925a7cf65fec3e3affeca90c955c2d90c81a434289` |
| `/home/kiseon/CLOUDY/c07.02.02/source/atmdat_char_tran.cpp` | `ccf60d9ba3aedcf2c6c9ecb6106ea90b2d5117085eabb258bccec58fa7abfd6c` |

The ADS API was attempted before web search, but DNS lookup failed in this
shell. The original final journal PDF was subsequently obtained from the
[IOP publisher](https://iopscience.iop.org/article/10.1086/312008/pdf), and its
three pages were read with `pdftotext -layout`. The nine Table 1 entries and
Section 2 method were checked against that PDF, not inferred from the 1996
rate-fit paper or Cloudy's comments. The paper does not provide a row-by-row
source for the nine individual defects, so their specific product-channel
derivations remain unverified. No production simulation, production data, or
external source file was changed in this provenance investigation.
