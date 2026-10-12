# Ne I source adoption and validation plan

This directory separates the adopted CHIANTI numerical collision model from
unverified atomic accuracy. No production Ne I atom or cooling coefficients
are installed by this plan. The published 2012 collision Methods and the
2024 CHIANTI conversion method were not verified from the final article.

The original `ne_1.scups` contains 465 type-2 Burgess--Tully transitions,
15 scaled nodes each, with abscissae zero and one at the two endpoints.
Its header cites Zatsarinny and Bartschat (2012), identifies Dufresne's
January 15, 2024 preparation, and says the last three strengths were used
for the infinite-temperature extrapolation. It declares no physical
Maxwellian-temperature accuracy interval. Inversion of its rounded first
finite nodes gives 999.410754--1000.355160 K; the last finite nodes give
998530.651--1001955.908 K. These are inferred numerical nodes, not published
physical validity limits. A common adopted floor of 1000.36 K holds the
effective collision strength constant below that temperature while the
actual gas temperature enters Maxwellian rates and detailed balance. The
intended 1,000--100,000 K use is an explicit numerical approximation.

The raw CHIANTI level, radiation, and collision files have SHA-256 values
`7e813574e3157b9735a7d89660abcf463019e47f1df228ad41ad8c84b6167788`,
`92ae3ab4b0ce76f358617a0ba8d30c0758f0bcfbd7e0529d85ab4eca7509bd8f`,
and `d51fbb34fb71c533a53668e3425a213672ddf66eaf0c7fc51c0b72b6363e4b06`.
The original bounded atom has 31 states, 142 CHIANTI spontaneous rows, 465
collision rows, and 6,975 original scaled collision nodes. The corrected
private model has 144 spontaneous rows after the two identified decays are
added. The existing 401-temperature, four-density native profile is
historical because level 15 was falsely treated as radiatively stable;
the newly adopted temperature floor also changes the first temperature.

## Identified spontaneous radiation correction

CHIANTI supplies no downward spontaneous transition from source level 15,
although its own level list identifies that state as `3p 1S0` at
152970.733 cm^-1, with statistical weight one. Cloudy c25 Stout uses the
same observed energies, weights, and source indices for the first 31 states.
Its NIST-dated September 16, 2014 `.tp` list gives lower/upper level pairs
`3,15` with A = 9.00e5 s^-1 and `5,15` with A = 6.82e7 s^-1. Among its 66
positive A pairs within the first 31 states, exactly these two are absent
from CHIANTI; CHIANTI also has 78 positive pairs absent from Stout. The
original source radiation graph treated levels 1 and 15 as stable. The two
missing decay channels are restored selectively by
`tools/fitting/neon_spontaneous_probabilities.py`; all 142 CHIANTI A rows are
retained. Source hashes, first-31 identity, and exact Stout pair/rate guards
protect this combination. The corrected graph has only the ground state as a
radiatively closed state. This is a necessary physical correction, not a
source-node interpolation adjustment, and it makes the previous cooling fit
and full native profile physically stale. The Stout file identifies its A
source as NIST dated September 16, 2014; this does not verify the separate
2012 scattering Methods or 2024 CHIANTI conversion.

## Conditional calculation sequence

After review of the shared source edits, a persistent controller should run
the following scoped sequence with frozen source and input hashes, finite
CPU, memory, and wall limits, preserved partial products, and an explicit
stop on any failed output check:

1. Generate a private 31-state, 144-radiation-row source-node atom. Verify
   465 collision rows, 6,975 unchanged source nodes, the
   1000.36 K floor record, source hashes, level identities, and radiation
   direction. Keep the previous candidate as a historical baseline.
2. Repeat the affected native source-strength and population/cooling checks
   at 401 temperatures from 1,000 to 100,000 K and electron densities
   1e-6, 1e2, 1e6, and 1e10 cm^-3. Compare all original collision nodes,
   direct effective strengths, full level populations, cooling, and the
   complete curves. Require energy and population balance as applicable.
3. Fit only Ne I low-density coefficients and validate the compiled
   cooling reader at 801 temperatures over the same range. Inspect the full
   fit and native-difference profiles; do not reuse the older first-point
   result if the floor changes it.
4. With the private atom and coefficients in one staged directory, run the
   actual 401-by-161 registered cooling table using an isolated current-source
   optimized Fortran build. Check the 40 exact-grid/midpoint consumer cases,
   structural counts, the physical density-ratio flag, and the existing 0.5%
   midpoint and applied-cooling tolerances. Table temperatures outside the
   adopted 1,000--100,000 K interval test numerical interpolation only.

No unrelated ions, production data, or Fortran objects need alteration for
these checks. Installation follows review of the source and full-profile
results and the resulting numerical differences.
