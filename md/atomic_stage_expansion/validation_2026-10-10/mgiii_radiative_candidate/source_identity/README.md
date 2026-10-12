# Mg III observed-level identity check

The file `mgiii_nist_levels_2026-10-10.tsv` is the unmodified ASCII response
retrieved on 2026-10-10 from the [NIST Atomic Spectra Database Mg III level
query](https://physics.nist.gov/cgi-bin/ASD/energy1.pl?spectrum=Mg%20III&units=0&format=3&output=0&page_size=500&multiplet_ordered=0&conf_out=on&term_out=on&level_out=on&j_out=on).
Its SHA-256 is
`2ef56d303565c720ecacc1254f0402850d806229c4d6e79c7dc6398633506084`.
The source CHIANTI `mg_3.elvlc` SHA-256 and row-level comparisons are in
`result.json`. The script only reads these files and does not change either
source.

NIST lists 111 Mg III levels. Of 68 CHIANTI levels above 27 with an observed
energy, 67 have a NIST energy within 1 cm^-1 with the same `J` and the same
outer orbital configuration. NIST uses intermediate-coupling labels for many
of those levels, so this comparison does not equate the term labels or verify
wavefunction mixing.

CHIANTI level 167 is `2s2 2p5 7p 1S`, `J = 0`, with theoretical energy
634430 cm^-1 but an assigned observed energy of 628148.375 cm^-1. The NIST
energy 628149 cm^-1 belongs to a `7d`, `J = 1` level. Its parity and
statistical weight differ from CHIANTI level 167. The observed-energy
assignment is therefore physically inconsistent with that level identity;
the earlier private Mg III candidate using it is retained as a stale
conditional numerical diagnostic.

CHIANTI level 127 has no observed energy. Its selected energy of 620693
cm^-1 is theoretical. Level 143 has an observed energy of 620701.688 cm^-1,
matching the NIST `6d`, `J = 1` level at 620702 cm^-1. Their selected energy
gap of 8.688 cm^-1 is not an independently observed interval. Neither this
NIST comparison nor the numerical balance ledger establishes the accuracy of
the collision strengths or dipole matrix elements.
