# Mg III radiative identity diagnostic

The full bound CHIANTI source has 167 levels, 3,364 positive spontaneous
transition rows, and 13,861 electron-collision pairs. Direct execution of
`audit_mgiii.py` finds 16 spontaneous transitions whose source index order
is opposite to the selected observed/theoretical energy order. All 16 have
opposite-parity configurations, positive source `gf`, and the electric-dipole
angular-momentum selection rule, including exclusion of `J=0` to `J=0`.
The frequency reconstructed from their source A/gf agrees with the original
theoretical gap to at most `0.1798105100310421%`. The JSON preserves each
transition and all three raw-source hashes.

These measurements identify a possible line-strength-preserving correction,
analogous to the [oxygen investigation](../../oxygen_valence_convergence.md).
They do not constitute an adopted full Mg III atom. Selection rules are
necessary conditions for electric-dipole character; source classification and
the convention used for the original radiative coefficients must also be
checked before correcting all rows. A correction must preserve the original
collision scaling coordinates, use one physical level-energy vector in
excitation and emitted photons, and use the physical emitting level's
statistical weight in the spontaneous probability.

The first 27 observed-energy levels have separately sourced Jönsson et al.
(2014) radiative data; they must not be silently replaced by the higher-level
Liang and Badnell (2010) convention. Forbidden noninverted rates need either
an identified multipole correction or an explicit effective-rate approximation
with a stated range. No source or production file was modified by this
diagnostic. It does not measure finite-density populations or the accuracy of
the collision strengths.
