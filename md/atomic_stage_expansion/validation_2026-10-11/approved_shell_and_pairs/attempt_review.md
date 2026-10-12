# Preserved diagnostic attempt

The full production build passed. The consumer attempt correctly rejected
20 eV when inner-shell absorption was explicitly enabled: Fe I has a 14.7 eV
inner edge. The diagnostic incorrectly assumed that all first three test
energies were below every inner edge. This was a fixture error, not a failure
of the production guard. The corrected diagnostic uses 10, 12, and 14 eV,
which are below the measured minimum inner edge in both rate rosters.

The corrected fixture uses the independently checked five-level S II atom
with unchanged charge-transfer term identities. This avoids repeatedly
building an unrelated 70-level population table; the measured cross sections
do not depend on those populations. Failed logs and driver sources are
preserved here. The continuation uses a new output directory.
