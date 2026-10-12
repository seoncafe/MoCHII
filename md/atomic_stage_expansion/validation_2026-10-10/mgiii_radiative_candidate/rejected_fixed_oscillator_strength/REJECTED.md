# Rejected Mg III radiative candidate

The preserved candidate used `gf` as if it were an energy-independent dipole line strength. `gf` is an oscillator strength and is proportional to transition frequency times the dipole line strength. Consequently, the candidate scaled E1 spontaneous rates as frequency squared instead of frequency cubed. The numerical parser and source-node checks did not test this physical assumption. This directory is retained only for diagnosis; its atom and results must not be used for cooling or production acceptance.
