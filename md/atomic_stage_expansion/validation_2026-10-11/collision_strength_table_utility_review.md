# Generic collision-strength temperature table correction

The internal Python table function formerly named `clvii_effective_strength`
was used for other ions too. It is now
`effective_collision_strength_temperature_table`; both callers in
`chianti_cooling.py` and the fitting tests use the physical generic name.
No production input parameter or Fortran public routine changed.

The function validates a finite, positive, ordered grid and matching positive
strengths. It permits only eight floating-point spacings of upper-end rounding,
then evaluates the exact endpoint; a genuine extension such as 100000.1 K is
rejected. This does not authorize extrapolation. Five focused unittest cases
passed (0.643 seconds in the recorded direct execution). Unrelated products
were not regenerated. The existing generic error text was already correct;
this change does not claim to correct that text.
