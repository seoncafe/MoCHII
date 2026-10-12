# Final gas solve after the final pass (2026-10-12)

`src/main.f90`: after the loop the last solved state is transported once
more (the imaging pass under `par%ion_peel`) and the gas is solved once more
in that field, so the written H, He, T_e, n_e and metal stages are the
solution of the written field. Header: `STATESOL` (replaces `FLDCONS`),
`FINALDX`/`FINALDTE` = the final solve's step.

- `bb1e5_it2.in`: the Phase 5 attribution input `it2/d` (1e5 K blackbody,
  HII40 geometry, C N O Ne S, 2 iterations from the initial guess); lart3,
  16 ranks, collectives pinned (`run.sh`). Its `atomic_dir`/`data_dir`
  pointed at the attribution data snapshot (19:33 KST, 2026-10-11), deleted
  after this run.
  - Iterations 1-2 print the same steps as `it2/d` to all digits.
  - Final solve: max |dx_HII| 8.60e-2, max |dTe|/Te 0.516, vol |dx_HII|
    1.23e-2 (the run is not converged; this is how far the written state is
    from the state the written field was transported through).
  - `check_written_charge.py`: written n_e against the charge of the
    written H, He and metal stages (n_H = 100). Before (`it2/d`): max
    3.67e-4, volume mean 1.37e-4, 136380 of 136464 gas leaves above 1e-6.
    After: max 2.49e-11, mean 1.54e-11, none above 1e-6.
  - Header: `STATESOL` 1, `FINALDX` 0.08597, `FINALDTE` 0.5159,
    `NITERDON` 2.
- `fixed_te/id_fixed_te.in`: `tests/cooling_ne/id_new.in` with
  `solve_te = .false.`, lart1, 16 ranks: the fixed-Te final solve
  (`gas_equilibrium_update`) runs and prints its step (exit 0).
- `MoCHII_final_solve.x`: the tree binary of 2026-10-12 07:00 KST.
