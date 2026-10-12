# Node-shared rate and diffuse tables (2026-10-12)

Binary: the tree after the change; reference outputs (`ref_*`) are from the
binary before it (the `rates_iterate_check` a2 run and the
`refine_rates_check` run). lart1, 16 ranks, `I_MPI_ADJUST_*=1`.

| check | result |
|---|---|
| `a2.in` (id_new, Sobol, 2^3S coupling, 4 iterations) | `_lines.txt`, `_nebcont.txt` byte-identical to `ref_a2_*` |
| `refine.in` (re-refinement at iteration 2) | byte-identical to `ref_refine_*` |
| `emis.in` (`emis_output`), `check_emis_sum.py emis` | 15 metal stages, 66 lines: max \|sum(emis V)/L - 1\| = 3.8e-6 (the 5-digit print of `_lines.txt`) |
| `emis_chunk.in`, a build with `EMIS_CHUNK_VALUES = 1000` (many chunks) | `_emis.h5` 61/61 datasets and `_lines.txt` identical to `emis` |
| `pn150_mem.in` (PN150 L7, 2097152 leaves, 1 iteration), mocafe | node peak 78.8 GB at 40 ranks (`mem_np40.txt`), 140 GB at 72 (`mem_np72.txt`); thermal solve 602 s / 655 s; `pn150_mem_np40_lines.txt` and `pn150_mem_np72_lines.txt` identical in every printed value |

Before the change the PN150 worker measured 191 GB at 40 ranks on lart3 and
an out-of-memory stop at 56. Rates and emissivity files were deleted after
the comparison.
