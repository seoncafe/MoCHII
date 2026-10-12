# par%rates_write_interval and the resume it serves (2026-10-12)

Inputs are `tests/cooling_ne/id_new.in` with `gas_niter` and the options
below; tree binary, lart1, 16 ranks, `I_MPI_ADJUST_{ALLREDUCE,REDUCE,BCAST,GATHER}=1`
(`run.sh`, `run2.sh`). Rates files were deleted after the comparison.

| run | input | result |
|---|---|---|
| a | `gas_niter = 4`, interval 1, `launch_sequence = 'random'`, coupling off (the input's pins) | iterate written after iterations 1-3, header NITERDON 3, CONVERGD 0, STATESOL 0; no temporary file left; `_rates` at the end, STATESOL 1 |
| b | restart from a's iterate, `gas_niter = 1` | runs; its pass differs from a's iteration 4 by 1.1e-4 in absorbed energy: the random generator state is not stored, so no restart continues a `'random'` run |
| n | interval -1 | no iterate written, normal exit |
| a2 | as a with the default Sobol launch and He I 2^3S coupling | iterate after 1-3, carries `Phi_HeI_23S` |
| b2 | restart from a2's iterate, `gas_niter = 1` | pass and step identical to a2's iteration 4 in every printed digit; `_rates` 33/33 datasets bit-identical, `_lines.txt`, `_nebcont.txt` byte-identical |

b2 needed `Phi_HeI_23S` in the file: before it the restart took Phi_3 = 0
for the 2^3S split of the first pass's He I diffuse channel.
