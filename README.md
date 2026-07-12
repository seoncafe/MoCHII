# MoCHII

**MOnte Carlo for H II regions** — Monte Carlo radiative transfer in dusty
photoionized nebulae on adaptive (octree AMR) grids.

One self-consistent radiation field drives both the gas and the dust:
H/He photoionization and thermal balance, metal lines added one ion at a
time (CHIANTI-fitted rates), recombination lines, the nebular continuum
(free-bound / free-free / two-photon), and dust + PAH emission.

Built on the validated AMR engine of MoCafe v2.00; gas microphysics follows
the analytic-fit approach of EXHALE; MOCASSIN serves as the 3D validation
reference. Author: Kwang-il Seon (KASI).

- Design and staged plan: `docs/PLAN.md`
- Module provenance: `src/PORTING.md`
- Guidance for coding sessions: `CLAUDE.md`

Status: skeleton (2026-07-11). First milestones: G0 (ionizing bins + rate
integrals) and the CHIANTI fitting pipeline under `tools/fitting/`.
