# tests

Regression cases + baselines, MOCASSIN-suite style.  Each directory is
self-contained: grid (or its builder), input file(s), and a `check_*.py`
gate script.  Gates per stage are defined in docs/PLAN.md section 5;
results are recorded in CLAUDE.md and docs/MoCHII_physics.pdf.

| directory | gate |
|---|---|
| `g0_gamma` | rate integrals vs analytic attenuation |
| `g1_stromgren` | Stromgren sphere vs the 1D Gauss-Seidel reference; AMR = uniform |
| `g2a_thermal` | H/He thermal balance vs the 1D reference |
| `g2_hii` | MOCASSIN Lexington HII20/HII40 + line fluxes; element smokes (`smoke_*.in`), emissivity output (`smoke_emis.in`), He I excited channel (`hii40_hei.in`) |
| `g4_refine` | I-front re-refinement vs native level 7 (+ window recycling) |
| `g4_tng` | Illustris-TNG post-processing demo |
| `d_dusty` | dusty Stromgren, scattering bracket, T_dust/IR, SEDust smokes (`sedust_smoke.in`, `dustemis_smoke.in`, `pahlive_smoke.in`), FUV option (`d_fuv_*.in`, `check_fuv.py`) |
| `peel` | peel-off imaging: optically-thin analytic direct gate (`peel_thin.in`, `check_peel.py`), scattered morphology (`peel_scat.in`), bin cubes (`peel_cube.in`) |
| `uni_dda` | Cartesian (car) grid vs single-level octree: shared walk bit-identical; incremental DDA walk (`car_walk='dda'`) to rounding + timing (`check_uni.py`); namelist-built car grid (`car_namelist.in`) + `nH_const`/`rmax` density model reproduce the file-built sphere (`amr` override `car_amr_ovr.in`) |
| `pdr` | PDR physics: metal electrons, photoelectric heating, H-impact cooling (`pdr_L5.in`, `check_pdr.py`) |
