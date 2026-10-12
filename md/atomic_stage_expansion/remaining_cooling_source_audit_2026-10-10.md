# Remaining Ar I, Ar II, Fe I, Ca III, and Ca IV cooling source audit

Date: October 10, 2026. This audit inspected local source files and the code that reads and evaluates them. A file containing radiative rates alone does not supply electron-impact cooling. At the initial file check, `data/atomic/cooling_ca_4.txt` existed; `cooling_ar_1.txt`, `cooling_ar_2.txt`, `cooling_fe_1.txt`, and `cooling_ca_3.txt` did not. Ar II and Ca IV were subsequently installed and validated as recorded below; this initial inventory is historical.

## Directly usable local electron data

| Ion | Local Stout source, under `/home/kiseon/CLOUDY/c25.00/data/stout/` | Parsed levels / positive A rows / electron collision rows | Measured source temperature nodes | Assessment for 1,000–100,000 K numerical use |
|---|---|---:|---|---|
| Ar I | `ar/ar_1/ar_1.{nrg,tp,coll}` | 503 / 428 / 0 | None | `coll` says `g-bar` only. No tabulated electron effective collision strengths. |
| Ar II | `ar/ar_2/ar_2.{nrg,tp,coll}` | 418 / 307 / 1 | 1,580–158,000 K, 11 nodes | A coherent two-level fine-structure atom is feasible. The adopted constant first-node Upsilon would cover 1,000–1,580 K; that interval is an approximation, not a source validity claim. |
| Fe I | `fe/fe_1/fe_1.{nrg,tp,coll}` | 846 / 2,542 / 0 | The `coll` file has nine H and proton rate nodes, 30–3,000 K, but no electron collision row. | H/proton data cannot be substituted for electron effective collision strength. |
| Ca III | `ca/ca_3/ca_3.{nrg,tp,coll}` | 181 / 530 / 0 | None | No tabulated electron effective collision strengths. |
| Ca IV | `ca/ca_4/ca_4.{nrg,tp,coll}` | 70 / 3 / 1 | 10,000–1,000,000 K, 11 nodes | A coherent two-level fine-structure atom is feasible with the adopted constant first-node Upsilon for 1,000–10,000 K. The cold interval is an approximation; resonance accuracy there is unverified. |

Counts above are direct parses: numeric `.nrg` records, `A` `.tp` records, and `CS ELECTRON` `.coll` records. They are source counts, not selected bound atom sizes. The Ar II 1→2 levels have statistical weights 4 and 2, energies 0 and 1431.583 cm^-1, and A(2→1) = 0.0530 s^-1. Their 1/1431.583 cm wavelength is 6.9853 micrometers, identifying the approximately 6.98-micrometer line. Ca IV's corresponding weights are 4 and 2, energies 0 and 3118.200 cm^-1, and A = 0.543 s^-1. For each, the sole electron row is exactly 1→2. The full higher-level radiative set cannot be excited with this collision source; a larger atom would require additional electron data. The published calculation restricts each ion to its own reliable temperature range in Table 3; the local tabulated intervals must not be interpreted as a broader accuracy claim.

The local CHIANTI v11.0.2 tree has no `.elvlc/.wgfa/.scups` excitation triplet for any of these five ions. Ca III and Ca IV have `.easplups/.easplom` files, which do not provide ordinary bound-bound electron excitation rows. Ar I and Ar II have only auxiliary CHIANTI files; Fe I has only `.diparams`. This is a file inventory, not a claim that no published collision calculation exists.

## Independently invoked Ar II and Ca IV source

The local MOCASSIN v2.02.73.2 files `/home/kiseon/RT_Codes/MOCASSIN/mocassin-mocassin.2.02.73.2/data/arii.dat` and `data/caiv.dat` are actual two-level excitation tables, not unused reference files. Their headers identify Pelan and Berrington (1995, A&AS 110, 209) for Upsilon; Ar II has 11 nodes from 1,585 to 158,489 K and Ca IV has 11 nodes from 10,000 to 1,000,000 K. The table headers attribute A values and level energies separately to NIST. The corresponding files under `MOCASSIN-3.0/data/` have identical SHA-256 values, but this does not establish that the later tree invokes them.

For the v2 tree, `data/fileNames.dat` lines 131 and 153 select `data/arii.dat` and `data/caiv.dat`. `source/hydro_mod.f90:499–672` reads those names, confirms file existence, allocates the level/temperature/rate arrays, and loads the tables. `source/hydro_mod.f90:711–730` marks ions with extant files as available. `source/update_mod.f90:1268–1312` calls `equilibrium` for available ions; lines 1100–1109 add the resulting forbidden-line power to the collisional cooling sum. `source/emission_mod.f90:1913–2122` interpolates collision strengths on log10(T), constructs electron excitation/de-excitation rates and solves populations. Its interpolation/extrapolation is a property of MOCASSIN, not evidence that the Pelan–Berrington scattering calculation is accurate outside its tabulated interval. `ari.dat`, `fei.dat`, and `caiii.dat` are named in `fileNames.dat` but absent from this local v2 `data` directory, so its file-existence gate does not activate them.

The Stout Ar II and Ca IV `.coll` records independently contain the same explicitly attributed 11-value fine-structure Upsilon sequences. Direct numeric parsing found exact equality of the tabulated values in the two local formats; their Ar II temperature coordinates differ slightly by rounding (for example, 1,580 versus 1,585 K). These two local formats make a two-level implementation possible without importing MOCASSIN's interpolation behavior. Preserve one internally consistent energy/A/collision set and check the temperature coordinates and level identities before accepting a translation.

## Reproducibility hashes

All hashes below are SHA-256 of the locally inspected raw files. The triplets are listed in `.nrg`, `.tp`, `.coll` order.

| Source | SHA-256 triplet |
|---|---|
| Stout `ar_1` | `f2e7ee736ee512248777f99b400309365ce7c5d92c198e7b2d80d79f14430cce`; `76f764840c39787d265485f4b367df0ba7da2e3489600ecea99b0a997649fbda`; `0c5de44555d407e01cb9c716105d45309258ea59ed83d300b4280068e993e0d0` |
| Stout `ar_2` | `2f8151e6d5b081488478f3252f37134a72cfafb534a846d205575c1d0571ff30`; `0dc0d90620687c7f67e1b9369b38158d1d1778956796ae335812951d8cd5f655`; `c6f860bc6d8e2f88c5a5e43e487fbf440c56431a593bc9feca22ba0828624e0e` |
| Stout `fe_1` | `aad0a606c0026b8ae8ccfbddebacb0ce966ba4f7543f25389e35d408aaa23f15`; `ad6a2add392db90ce05e1e699eddcf5c1946d721e870eb9467e1c836ee01220a`; `5d3c48ff45822bfb59e6b1dea6b29d11bf403e5e0bcc9fd9e3151e4208533cc9` |
| Stout `ca_3` | `5c4b24667182d73cee0fcf81e8afcdc6ddd1eb234e0830a925db4234b3c2db20`; `0de13d2765627ac1336b2ab030ace55de8c64d2efe6e21eab505c0f02796fcf6`; `4d0fa2cde2fd8ac25b7f3638ea63c7d71b996e8a131f505d9c1cc7a48912eaf3` |
| Stout `ca_4` | `333ea9ead15655b240d573036201aea288c0b4f898d05b96599f816bb7aefaaf`; `fc71397bc3cc9b9225c180a987143fc5ec82bdf1a77d630e0563032dfea77e32`; `e68cce06232964d061e6ef945b353ad7aaaff68ea5522a10aee767b33828587f` |
| MOCASSIN v2 `data/arii.dat` | `c80cfcafdabd9be0b7076ff0bdde749affa711cdbd7495ac8ec260ef7e1d7045` |
| MOCASSIN v2 `data/caiv.dat` | `48819a5cad457867eb08eecea9d6f0c50334a6ca60bc0b6c88e5c17489b44c25` |

## Implementation requirements and unresolved evidence

1. Ar II and Ca IV: retain just the ground fine-structure pair, select one coherent local energy/A/Upsilon set, and represent the 11 physical temperature nodes directly. `TEMPERATURE_TABLE` with `LOG_PCHIP_T` is an implementable interpolation choice already present in MoCHII; it reproduces the table nodes but does not reproduce MOCASSIN's spline of Upsilon against log10(T) between nodes. Quantify that difference before accepting this choice. Hold Upsilon constant only below the first node under the already adopted cold approximation, while the actual T enters Maxwellian rates. Verify source strength, two-level balance, cooling, and line emissivity on 1,000–100,000 K, including each node and first-node neighborhoods. Derive any low-density cooling product from this accepted atom. Higher levels require additional collisional data and a separate convergence judgment.
2. Ar I, Fe I, and Ca III: locate an identified electron-impact source with level-resolved effective collision strengths, radiative rates, source temperature interval, and a physically consistent level/ionization boundary. The present Stout files do not provide such an electron source. Do not set their missing cooling to zero, infer electron rates from H/proton rates, or generate an atom from radiative rows alone. Any discovered source needs direct parser and native profile validation before production use.
3. Source citation: the local Stout and MOCASSIN headers attribute the Ar II/Ca IV strengths to [Pelan and Berrington, *A&AS* **110**, 209 (1995)](https://ui.adsabs.harvard.edu/abs/1995A%26AS..110..209P/abstract), on the chlorine-like ground-state fine-structure transition. A scan of the final journal article is present at `/nfs/mocafe/kiseon/MoCHII/references/Pelan_Berrington_1995_AAS_110_209.pdf`. Its Section 2 describes an LS-coupled R-matrix calculation with 14 terms and an algebraic transformation to fine-structure collision strengths. Section 3 defines Maxwellian effective collision strengths; Table 3 restricts Ar II to `log10(T/K)=3.2–5.2` and Ca IV to `4.0–6.0`, which correspond to the 11 local temperature nodes for each ion after rounding. Its conclusion estimates accuracy better than 10% except for Ar II. These published statements do not validate extrapolation below the first local table node or establish a 10% Ar II error bound. The scan has image pages, so the methods and conclusions were checked against OCR and page images; numerical values still come from the parsed local tables.

This audit did not run a plasma model or install atomic data. MAPPINGS v5, TORUS, CMacIonize, and Wood code were searched for directly usable replacements; no identified, invoked level-resolved electron tables for the three unresolved ions were established by this inspection. Their absence from this audit is not proof of absence from the wider literature.

## Current implementation following this audit

Ar II and Ca IV ground fine-structure atoms and their fits are now installed.
Ar II was accepted only after the cooling registry acquired source-limited
temperature rows and actual rebuilt Intel validation. Its source ceiling is
158,000 K and higher temperature requests are rejected. Original numerical
products were preserved during installation. Source limitation and numerical
validation are recorded in
[the registry report](validation_2026-10-10/source_temperature_registry/README.md).
The current fit inventory is
[35 active metal fits](validation_2026-10-10/production_cooling_coverage_current.json),
not a complete cooling model for all 43 active stages.
