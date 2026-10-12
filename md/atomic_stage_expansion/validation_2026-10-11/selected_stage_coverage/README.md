# Current selected atom and active-stage coverage

`tools/fitting/audit_selected_stage_coverage.py` measures the current selected
model, verifies its manifest hashes, follows the original MOCASSIN v2 file
mapping, counts installed levels/radiation/collision rows, checks positive
radiation energy gaps, and reads the production ELEMENT/TRANSITION declarations.
It does not certify atomic accuracy or numerical consumer validation.

| Selection | Installed atoms | Explicit levels | Missing supplied atoms | Active metal stages | Active adjacent transitions |
| --- | ---: | ---: | --- | ---: | ---: |
| mocassin_v2 | 55 | 571 | Ca VII | 43 | 32 |
| extended after Mg I correction | 62 | 1951 | Ca VII | 43 | 32 |

`current_reference.json` and `extended_after.json` are current measurements.
`extended_before.json` preserves the confirmed Mg I omission. All installed
supplied models have at least the reference level count, positive radiation
row count, and positive collision pair count. No checked installed radiation
row has a nonpositive photon energy. These structural conditions do not prove
physical completeness: configuration identity, collision strength validity,
and the Ca VII correction still require physical review.

Production activation is incomplete: the target is 77 metal stages and 66
adjacent transitions. Prepared rate files and completed staged validation do
not substitute for the declarations read by production species setup.
