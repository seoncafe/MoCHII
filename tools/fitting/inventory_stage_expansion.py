#!/usr/bin/env python3
"""Read-only source inventory for the proposed sub-100 eV metal stages.

Only the two requested report files are written. The physical-energy reader
does not generate production files on import.
"""

import argparse
import ast
import hashlib
import json
import re
from pathlib import Path
from charge_exchange_expansion import (KF99_ADDED_HEAT, kf96_added_rates,
    kf96_recombination_total_rates, KF99_RECOMBINATION_MEAN_HEAT)
from charge_exchange_capture import langevin_recombination_rates, avoided_crossing_mean_heat


ROOT = Path(__file__).resolve().parents[2]
HOME = Path.home()
GENERATOR = ROOT / "tools/fitting/make_element_data.py"
VERNER = HOME / "RT_Codes/Verner_Fortran_sub"
CHIANTI = HOME / "RT_Codes/CHIANTI/dbase"
CLOUDY = HOME / "CLOUDY/c25.00"
HC_EV_CM = 1.2398419843320026e-4
SHELL_L = (0, 0, 1, 0, 1, 2, 0)


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def source(path):
    return {"path": str(path), "exists": path.is_file(),
            "sha256": digest(path) if path.is_file() else None}


def constants():
    tree = ast.parse(GENERATOR.read_text())
    found = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "ELEMENTS":
                    found[target.id] = ast.literal_eval(node.value)
    if "ELEMENTS" not in found:
        raise ValueError("Generator ELEMENTS declaration is unavailable")
    from ionization_potentials import read_ionization_potentials
    return found["ELEMENTS"], read_ionization_potentials()


def fortran_array(text, label, size):
    match = re.search(r"DATA \(" + label + r"\(I\),I=1," + str(size) +
                      r"\)\s*/([^/]+)/", text)
    if not match:
        raise ValueError(f"Missing Fortran array {label}")
    return [int(x) for x in match.group(1).split(",")]


def verner_tables():
    ph = re.sub(r"\n     \S", " ", (VERNER / "phfit2.f").read_text())
    ci = re.sub(r"\n     \S", " ", (VERNER / "cfit.f").read_text())
    out = {}
    for label, data, n, tail in (("PH1", ph, 6, r",\s*(\d+)"),
                                 ("PH2", ph, 7, ""), ("CF", ci, 5, "")):
        pattern = (r"DATA\s*\(" + label + r"\(I,\s*(\d+),\s*(\d+)" +
                   tail + r"\),I=1," + str(n) + r"\)\s*/([^/]+)/")
        table = {}
        for match in re.finditer(pattern, data):
            *keys, values = match.groups()
            table[tuple(map(int, keys))] = [float(x) for x in values.split(",")]
        out[label] = table
    out["NINN"] = fortran_array(ph, "NINN", 30)
    out["NTOT"] = fortran_array(ph, "NTOT", 30)
    return out


def badnell(path, paired=False):
    sections = {"C": {}, "E": {}}
    section = "C"
    for line in path.read_text().splitlines():
        fields = line.split()
        if paired and len(fields) >= 5 and fields[:2] == ["Z", "N"]:
            section = "C" if fields[4].startswith("C") else "E"
        elif len(fields) >= (5 if paired else 8) and fields[0].isdigit():
            z, ne, m = map(int, fields[:3])
            if m == 1:
                sections[section][(z, ne)] = [float(v) for v in fields[4:]]
    if not paired:
        return sections["C"]
    return {key: {"terms": len(c), "min_excitation_K": min(sections["E"][key])}
            for key, c in sections["C"].items()
            if key in sections["E"] and len(c) == len(sections["E"][key])}


def nist_table():
    result = {}
    for line in (CLOUDY / "data/ionization_potentials.dat").read_text().splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[0].isdigit() and fields[1].isdigit():
            result[(int(fields[0]), int(fields[1]))] = {
                "threshold_eV": float(fields[2]) * HC_EV_CM,
                "uncertainty_eV": float(fields[3]) * HC_EV_CM}
    return result


def numeric_rows(path, style):
    if not path.is_file():
        return None
    count = 0
    for index, line in enumerate(path.read_text(errors="replace").splitlines()):
        fields = line.split()
        if style == "chianti_level":
            valid = bool(fields and fields[0].isdigit())
        elif style == "chianti_radiative":
            valid = len(fields) >= 5 and fields[0].isdigit() and fields[1].isdigit()
        elif style == "chianti_collision":
            valid = (len(fields) == 8 and fields[0].isdigit() and
                     fields[1].isdigit() and fields[5].isdigit() and
                     fields[6].isdigit())
        elif style == "stout_level":
            valid = index > 0 and bool(fields and fields[0].isdigit())
        elif style == "stout_radiative":
            valid = bool(fields and fields[0] == "A")
        else:
            valid = len(fields) >= 3 and fields[:2] == ["CS", "ELECTRON"]
        count += valid
    return count


def stout_collision_metadata(path):
    grid = None
    partners = {}
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split()
        if fields and fields[0] == "TEMP":
            grid = [float(x) for x in fields[1:]]
        elif len(fields) >= 3 and fields[0] in ("CS", "RATE"):
            partner = " ".join(fields[:2])
            partners[partner] = partners.get(partner, 0) + 1
    return {"temperature_grid_K": grid, "collision_partner_row_counts": partners}


def line_sources(element, stage):
    ion = f"{element}_{stage}"
    chi = CHIANTI / element / ion
    stout = CLOUDY / "data/stout" / element / ion
    kinds = {"levels": ("elvlc", "nrg", "chianti_level", "stout_level"),
             "radiative": ("wgfa", "tp", "chianti_radiative", "stout_radiative"),
             "electron_collisions": ("scups", "coll", "chianti_collision", "stout_collision")}
    result = {}
    for label, (cext, sext, cstyle, sstyle) in kinds.items():
        cp = chi / f"{ion}.{cext}"
        candidates = sorted(stout.glob(f"{ion}*.{sext}"))
        candidates_records = []
        for p in candidates:
            item = {**source(p), "rows": numeric_rows(p, sstyle)}
            if label == "electron_collisions":
                item.update(stout_collision_metadata(p))
            candidates_records.append(item)
        result[label] = {
            "chianti": {**source(cp), "rows": numeric_rows(cp, cstyle)},
            "cloudy_stout": candidates_records}
    for label in ("rrparams", "drparams"):
        p = chi / f"{ion}.{label}"
        result[label] = source(p)
        if p.is_file():
            first = p.read_text().splitlines()[0].split()
            result[label]["fit_type"] = int(first[0]) if first else None
    return result


def current_stage_count(element):
    path = ROOT / f"data/atomic/element_{element}.txt"
    match = re.search(r"^ELEMENT\s+\w+\s+\d+\s+(\d+)\s*$", path.read_text(), re.M)
    if not match:
        raise ValueError(f"No ELEMENT declaration: {path}")
    return int(match.group(1)), source(path)


def inventory():
    elements, eth = constants()
    tables = verner_tables()
    rr = badnell(ROOT / "data/atomic/badnell_rr.dat")
    dr = badnell(ROOT / "data/atomic/badnell_dr.dat", paired=True)
    nist = nist_table()
    records = []
    element_records = []
    input_files = [GENERATOR,
                   ROOT / "tools/fitting/charge_exchange_expansion.py",
                   ROOT.parent / "references/Kingdon_Ferland_1999_ApJ_516_L107.pdf",
                   ROOT.parent / "references/Ferland_Korista_Verner_Dalgarno_1997_ApJ_481_L115.pdf",
                   ROOT.parent / "references/Zhang_2026_Fe3plus_to_Fe6plus_supplement.zip",
                   ROOT.parent / "references/Zhang_2026_ApJS_285_50.pdf",
                   ROOT.parent / "references/Fernandez-Menchero_2017JPhysB_50_065203.pdf",
                   ROOT / "tools/fitting/charge_exchange_capture.py",
                   ROOT / "tools/fitting/zhang_fe_recombination.py",
                   CLOUDY / "data/ctrecombdata.dat", VERNER / "phfit2.f", VERNER / "cfit.f",
                   ROOT / "data/atomic/badnell_rr.dat",
                   ROOT / "data/atomic/badnell_dr.dat",
                   CHIANTI / "VERSION", CLOUDY / "data/ionization_potentials.dat",
                   CLOUDY / "source/ion_recomb_Badnell.cpp",
                   CLOUDY / "source/atmdat_dielrec_fe.cpp",
                   CLOUDY / "source/atmdat_adfa.cpp"]
    for element, z, declared in elements:
        current, element_source = current_stage_count(element)
        if current != declared:
            raise ValueError(f"Current file and generator disagree for {element}")
        input_files.append(ROOT / f"data/atomic/element_{element}.txt")
        target = current
        while target <= z:
            ne = z - target + 1
            potential = eth.get((z, ne))
            if potential is None:
                raise ValueError(f'Unidentified physical ionization energy: Z={z}, N={ne}')
            if potential >= 100:
                break
            target += 1
        element_records.append({"element": element, "atomic_number": z,
                                "current_stages": current, "target_stages": target,
                                "added_stages": list(range(current + 1, target + 1)),
                                "current_file": element_source})
        for i in range(current, target):
            ne = z - i + 1
            nrec = z - i
            shells = []
            nout = tables["NTOT"][ne - 1]
            if z == ne and z > 18 or z == ne + 1 and z in (20, 21, 22, 25, 26):
                nout = 7
            for shell in range(1, nout + 1):
                fit = tables["PH1"].get((z, ne, shell))
                if fit is not None and fit[2] > 0:
                    shells.append({"shell_index": shell, "orbital_l": SHELL_L[shell - 1],
                                   "threshold_eV": fit[0]})
            if z in (15, 17, 19) or (z > 20 and z != 26):
                inner = 0.0
            elif ne < 3:
                inner = 1e30
            else:
                inner = tables["PH1"][(z, ne, tables["NINN"][ne - 1])][0]
            cp = tables["CF"].get((z, ne))
            np = nist.get((z, i - 1))
            fit_rr = rr.get((z, nrec))
            fit_dr = dr.get((z, nrec))
            lines = line_sources(element, i + 1)
            for group in ("levels", "radiative", "electron_collisions"):
                item = lines[group]
                input_files.extend([Path(item["chianti"]["path"])])
                input_files.extend(Path(x["path"]) for x in item["cloudy_stout"])
            for group in ("rrparams", "drparams"):
                input_files.append(Path(lines[group]["path"]))
            issues = []
            if (z, ne) not in eth:
                issues.append("Generator ETH has no entry; threshold policy is unresolved")
            if np and abs(min(x["threshold_eV"] for x in shells) - np["threshold_eV"]) > 0.1:
                issues.append("PH1 edge differs from Cloudy NIST ionization potential by more than 0.1 eV")
            if cp and abs(cp[0] - min(x["threshold_eV"] for x in shells)) > 0.1:
                issues.append("Voronov CI fit energy differs from PH1 edge by more than 0.1 eV")
            if fit_rr is None and element != 'fe':
                issues.append("No local Badnell M=1 RR record; evaluate CHIANTI and Cloudy candidates")
            if fit_dr is None and element != 'fe':
                issues.append("No local Badnell M=1 DR record; low-temperature DR is unresolved")
            if lines["electron_collisions"]["chianti"]["rows"] is None and not any(
                    x["rows"] for x in lines["electron_collisions"]["cloudy_stout"]):
                issues.append("No identified electron collision rows for added stage")
            if element == "cl" and i + 1 == 7:
                issues.append("Stout Cl VII collision grid begins at 9800 K; adopted constant first-node Upsilon below that temperature, with actual T in Maxwellian rates; cold resonance accuracy unverified")
            if element == "fe":
                issues.append("Zhang 2026 final-publisher ground-level RR and S1 DR fits are verified and selected; compact DR differs from author ADF09 at some nodes by more than the stated 2%; linked calculations passed 1204 rate checks")
            cx_key=(element,i)
            cx_rate={**kf96_added_rates(), **kf96_recombination_total_rates()}.get(cx_key)
            cx_heat={**KF99_ADDED_HEAT, **KF99_RECOMBINATION_MEAN_HEAT}.get(cx_key)
            capture=langevin_recombination_rates().get(cx_key)
            if capture is not None:
                cx_heat=avoided_crossing_mean_heat(i)
                issues.append("User adopted Ferland 1997 unit-probability Langevin capture upper-limit approximation with reduced collision mass and KF99 avoided crossing with exit charge, mean heat 2.86*(q-1) eV; no certified temperature interval, product branching or reverse rate")
            elif cx_rate is None:
                issues.append("Higher-charge transfer remains unselected: the verified Ferland 1997 Langevin prescription is a unit-probability capture upper limit, not an independently established species-specific reaction rate")
            else:
                issues.append("KF96 total recombination with KF99 mean heat is available for selective staging; ground-state reactants and low-density radiative relaxation are assumed, product branching and diffuse cascade reabsorption are not modeled")
            records.append({"element": element, "atomic_number": z,
                            "absorbing_stage": i, "added_stage": i + 1,
                            "absorbing_charge": i - 1, "absorbing_electrons": ne,
                            "recombining_electrons": nrec,
                            "thresholds_eV": {"ph1_lowest_occupied": min(x["threshold_eV"] for x in shells),
                                              "nist_cloudy": np, "generator_eth": eth.get((z, ne)),
                                              "ci_fit_parameter": cp[0] if cp else None},
                            "photoionization": {"ph1_occupied_shells": shells,
                                                "ph2_present": (z, ne) in tables["PH2"],
                                                "ph2_required_below_inner_edge": inner > 0 and inner < 1e30,
                                                "ph1_switch_eV": inner},
                            "collisional_ionization": {"voronov_fit_present": cp is not None,
                                                       "dere_hybrid_factor_review_required": True},
                            "recombination": {"badnell_rr_m1_present": fit_rr is not None,
                                              "badnell_dr_m1": fit_dr,
                                              "chianti_rrparams": lines["rrparams"],
                                              "chianti_drparams": lines["drparams"]},
                            "charge_exchange": {
                                "with_H_plus": "excited-product reverse not modeled" if cx_rate is not None or capture is not None else "unverified",
                                "with_H_neutral": "User-adopted Langevin capture upper limit" if capture is not None else "unverified" if cx_rate is None else "KF96 published total rate",
                                "capture_rate_cm3_s":capture,
                                "heat_defects_eV": {"recombination":cx_heat} if cx_heat is not None else None,
                                "heat_source": "KF99 avoided crossing with exit charge, 2.86*(q-1) eV" if capture is not None else None if cx_rate is None else "KF99 Table 1, DOI 10.1086/312008",
                                "temperature_interval_K": None if cx_rate is None else [cx_rate.min_temperature,cx_rate.max_temperature],
                                "product_branching_verified": False},
                            "line_sources": {k: lines[k] for k in ("levels", "radiative", "electron_collisions")},
                            "selected_level_count": None,
                            "selected_radiative_transition_count": None,
                            "selected_collision_transition_count": None,
                            "unresolved": issues})
    return {"scope": "metal stages reached by successive physical ground-state ionization potentials below 100 eV",
            "units": {"energy": "eV", "temperature": "K"},
            "database_versions": {"chianti": (CHIANTI / "VERSION").read_text().strip(),
                                  "cloudy": "c25.00 (local directory label)",
                                  "badnell_rr_header": (ROOT / "data/atomic/badnell_rr.dat").open().readline().strip(),
                                  "badnell_dr_header": (ROOT / "data/atomic/badnell_dr.dat").open().readline().strip()},
            "sources": [source(p) for p in sorted(set(input_files))],
            "elements": element_records, "transitions": records,
            "totals": {"current_metal_stages": sum(x["current_stages"] for x in element_records),
                       "target_metal_stages": sum(x["target_stages"] for x in element_records),
                       "added_stages": len(records),
                       "target_adjacent_transitions": sum(x["target_stages"] - 1 for x in element_records)}}


def markdown(data):
    lines = ["# Atomic stage expansion source preflight", "",
             "Generated by `python3 tools/fitting/inventory_stage_expansion.py`. This is a source inventory, not generated atomic data or a validated physical model.", "",
             "The target criterion uses the physical ionization potential from Cloudy's identified NIST ASD snapshot dated 2022-11-15. The occupied-valence edge and ETH use that potential; inner PH1 edges and Voronov CI fit parameters retain their independent source values. Original PH1 edges remain in this audit for comparison.", "",
             f"Measured totals: {data['totals']['current_metal_stages']} current metal stages; {data['totals']['added_stages']} proposed additions; {data['totals']['target_metal_stages']} target stages; {data['totals']['target_adjacent_transitions']} target adjacent transitions.", "",
             "| Element | Current | Target | Added |", "|---|---:|---:|---:|"]
    for el in data["elements"]:
        lines.append(f"| {el['element'].capitalize()} | {el['current_stages']} | {el['target_stages']} | {len(el['added_stages'])} |")
    lines += ["", "## Proposed transitions", "",
              "PH1/PH2 and CI indicate source fit presence. RR/DR mean local Badnell M=1 records. CHIANTI and Stout counts are parsed source rows, not selected model sizes or validated data. `?` means unresolved.", "",
              "| Transition | PH1/NIST/CI [eV] | Shells | PH2 | RR/DR | CHIANTI levels/A/CS | Stout levels/A/CS |",
              "|---|---|---:|---|---|---|---|"]
    for rec in data["transitions"]:
        t = rec["thresholds_eV"]
        ph = rec["photoionization"]
        ls = rec["line_sources"]
        c = [ls[k]["chianti"]["rows"] for k in ("levels", "radiative", "electron_collisions")]
        s = [max((x["rows"] or 0 for x in ls[k]["cloudy_stout"]), default=0)
             for k in ("levels", "radiative", "electron_collisions")]
        ctext = "/".join("-" if x is None else str(x) for x in c)
        stext = "/".join(str(x) for x in s)
        nist = t["nist_cloudy"]["threshold_eV"] if t["nist_cloudy"] else float("nan")
        lines.append(f"| {rec['element'].capitalize()} {rec['absorbing_stage']} -> {rec['added_stage']} | {t['ph1_lowest_occupied']:.3f}/{nist:.3f}/{t['ci_fit_parameter']:.3f} | {len(ph['ph1_occupied_shells'])} | {'yes' if ph['ph2_present'] else 'no'} | {'yes' if rec['recombination']['badnell_rr_m1_present'] else '?'} / {'yes' if rec['recombination']['badnell_dr_m1'] else '?'} | {ctext} | {stext} |")
    lines += ["", "## Unresolved inputs and interpretation", "",
              "Eleven added boundaries have KF96 total recombination fits and KF99 mean heat under ground-state, low-density assumptions. Fe4 is a Landau-Zener estimate with a radiative-transfer caveat, and Fe V can now be generated with verified final-journal RR/DR. All Fe target stages now have selected recombination and capture inputs. N4 uses the original KF96 fit instead of Cloudy's unidentified revision. The user adopted the Langevin unit-transfer capture upper limit for the ten higher-charge boundaries, with reduced collision mass and KF99 avoided crossing with exit charge, mean heat 2.86*(q-1) eV. Their accuracy and product branching are not established. Existing Cl3/4 and Ca3/4 omissions have been corrected without changing stage counts. No missing rate is replaced by an inferred zero. Product branching and diffuse cascade reabsorption remain unverified. Selected level and transition counts remain unset until convergence tests; raw source row counts do not establish a usable atom.", ""]
    for rec in data["transitions"]:
        label = f"{rec['element'].capitalize()} {rec['absorbing_stage']} -> {rec['added_stage']}"
        special = [x for x in rec["unresolved"] if not x.startswith("Added-boundary")]
        if special:
            lines.append(f"- {label}: " + "; ".join(special) + ".")
    lines += ["", "## Source identification", "",
              f"CHIANTI version: `{data['database_versions']['chianti']}`. Cloudy directory: `c25.00`. Badnell table dates are recorded in the JSON headers.", "",
              "| Source | SHA-256 |", "|---|---|"]
    keys = ("phfit2.f", "cfit.f", "badnell_rr.dat", "badnell_dr.dat",
            "ionization_potentials.dat", "VERSION", "ion_recomb_Badnell.cpp",
            "charge_exchange_expansion.py", "ctrecombdata.dat",
            "Kingdon_Ferland_1999_ApJ_516_L107.pdf",
            "Ferland_Korista_Verner_Dalgarno_1997_ApJ_481_L115.pdf",
            "Zhang_2026_Fe3plus_to_Fe6plus_supplement.zip",
            "Zhang_2026_ApJS_285_50.pdf", "Fernandez-Menchero_2017JPhysB_50_065203.pdf",
            "charge_exchange_capture.py", "zhang_fe_recombination.py")
    for item in data["sources"]:
        if Path(item["path"]).name in keys:
            lines.append(f"| `{item['path']}` | `{item['sha256']}` |")
    lines += ["", "The Stout table shows the largest count among local alternatives for each file type; those maxima do not define one atom. The complete machine-readable record includes each occupied shell edge, source paths and SHA-256 hashes, source fit types, all measured line row counts, and explicit nulls for unselected atom sizes. Cloudy `c25.00` is a local directory label; this check does not independently verify a release tag or published provenance. CHIANTI and Stout electron collision source counts do not establish matching collision partners, temperature validity, or fit accuracy. These require line-by-line checks before production generation.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "md/atomic_stage_expansion")
    args = parser.parse_args()
    data = inventory()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "source_inventory.json").write_text(json.dumps(data, indent=2) + "\n")
    (args.output_dir / "source_inventory.md").write_text(markdown(data))
    print(json.dumps(data["totals"], sort_keys=True))


if __name__ == "__main__":
    main()
