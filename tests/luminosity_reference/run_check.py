#!/usr/bin/env python3
"""Check full-spectrum source normalization through the production ion_setup path."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "tests" / "ion_stages"))
from run_checks import link

KB_EV = 8.617333262145e-5
EV2ERG = 1.602176634e-12
H_CGS = 6.62607015e-27
HC_ANG = 12398.42
EION = 13.598
EFUV = 6.0
EMAX = 100.0


def near(actual, expected, label, rtol=3e-8, atol=1e-10):
    if not math.isclose(actual, expected, rel_tol=rtol, abs_tol=atol):
        raise AssertionError(f"{label}: got {actual:.16g}, expected {expected:.16g}")


def table(path, energy, columns, units="per_ev"):
    rows = []
    for i, e in enumerate(energy):
        if units in ("per_ev", "shape_ev"):
            x, factor = e, 1.0
        elif units == "per_hz":
            x, factor = e * EV2ERG / H_CGS, H_CGS / EV2ERG
        elif units == "per_ang":
            x, factor = HC_ANG / e, e * e / HC_ANG
        elif units == "per_um":
            x, factor = HC_ANG / (1e4 * e), e * e * 1e4 / HC_ANG
        else:
            raise ValueError(units)
        rows.append((x, *(column[i] * factor for column in columns)))
    if units in ("per_ang", "per_um"):
        rows.reverse()
    path.write_text("".join(" ".join(f"{v:.17e}" for v in row) + "\n" for row in rows))


def linear_integral(energy, density, low, high):
    total = 0.0
    for e0, e1, f0, f1 in zip(energy[:-1], energy[1:], density[:-1], density[1:]):
        lo, hi = max(e0, low), min(e1, high)
        if hi <= lo:
            continue
        slope = (f1 - f0) / (e1 - e0)
        total += f0 * (hi - lo) + 0.5 * slope * ((hi - e0) ** 2 - (lo - e0) ** 2)
    return total


def planck_integral(low, high, temperature):
    kt = KB_EV * temperature
    value, error = quad(lambda x: x**3 / math.expm1(x), low / kt, high / kt,
                        epsabs=1e-12, epsrel=1e-12)
    return value * kt**4


def config(path, **params):
    defaults = dict(nsource=1, nnu_ion=24, eion_min=EION, eion_max=EMAX,
                    nphotons=120, nphotons_fuv=30, ion_align_edges=False,
                    luminosity_reference="full_spectrum")
    defaults.update(params)
    def literal(value):
        if isinstance(value, bool):
            return ".true." if value else ".false."
        if isinstance(value, str):
            return "'" + value.replace("'", "''") + "'"
        if isinstance(value, (list, tuple)):
            return ", ".join(literal(x) for x in value)
        return str(value)
    path.write_text("&parameters\n" + "".join(f" par%{k} = {literal(v)}\n" for k, v in defaults.items()
                                              if v is not None) + "/\n")


def invoke(exe, directory, name, **params):
    input_file = directory / (name + ".in")
    config(input_file, **params)
    command = ["mpirun", "-np", "2", str(exe), str(input_file)]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                            env=dict(os.environ, OMP_NUM_THREADS="1"), timeout=120)
    (directory / (name + ".log")).write_text(result.stdout + result.stderr)
    if result.returncode != 0:
        raise RuntimeError(f"{name}: exit {result.returncode}: {(result.stdout + result.stderr)[-2000:]}")
    values = {}
    for line in result.stdout.splitlines():
        if not line.startswith("RESULT "):
            continue
        fields = line.split()
        if fields[1] == "source":
            values.setdefault("source", {})[int(fields[2])] = list(map(float, fields[3:]))
        elif fields[1] in ("packet_counts", "source_count"):
            values[fields[1]] = [int(x) for x in fields[2:]] if fields[1] == "packet_counts" else int(fields[2])
        else:
            data = list(map(float, fields[2:]))
            values[fields[1]] = data if len(data) > 1 else data[0]
    if "Ltot" not in values:
        raise AssertionError(f"{name}: no driver result")
    return values


def check_source(result, index, energy, density, scale, name, fuv=False):
    row = result["source"][index]
    full = linear_integral(energy, density, energy[0], energy[-1])
    lion = scale * linear_integral(energy, density, EION, EMAX) / full
    lfuv = scale * linear_integral(energy, density, EFUV, EION) / full if fuv else 0.0
    near(row[0], energy[0], name + " Emin")
    near(row[1], energy[-1], name + " Emax")
    near(row[2], scale, name + " Lfull")
    near(row[3], lion, name + " Lion")
    near(row[4], lfuv, name + " Lfuv")
    near(row[5], scale - lion - lfuv, name + " Loutside")
    near(row[6], lion + lfuv, name + " Ltransport")
    return lion + lfuv


def check_packets(result, name):
    total, ionizing, fuv = result["packet_luminosities"]
    nwhole, nion, nfuv = result["packet_counts"]
    near(total * nwhole, result["Ltot"], name + " whole packets")
    if nion:
        near(ionizing * nion, result["Lion"], name + " ionizing packets")
    if nfuv:
        near(fuv * nfuv, result["Lfuv"], name + " FUV packets")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path)
    args = parser.parse_args()
    context = tempfile.TemporaryDirectory(prefix="luminosity_reference_") if args.work_dir is None else None
    work = args.work_dir or Path(context.name)
    work.mkdir(parents=True, exist_ok=True)
    report = {"passed": False, "checks": [], "work_dir": str(work)}
    try:
        exe = link(work, "tests/luminosity_reference/check_luminosity_reference.f90")
        energy = [1.0, 6.0, EION, 20.0, 60.0, 100.0, 200.0]
        flat = [2.0] * len(energy)
        power = [0.1 * e for e in energy]
        for label, density in (("flat", flat), ("power_law", power)):
            for units in ("per_ev", "per_hz", "per_ang", "per_um"):
                name = label + "_" + units
                sed = work / (name + ".sed")
                table(sed, energy, [density], units)
                result = invoke(exe, work, name, ion_spectrum=str(sed), spectrum_type=units,
                                luminosity=1e39, add_fuv=True, efuv_min=EFUV,
                                fuv_sampling="separate")
                expected = check_source(result, 1, energy, density, 1e39, name, True)
                near(result["Ltot"], expected, name + " total")
                near(result["binned"], expected, name + " diagnostic bins")
                near(result["Lion"] + result["Lfuv"], expected, name + " split")
                check_packets(result, name)
                report["checks"].append(name)

        sed = work / "flat_per_ev.sed"
        implicit = invoke(exe, work, "default_reference", ion_spectrum=str(sed),
                          spectrum_type="per_ev", luminosity=1e39,
                          luminosity_reference=None)
        explicit = invoke(exe, work, "explicit_ionizing", ion_spectrum=str(sed),
                          spectrum_type="per_ev", luminosity=1e39,
                          luminosity_reference="ionizing_band")
        for key in ("Ltot", "Lion", "Lfuv", "binned", "packet_luminosities"):
            a, b = implicit[key], explicit[key]
            if a != b:
                raise AssertionError(f"default and explicit ionizing differ for {key}: {a}, {b}")
        if implicit["source_count"] or explicit["source_count"]:
            raise AssertionError("source diagnostics should exist only in full_spectrum mode")
        report["checks"].append("default_equals_explicit_ionizing")

        unset = invoke(exe, work, "absolute_unset", ion_spectrum=str(sed),
                       spectrum_type="per_ev", luminosity=-999.0)
        integral = linear_integral(energy, flat, EION, EMAX)
        near(unset["Ltot"], integral, "absolute unset")
        near(unset["parameter_luminosity"],
             linear_integral(energy, flat, energy[0], energy[-1]),
             "absolute unset recorded full luminosity")
        check_source(unset, 1, energy, flat, linear_integral(energy, flat, energy[0], energy[-1]),
                     "absolute unset")
        check_packets(unset, "absolute unset")
        report["checks"].append("absolute_unset")

        narrow = invoke(exe, work, "narrow_transport", ion_spectrum=str(sed),
                        spectrum_type="per_ev", luminosity=1e39, eion_max=80.0)
        near(narrow["source"][1][2], 1e39, "narrow source full integral")
        near(narrow["Ltot"], 1e39 * linear_integral(energy, flat, EION, 80.0)
             / linear_integral(energy, flat, energy[0], energy[-1]), "narrow transport")
        near(narrow["Ltot"] / (80.0 - EION),
             unset["Ltot"] / (EMAX - EION) * 1e39
             / linear_integral(energy, flat, energy[0], energy[-1]),
             "transport interval spectral density")
        report["checks"].append("transport_band_preserves_density")

        planck = invoke(exe, work, "planck", tstar=30000.0, luminosity=1e39,
                        add_fuv=True, efuv_min=EFUV, fuv_sampling="separate")
        full = math.pi**4 / 15 * (KB_EV * 30000.0)**4
        lion = 1e39 * planck_integral(EION, EMAX, 30000.0) / full
        lfuv = 1e39 * planck_integral(EFUV, EION, 30000.0) / full
        near(planck["Ltot"], lion + lfuv, "Planck transport", rtol=2e-7)
        near(planck["source"][1][2], 1e39, "Planck Lfull")
        near(planck["source"][1][3], lion, "Planck Lion", rtol=2e-7)
        near(planck["source"][1][4], lfuv, "Planck Lfuv", rtol=2e-7)
        near(planck["source"][1][1], -1.0, "Planck infinity sentinel")
        check_packets(planck, "Planck")
        report["checks"].append("planck_full_integral")

        multi_sed = work / "multi.sed"
        table(multi_sed, energy, [flat, power])
        multi = invoke(exe, work, "multi_precedence", nsource=2, src_spectrum_file=str(multi_sed),
                       ion_spectrum=str(sed), src_tstar=[20000.0, 80000.0],
                       src_lum=[2e38, 3e38], spectrum_type="per_ev")
        expected = (check_source(multi, 1, energy, flat, 2e38, "multi 1") +
                    check_source(multi, 2, energy, power, 3e38, "multi 2"))
        near(multi["Ltot"], expected, "multi total")
        check_packets(multi, "multi")
        report["checks"].append("multi_source_precedence")

        mixed_spectra = invoke(exe, work, "temperature_over_global_file",
                               nsource=2, ion_spectrum=str(sed),
                               src_tstar=[30000.0, -999.0], src_lum=[2e38, 3e38],
                               spectrum_type="per_ev")
        first_planck = 2e38 * planck_integral(EION, EMAX, 30000.0) / (
            math.pi**4 / 15 * (KB_EV * 30000.0)**4)
        second_table = check_source(mixed_spectra, 2, energy, flat, 3e38,
                                    "global file source")
        near(mixed_spectra["source"][1][0], 0.0, "source Planck Emin")
        near(mixed_spectra["source"][1][1], -1.0, "source Planck infinity")
        near(mixed_spectra["source"][1][3], first_planck,
             "source temperature over global file", rtol=2e-7)
        near(mixed_spectra["Ltot"], first_planck + second_table,
             "mixed spectrum precedence", rtol=2e-7)
        report["checks"].append("component_temperature_over_global_file")

        external = dict(ion_spectrum=str(sed), spectrum_type="per_ev",
                        luminosity=1e39, src_lum=1e39, ext_tstar=40000.0,
                        ext_intensity=2e37, distance2cm=1.0)
        mixed_band = invoke(exe, work, "mixed_ionizing", **external,
                            luminosity_reference="ionizing_band")
        mixed_full = invoke(exe, work, "mixed_full", **external)
        external_luminosity = math.pi * 6.0 * external["ext_intensity"]
        internal_full = 1e39 * linear_integral(energy, flat, EION, EMAX) / linear_integral(
            energy, flat, energy[0], energy[-1])
        near(mixed_band["Ltot"], 1e39 + external_luminosity,
             "mixed ionizing total")
        near(mixed_full["Ltot"], internal_full + external_luminosity,
             "mixed full total")
        near(mixed_band["Ltot"] - 1e39,
             mixed_full["Ltot"] - mixed_full["source"][1][6],
             "external luminosity invariant")
        check_packets(mixed_band, "mixed ionizing")
        check_packets(mixed_full, "mixed full")
        report["checks"].append("mixed_internal_external")

        fuv_energy = [1.0, 6.0, 8.0, 11.0, 13.0]
        fuv_density = [1.0] * len(fuv_energy)
        fuv_sed = work / "fuv_only.sed"
        table(fuv_sed, fuv_energy, [fuv_density])
        fuv_only = invoke(exe, work, "fuv_only", ion_spectrum=str(fuv_sed),
                          spectrum_type="per_ev", luminosity=1e39,
                          add_fuv=True, efuv_min=EFUV, fuv_sampling="separate")
        near(fuv_only["Ltot"], check_source(fuv_only, 1, fuv_energy, fuv_density,
                                            1e39, "FUV only", True), "FUV total")
        if fuv_only["packet_counts"][1] != 0 or fuv_only["packet_counts"][2] != 30:
            raise AssertionError("FUV-only source launched wrong strata")
        check_packets(fuv_only, "FUV only")
        report["checks"].append("fuv_only_overlap")
        report["passed"] = True
    finally:
        (work / "result.json").write_text(json.dumps(report, indent=2) + "\n")
        if context is not None:
            print(json.dumps(report, indent=2))
            context.cleanup()
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
