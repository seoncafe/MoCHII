#!/usr/bin/env python3
"""Source-based terminal ionization diagnostic for the approved 64-stage roster.

This is an electron-only, zero-radiation stationary calculation. The extra
transition is evaluated from local atomic source fits, not loaded by MoCHII.
No unsupported terminal charge-transfer rate is set to zero in a physical
prediction: hydrogen is fully ionized and charge transfer is excluded from
both compared cascades by construction.
"""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import make_element_data as sources


ROSTER = {"c": (6, 5), "n": (7, 6), "o": (8, 5), "ne": (10, 5),
          "s": (16, 7), "ar": (18, 7), "mg": (12, 4), "fe": (26, 7),
          "si": (14, 5), "cl": (17, 7), "ca": (20, 6)}
KB_EV = 8.617333262145e-5


def read_serialized(directory, element):
    path = directory / f"element_{element}.txt"
    rows = {}
    stage = None
    for line in path.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if fields[0] == "TRANSITION":
            stage = int(fields[1])
            rows[stage] = {}
        elif fields[0] in {"ETH", "CI", "RR", "RR2", "DR", "DR2"}:
            rows[stage][fields[0]] = [float(x) for x in fields[1:]]
    return rows


def hybrid_factors():
    source = Path("src/recomb_mod.f90").read_text()
    section = source.split("function ci_dere_ratio", 1)[1]
    section = section.split("end function ci_dere_ratio", 1)[0]
    factors = {}
    for z, values in re.findall(r"case \((\d+)\); tab\(1:\d+\) = \[([^\]]+)\]", section):
        factors[int(z)] = [float(x.replace("_wp", "")) for x in values.split(",")]
    if len(factors) < 11:
        raise ValueError("Could not parse runtime HYBRID factors")
    return factors


def recombination(row, temperature):
    if "RR2" in row:
        a, eta = row["RR2"]
        rate = a * (temperature / 1e4) ** (-eta)
    else:
        a, b, t0, t1, c, t2 = row["RR"]
        s0, s1 = np.sqrt(temperature / t0), np.sqrt(temperature / t1)
        bt = b + c * np.exp(-t2 / temperature)
        rate = a / (s0 * (1 + s0) ** (1 - bt) * (1 + s1) ** (1 + bt))
    if "DR" in row:
        fields = row["DR"]
        count = int(fields[0])
        for c, energy in zip(fields[1:1 + count], fields[1 + count:1 + 2 * count]):
            rate += temperature ** -1.5 * c * np.exp(-energy / temperature)
    if "DR2" in row:
        a, b, t0, t1 = row["DR2"]
        rate += a * temperature ** -1.5 * np.exp(-t0 / temperature) * (1 + b * np.exp(-t1 / temperature))
    return rate


def ionization_log(row, temperature):
    energy, p, a, x, k = row["CI"]
    u = energy / (KB_EV * temperature)
    return np.log(a) + np.log1p(p * np.sqrt(u)) + k * np.log(u) - u - np.log(x + u)


def source_terminal(element, charge, stage, rr_table, dr_table, ci_table):
    ne = charge - stage + 1
    key = (charge, ne)
    row = {"ETH": [sources.ETH[key]], "CI": ci_table[key]}
    rec_key = (charge, charge - stage)
    rr_origin = "Badnell M=1"
    dr_origin = "Badnell M=1"
    if rec_key in rr_table:
        row["RR"] = rr_table[rec_key]
    else:
        kind, values = sources.read_rrparams(element, stage + 1)
        row["RR" if kind == "bad" else "RR2"] = values
        rr_origin = "CHIANTI"
    if rec_key in dr_table:
        c, energy = dr_table[rec_key]
        row["DR"] = [len(c), *c, *energy]
    else:
        dr = sources.read_drparams(element, stage + 1)
        if dr is None:
            raise ValueError(f"Missing terminal DR source for {element}")
        if dr[0] == "t1":
            row["DR"] = [len(dr[1]), *dr[1], *dr[2]]
        else:
            row["DR2"] = dr[1]
        dr_origin = "CHIANTI"
    return row, rr_origin, dr_origin


def log_fractions(log_ion, rec):
    log_pop = np.zeros(len(log_ion) + 1)
    for i, (log_ci, alpha) in enumerate(zip(log_ion, rec)):
        if alpha <= 0:
            raise ValueError("Nonpositive recombination rate in sampled domain")
        log_pop[i + 1] = log_pop[i] + log_ci - math.log(alpha)
    peak = max(log_pop)
    return log_pop - (peak + math.log(sum(np.exp(log_pop - peak))))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--atomic-dir", type=Path,
                        default=Path("md/atomic_stage_expansion/validation_2026-10-10/"
                                     "langevin64_exit_heat/products/atomic"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    temps = np.geomspace(1e3, 1e5, 801)
    rr_table, dr_table, ci_table = (sources.parse_badnell_rr(),
                                      sources.parse_badnell_dr(), sources.parse_cf())
    factors = hybrid_factors()
    results = []
    full_profiles = {"temperature_K": temps}
    fig, axes = plt.subplots(4, 3, figsize=(12, 13), sharex=True)
    efig, eaxes = plt.subplots(4, 3, figsize=(12, 13), sharex=True)
    for axis, (element, (charge, terminal)) in zip(axes.flat, ROSTER.items()):
        rows = read_serialized(args.atomic_dir, element)
        if len(rows) != terminal - 1:
            raise ValueError(f"Unexpected serialized stage count: {element}")
        extra, rr_origin, dr_origin = source_terminal(element, charge, terminal,
                                                       rr_table, dr_table, ci_table)
        if (element, terminal) in sources.CX_HIGH:
            raise ValueError(f"Terminal charge transfer unexpectedly present: {element}")
        native = np.zeros((len(temps), terminal))
        extended = np.zeros((len(temps), terminal + 1))
        native_log = np.zeros_like(native)
        extended_log = np.zeros_like(extended)
        ratio = np.zeros(len(temps))
        extra_energy_fraction = np.zeros(len(temps))
        factor_row = factors[charge]
        if len(factor_row) != terminal - 1:
            raise ValueError(f"Unexpected HYBRID factor count: {element}")
        for j, temperature in enumerate(temps):
            log_ci = [ionization_log(rows[i], temperature) + math.log(factor_row[i - 1])
                  for i in range(1, terminal)]
            ci = [math.exp(x) for x in log_ci]
            alpha = [recombination(rows[i], temperature) for i in range(1, terminal)]
            native_log[j] = log_fractions(log_ci, alpha)
            native[j] = np.exp(native_log[j])
            terminal_log_ci = ionization_log(extra, temperature)
            terminal_ci = math.exp(terminal_log_ci)
            terminal_alpha = recombination(extra, temperature)
            ratio[j] = terminal_ci / terminal_alpha
            extended_log[j] = log_fractions(log_ci + [terminal_log_ci], alpha + [terminal_alpha])
            extended[j] = np.exp(extended_log[j])
            energy_extended = sum(extended[j, i - 1] * ci[i - 1] * rows[i]["ETH"][0]
                                  for i in range(1, terminal))
            added = extended[j, terminal - 1] * terminal_ci * extra["ETH"][0]
            extra_energy_fraction[j] = added / (energy_extended + added) if energy_extended + added > 0 else 0
        full_profiles[f"{element}_truncated_fraction"] = native
        full_profiles[f"{element}_extended_fraction"] = extended
        full_profiles[f"{element}_truncated_log_fraction"] = native_log
        full_profiles[f"{element}_extended_log_fraction"] = extended_log
        full_profiles[f"{element}_added_ci_energy_fraction"] = extra_energy_fraction
        for stage in range(terminal - 1):
            axis.semilogx(temps, extended_log[:, stage] / math.log(10),
                          color="0.75", linewidth=0.7)
        axis.semilogx(temps, native_log[:, -1] / math.log(10), color="black", label="truncated terminal")
        axis.semilogx(temps, extended_log[:, -2] / math.log(10),
                      color="tab:blue", linestyle="--", label="extended terminal")
        axis.semilogx(temps, extended_log[:, -1] / math.log(10), color="tab:red", label="next stage")
        axis.set_title(element.upper())
        axis.set_ylim(-100, 0.05)
        axis.set_ylabel("Log10 ion fraction")
        axis.set_xlabel("Temperature [K]")
        energy_axis = eaxes.flat[len(results)]
        energy_axis.loglog(temps, np.where(extra_energy_fraction > 0,
                                          extra_energy_fraction, np.nan))
        energy_axis.set_title(element.upper())
        energy_axis.set_ylabel("Added CI energy fraction")
        energy_axis.set_xlabel("Temperature [K]")
        energy_axis.set_ylim(1e-35, 1)
        result = {
            "element": element, "terminal_stage": terminal,
            "next_stage": terminal + 1, "edge_eV": extra["ETH"][0],
            "rr_source": rr_origin, "dr_source": dr_origin,
            "terminal_hybrid_factor": 1.0,
            "terminal_hybrid_factor_status": "No tabulated factor; runtime function returns unity",
            "terminal_charge_transfer_status": "No validated rate in CX_HIGH",
            "max_next_stage_fraction": float(max(extended[:, -1])),
            "max_fraction_difference": float(np.max(np.abs(native - extended[:, :-1]))),
            "max_added_ci_energy_fraction": float(max(extra_energy_fraction)),
            "ratio_at_100000_K": float(ratio[-1]),
            "next_fraction_at_100000_K": float(extended[-1, -1]),
        }
        results.append(result)
    for axis in axes.flat[len(ROSTER):]:
        axis.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    axes.flat[-1].legend(handles, labels, loc="center")
    fig.tight_layout()
    fig.savefig(args.output / "ion_fractions.pdf")
    plt.close(fig)
    for axis in eaxes.flat[len(ROSTER):]:
        axis.axis("off")
    efig.tight_layout()
    efig.savefig(args.output / "ionization_energy.pdf")
    plt.close(efig)
    np.savez_compressed(args.output / "full_profiles.npz", **full_profiles)
    result = {
        "case": "stationary electron-only adjacent stages, no photons, fully ionized hydrogen, no charge transfer",
        "density": "density cancels from electron-only stationary fractions; energy result is fraction of CI ionization-energy loss",
        "temperature_K": [float(temps[0]), float(temps[-1])],
        "sample_count": len(temps), "comparison": results,
        "atomic_dir": str(args.atomic_dir),
        "source_sha256": {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [Path(__file__), Path("src/recomb_mod.f90"),
                         *sorted(args.atomic_dir.glob("element_*.txt"))]
        },
    }
    (args.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({r["element"]: {"fraction": r["max_next_stage_fraction"],
                                     "energy": r["max_added_ci_energy_fraction"]}
                      for r in results}, indent=2))


if __name__ == "__main__":
    main()
