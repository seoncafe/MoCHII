"""Compare private oxygen atoms with independent collision/radiation balance."""

import json
import hashlib
import subprocess
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools" / "fitting"))
import chianti_cooling as C  # noqa: E402
from validate_excitation_density_profiles import full_excitation_balance  # noqa: E402

HERE = Path(__file__).parent
STATE = HERE / "controller_state_v2" / "state.json"
DEADLINE = time.monotonic() + 7800
while True:
    if STATE.exists():
        state = json.loads(STATE.read_text())
        if state["status"] == "completed":
            break
        if state["status"] == "failed":
            raise RuntimeError("Native predecessor failed; preserve its outputs")
    if time.monotonic() >= DEADLINE:
        raise TimeoutError("Native predecessor did not complete within 7800 s")
    time.sleep(30)

native_result = json.loads((HERE / "native_result.json").read_text())
manifest = json.loads((HERE / "candidate_manifest.json").read_text())
if not native_result.get("passed") or not manifest.get("passed"):
    raise ValueError("Incomplete predecessor results")
TEMPERATURE = np.geomspace(1000.0, 100000.0, 401)
DENSITY = (1e-6, 1e2, 1e6, 1e10)
saved = {"temperature_K": TEMPERATURE, "density_cm-3": np.asarray(DENSITY)}
report = {}
fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)

def atom_radiation(path):
    lines = path.read_text().splitlines()
    mark = next(i for i, line in enumerate(lines) if line.startswith("NRAD "))
    nrad = int(lines[mark].split()[1])
    parsed = []
    for line in lines[mark + 1:mark + 1 + nrad]:
        lower, upper, rate = line.split()
        parsed.append((int(lower), int(upper), 0.0, 0.0, float(rate)))
    return parsed

for axis, stage in zip(axes, (2, 3)):
    levels, _, collisions = C.bound_excitation_atom("o", stage)
    atom = Path(manifest["ions"][f"O{stage}"]["artifacts"]["base"]["path"])
    if hashlib.sha256(atom.read_bytes()).hexdigest() != manifest["ions"][f"O{stage}"]["artifacts"]["base"]["sha256"]:
        raise ValueError("Candidate atom changed after native validation")
    radiation = atom_radiation(atom)
    ordered = sorted(collisions, key=lambda tr: C.excitation_levels(levels, tr))
    atom_lines = atom.read_text().splitlines()
    start = next(i for i, line in enumerate(atom_lines) if line.startswith("NUPS "))
    if int(atom_lines[start].split()[1]) != len(ordered):
        raise ValueError("Candidate collision count differs from CHIANTI source")
    position = start + 1
    for tr in ordered:
        fields = atom_lines[position].split()
        count = int(fields[-1])
        lower, upper = C.excitation_levels(levels, tr)
        if (int(fields[0]), int(fields[1])) != (lower, upper):
            raise ValueError("Candidate collision identity differs from CHIANTI source")
        nodes = np.asarray([list(map(float, line.split()))
                            for line in atom_lines[position + 1:position + 1 + count]])
        if (not np.array_equal(nodes[:, 0], tr["xs"])
            or not np.array_equal(nodes[:, 1], tr["ys"])):
            raise ValueError("Candidate collision nodes differ from CHIANTI source")
        position += count + 1
    probe_temperature = np.geomspace(1000.0, 100000.0, 5)
    probe_grid = HERE / f"native_O{stage}_strength_grid.txt"
    probe_grid.write_text("1 5\n" + " ".join(format(t, ".17g") for t in probe_temperature) + "\n")
    probe_output = HERE / f"native_O{stage}_strengths.txt"
    subprocess.run([str(HERE / "native_source_check.x"), str(atom),
                    str(probe_grid), str(probe_output)], check=True, timeout=1200)
    strengths = np.loadtxt(probe_output)
    if strengths.shape != (5, len(ordered) + 1):
        raise ValueError("Native effective-strength output incomplete")
    reference_strength = np.column_stack([C.upsilon(tr, probe_temperature) for tr in ordered])
    if np.any(strengths[:, 1:][reference_strength == 0] != 0):
        raise ValueError("Native source zero strengths changed")
    positive = reference_strength > 0
    strength_error = float(np.max(abs(strengths[:, 1:][positive] / reference_strength[positive] - 1)))
    if strength_error > 1e-9:
        raise ValueError("Native effective strengths differ from CHIANTI source")
    raw = np.loadtxt(HERE / f"native_O{stage}_base.txt")
    nlev = len(levels)
    native_pop = raw[:, 3:].reshape(401, 4, nlev)
    native_cooling = raw[:, 2].reshape(401, 4)
    direct = full_excitation_balance(levels, radiation, collisions, TEMPERATURE, DENSITY)
    max_cooling = 0.0
    max_population = 0.0
    max_stationarity = 0.0
    max_energy = 0.0
    max_energy_gross = 0.0
    for j, density in enumerate(DENSITY):
        expected, direct_pop = direct[density]
        actual = native_cooling[:, j]
        maximum = float(np.max(np.abs(actual / expected - 1)))
        population = float(np.max(np.abs(native_pop[:, j, :] - direct_pop)))
        max_cooling = max(max_cooling, maximum)
        max_population = max(max_population, population)
        axis[0].loglog(TEMPERATURE, expected, label=f"n_e={density:.0e}")
        axis[1].semilogx(TEMPERATURE, 100 * (actual / expected - 1),
                            label=f"n_e={density:.0e}")
        saved[f"O{stage}_direct_cooling_ne{j}"] = expected
        saved[f"O{stage}_direct_pop_ne{j}"] = direct_pop
        saved[f"O{stage}_native_cooling_ne{j}"] = actual

        # Native populations are independently inserted into the unmodified
        # source collision strengths to check stationarity and thermal power.
        populations = native_pop[:, j, :]
        residual = np.zeros_like(populations)
        outgoing = np.zeros_like(populations)
        collision_power = np.zeros(401)
        gross_power = np.zeros(401)
        for tr in collisions:
            lower, upper = C.excitation_levels(levels, tr)
            gap = (levels[upper][1] - levels[lower][1]) * C.HC_ERG_CM
            if gap < 0:
                raise ValueError("Collision transition has negative physical gap")
            common = C.COLL_PREF * C.upsilon(tr, TEMPERATURE) / np.sqrt(TEMPERATURE)
            up = common / levels[lower][0] * np.exp(-gap / (C.K_B_ERG * TEMPERATURE))
            down = common / levels[upper][0]
            forward = density * populations[:, lower - 1] * up
            reverse = density * populations[:, upper - 1] * down
            net = forward - reverse
            residual[:, lower - 1] -= net
            residual[:, upper - 1] += net
            outgoing[:, lower - 1] += forward
            outgoing[:, upper - 1] += reverse
            collision_power += (populations[:, lower - 1] * up
                                - populations[:, upper - 1] * down) * gap
            gross_power += (populations[:, lower - 1] * up
                            + populations[:, upper - 1] * down) * gap
        for lower, upper, _, _, rate in radiation:
            flow = populations[:, upper - 1] * rate
            residual[:, lower - 1] += flow
            residual[:, upper - 1] -= flow
            outgoing[:, upper - 1] += flow
        stationarity = float(np.max(np.abs(residual) / np.maximum(np.max(outgoing, axis=1)[:, None], 1e-300)))
        energy = float(np.max(np.abs(collision_power / actual - 1)))
        energy_gross = float(np.max(np.abs(collision_power - actual)
                                    / np.maximum(gross_power, 1e-300)))
        max_stationarity = max(max_stationarity, stationarity)
        max_energy = max(max_energy, energy)
        max_energy_gross = max(max_energy_gross, energy_gross)
    axis[0].set_title(f"O {stage}: direct and native cooling")
    axis[1].set_title(f"O {stage}: native/direct")
    axis[0].set_ylabel("Cooling (erg cm$^3$ s$^{-1}$)")
    axis[1].set_ylabel("Difference (percent)")
    for panel in axis:
        panel.set_xlabel("Temperature (K)")
        panel.grid(True, alpha=.25)
        panel.legend(fontsize=7)
    report[f"O{stage}"] = {"source_collision_rows": len(ordered),
                            "source_collision_nodes": position - start - 1 - len(ordered),
                            "maximum_native_upsilon_relative_error": strength_error,
                            "maximum_native_cooling_relative_error": max_cooling,
                            "maximum_native_population_absolute_error": max_population,
                            "maximum_stationarity_residual_relative_to_outflow": max_stationarity,
                            "maximum_collision_to_radiation_relative_error": max_energy,
                            "maximum_energy_residual_relative_to_gross_collisions": max_energy_gross}
fig.tight_layout()
fig.savefig(HERE / "independent_balance_profiles.pdf")
np.savez_compressed(HERE / "independent_balance_profiles.npz", **saved)
passed = (all(x["maximum_native_cooling_relative_error"] < 1e-8
              and x["maximum_native_population_absolute_error"] < 1e-9
              and x["maximum_stationarity_residual_relative_to_outflow"] < 1e-8
              and x["maximum_energy_residual_relative_to_gross_collisions"] < 1e-8
              for x in report.values()))
result = {"passed": passed, "temperature_count": 401, "density_count": 4,
          "ions": report}
(HERE / "independent_balance_result.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
if not passed:
    raise ValueError("Independent balance checks failed")
