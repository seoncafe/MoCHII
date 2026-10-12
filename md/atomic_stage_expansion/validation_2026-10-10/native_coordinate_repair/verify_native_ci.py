"""Compare the staged C I native atom with its five-level source balance."""
from pathlib import Path
import json
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as source
from validate_added_charge_transfer import compile_physics_check
from validate_excitation_density_profiles import full_excitation_balance

HERE = Path(__file__).resolve().parent


def main():
    (HERE / "native_check").mkdir(exist_ok=True)
    executable = compile_physics_check("tests/ion_stages/check_excitation_atom.f90",
                                       HERE / "native_check")
    native_output = HERE / "atom_native_ci.txt"
    subprocess.run([str(executable),
                    str(ROOT / "data/atomic/cooling_c_1.txt"),
                    str(HERE / "nlevel_c_1.txt"), str(native_output)], check=True)
    rows = np.loadtxt(native_output)
    assert rows.shape == (2403, 9), rows.shape
    levels, radiation, collisions = source.bound_excitation_atom("c", 1)
    keep = set(sorted(levels)[:5])
    levels = {i: levels[i] for i in sorted(keep)}
    radiation = [r for r in radiation if r[0] in keep and r[1] in keep]
    collisions = [r for r in collisions if r["ll"] in keep and r["ul"] in keep]
    temperature = rows[::3, 0]
    density = [1e-6, 100.0, 1e6]
    reference = full_excitation_balance(levels, radiation, collisions,
                                        temperature, density)
    result = {"native_rows": len(rows), "temperature_count": len(temperature),
              "densities_cm3": density, "comparisons": []}
    for k, ne in enumerate(density):
        source_cooling, source_population = reference[ne]
        native_cooling = rows[k::3, 3]
        native_population = rows[k::3, 4:]
        difference = np.abs(native_cooling/source_cooling-1)
        population_difference = np.abs(native_population-source_population)
        result["comparisons"].append({
            "density_cm3": ne,
            "maximum_relative_cooling_difference": float(difference.max()),
            "worst_temperature_K": float(temperature[difference.argmax()]),
            "maximum_absolute_population_difference": float(population_difference.max())})
    (HERE / "native_ci_result.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
