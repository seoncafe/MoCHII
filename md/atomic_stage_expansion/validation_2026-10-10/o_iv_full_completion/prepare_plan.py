"""Freeze sources and finite successors for conditional full O IV diagnostics."""

import hashlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = HERE / "candidate"
FITTING = ROOT / "tools/fitting"
OXYGEN_BUILD = HERE.parent / "oxygen_valence_convergence"
COOLING_BUILD = HERE.parent / "bt_floor_cooling_products"
sys.path.insert(0, str(FITTING))
import chianti_cooling as C  # noqa: E402
from write_excitation_atom import _source_files, _selection_code_files  # noqa: E402


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    plan_path = HERE / "plan.json"
    if any(path.exists() for path in (plan_path, HERE / "controller", BASE)):
        raise FileExistsError("Preserve existing O IV plan, controller, and products")
    levels, radiation, collisions = C.bound_excitation_atom("o", 4)
    if (len(levels), len(radiation), len(collisions),
        sum(len(row["xs"]) for row in collisions)) != (147, 2319, 10731, 160965):
        raise ValueError("Current full O IV bound source size changed")
    oxygen_manifest = OXYGEN_BUILD / "optimized_build_manifest.json"
    cooling_manifest = COOLING_BUILD / "build_result.json"
    oxygen_build = json.loads(oxygen_manifest.read_text())
    cooling_build = json.loads(cooling_manifest.read_text())
    oxygen_files = [
        "src/define.f90", "src/nlevel_mod.f90", "src/nlevel_cooling_mod.f90",
        "tests/ion_stages/check_source_excitation_atom.f90",
        "tests/ion_stages/check_registered_excitation_atom.f90",
        *[f"md/atomic_stage_expansion/validation_2026-10-10/oxygen_valence_convergence/"
          f"optimized_build/{name}" for name in (
              "source_check.x", "registered_check.x", "define.o",
              "nlevel_mod.o", "nlevel_cooling_mod.o")],
    ]
    for relative in oxygen_files:
        if oxygen_build["sha256"].get(relative) != digest(ROOT / relative):
            raise ValueError(f"Optimized O IV executable or linked source changed: {relative}")
    cooling_executable = COOLING_BUILD / "build_mpi/check_cooling_fit.x"
    if (cooling_build.get("passed") is not True
        or cooling_build["executables"].get("check_cooling_fit.x") != digest(cooling_executable)
        or any(cooling_build["source_sha256"].get(str((ROOT / name).resolve()))
               != digest(ROOT / name)
               for name in ("src/define.f90", "src/cooling_mod.f90",
                            "tests/ion_stages/check_cooling_fit.f90"))):
        raise ValueError("Actual GNU coefficient reader changed")
    source_files = [HERE / name for name in (
        "prepare_plan.py", "run_stage.py", "check_energy_balance.py",
        "launch_controller.sh", "README.md")]
    source_files += [*_source_files("o", 4), *_selection_code_files("o", 4)]
    source_files += [FITTING / name for name in (
        "fit_cooling.py", "chianti_cooling.py", "write_excitation_atom.py",
        "valence_excitation_energies.py", "electric_dipole_probabilities.py",
        "atomic_excitation_energies.py", "ionization_potentials.py",
        "validate_source_excitation_atom.py", "validate_excitation_density_profiles.py",
        "validate_native_cooling_coefficients.py",
        "validate_registered_excitation_atom.py", "validate_stage_cooling.py",
        "atomic_output.py")]
    source_files += [ROOT / relative for relative in oxygen_files]
    source_files += [cooling_manifest, oxygen_manifest, cooling_executable]
    source_files += [COOLING_BUILD / "build_mpi" / name for name in (
        "define.o", "utility.o", "hydrogenic_recombination_mod.o",
        "recomb_mod.o", "gaunt_vh14_mod.o", "nlevel_mod.o",
        "nlevel_cooling_mod.o", "cooling_mod.o", "gaunt.o")]
    source_files += [Path(name) for name in cooling_build["source_sha256"]]
    source_files += [Path.home() / ".codex/scripts/calculation_workflow.py"]
    unique = list(dict.fromkeys(path.resolve() for path in source_files))
    if any(not path.is_file() for path in unique):
        missing = [str(path) for path in unique if not path.is_file()]
        raise FileNotFoundError(f"Missing frozen source inputs: {missing}")
    stages = (
        ("OIV_preflight", [sys.executable, str(HERE / "run_stage.py"), "preflight"],
         "preflight_result.json", 180),
        ("OIV_native_401x4", [sys.executable, str(HERE / "run_stage.py"), "native"],
         "native_result.json", 2400),
        ("OIV_energy_ledger", [sys.executable, str(HERE / "check_energy_balance.py")],
         "ledger_result.json", 900),
        ("OIV_fit_and_GNU_reader", [sys.executable, str(HERE / "run_stage.py"), "fit"],
         "fit_result.json", 3600),
        ("OIV_Intel_registry", [sys.executable, str(HERE / "run_stage.py"), "registry"],
         "result.json", 1800),
    )
    steps = []
    for identifier, command, result, limit in stages:
        target = BASE / result
        check = {"path": str(target), "equals": {"passed": True}}
        if result == "result.json":
            check["equals"]["ion"] = "OIV"
        steps.append({"id": identifier,
                      "command": ["/usr/bin/prlimit", "--as=8589934592",
                                  f"--cpu={limit - 60}", "--", *command],
                      "timeout_seconds": limit, "check": check})
    plan = {"cwd": str(ROOT),
            "environment": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                            "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
                            "I_MPI_FABRICS": "shm",
                            "MPLCONFIGDIR": "/tmp/mochii_oiv_full_fit_mpl"},
            "source_files": [str(path) for path in unique], "steps": steps}
    plan_path.write_text(json.dumps(plan, indent=2) + "\n")
    print(plan_path)


if __name__ == "__main__":
    main()
