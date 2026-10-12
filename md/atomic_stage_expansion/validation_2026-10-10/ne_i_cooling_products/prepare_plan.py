"""Freeze the approved Ne I product inputs immediately before launch."""

import json
import argparse
import hashlib
import math
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = HERE / "candidate"
PLAN = HERE / "plan.json"
STATE = HERE / "controller"
DRIVER_DIRECTORY = HERE.parent / "oxygen_valence_convergence/optimized_build"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    args = parser.parse_args()
    if PLAN.exists() or STATE.exists() or BASE.exists():
        raise FileExistsError("Existing Ne I plan or product; preserve it and inspect state")
    candidate = args.candidate.resolve()
    validation_path = args.validation.resolve()
    if not candidate.is_file() or not validation_path.is_file():
        raise FileNotFoundError("Corrected candidate and post-correction validation are required")
    from run_stage import STALE_CANDIDATE_SHA256, identified_radiative_rates
    candidate_sha = hashlib.sha256(candidate.read_bytes()).hexdigest()
    if candidate_sha == STALE_CANDIDATE_SHA256:
        raise ValueError("Former 31/142 Ne I atom is physically stale")
    identified_radiative_rates(candidate)
    validation = json.loads(validation_path.read_text())
    if (validation.get("passed") is not True
        or validation.get("candidate_sha256") != candidate_sha
        or validation.get("temperature_points") != 401
        or validation.get("densities_cm3") != [1e-6, 100.0, 1e6, 1e10]
        or validation.get("stout_added_radiative_pairs_checked") is not True
        or validation.get("ledger_numerical_checks_passed") is not True):
        raise ValueError("Post-correction source identity and native numerical checks are incomplete")
    errors = validation.get("maximum_relative_cooling_errors")
    if (not isinstance(errors, list) or len(errors) != 4
        or any(not isinstance(error, (float, int)) or not math.isfinite(error)
               or error < 0 or error >= 1e-12
               for error in errors)):
        raise ValueError("Corrected candidate has no accepted 401-temperature native profile")
    strength_error = validation.get("maximum_relative_strength_error")
    if (not isinstance(strength_error, (float, int))
        or not math.isfinite(strength_error) or strength_error < 0
        or strength_error >= 1e-12):
        raise ValueError("Corrected candidate has no accepted collision-strength profile")
    fitting = ROOT / "tools/fitting"
    source_files = [
        HERE / "prepare_plan.py", HERE / "run_stage.py", HERE / "launch_controller.sh",
        HERE / "run_controller_sync.sh",
        candidate, validation_path,
        DRIVER_DIRECTORY / "registered_check.x",
        HERE.parent / "bt_floor_cooling_products/build_result.json",
        HERE.parent / "bt_floor_cooling_products/build_mpi/check_cooling_fit.x",
    ]
    source_files += [HERE.parent / "bt_floor_cooling_products/build_mpi" / name
                     for name in ("define.o", "utility.o", "hydrogenic_recombination_mod.o",
                                  "recomb_mod.o", "gaunt_vh14_mod.o", "nlevel_mod.o",
                                  "nlevel_cooling_mod.o", "cooling_mod.o", "gaunt.o")]
    source_files += [HERE.parent / "oxygen_valence_convergence/optimized_build_manifest.json"]
    source_files += [DRIVER_DIRECTORY / name for name in (
        "define.o", "nlevel_mod.o", "nlevel_cooling_mod.o")]
    source_files += [ROOT / "tests/ion_stages" / name for name in (
        "check_cooling_fit.f90", "check_registered_excitation_atom.f90")]
    source_files += [fitting / name for name in (
        "fit_cooling.py", "chianti_cooling.py", "atomic_excitation_energies.py",
        "neon_spontaneous_probabilities.py", "write_excitation_atom.py",
        "ionization_potentials.py", "validate_stage_cooling.py",
        "validate_native_cooling_coefficients.py",
        "validate_registered_excitation_atom.py", "atomic_output.py")]
    source_files += [ROOT / "src" / name for name in (
        "define.f90", "nlevel_mod.f90", "nlevel_cooling_mod.f90", "cooling_mod.f90")]
    stem = Path.home() / "RT_Codes/CHIANTI/dbase/ne/ne_1/ne_1"
    source_files += [Path(str(stem) + "." + suffix)
                     for suffix in ("elvlc", "wgfa", "scups")]
    stout_stem = Path.home() / "CLOUDY/c25.00/data/stout/ne/ne_1/ne_1"
    source_files += [Path(str(stout_stem) + "." + suffix)
                     for suffix in ("nrg", "tp")]
    sys.path.insert(0, str(fitting))
    from ionization_potentials import SOURCE as ionization_source
    source_files.append(Path(ionization_source))
    for path in source_files:
        if not path.is_file():
            raise FileNotFoundError(f"Missing input before freeze: {path}")
    stages = (("NeI_atom", "atom", "atom_result.json", 600),
              ("NeI_fit", "fit", "fit_result.json", 3600),
              ("NeI_registry", "registry", "result.json", 1800))
    steps = []
    for identifier, stage, result, limit in stages:
        check = ({"path": str(BASE / result),
                  "equals": {"passed": True, "stage": stage}}
                 if stage != "registry" else
                 {"path": str(BASE / result),
                  "equals": {"passed": True, "ion": "NeI"}})
        command = ["/usr/bin/prlimit", "--as=8589934592",
                   f"--cpu={limit - 60}", "--", sys.executable,
                   str(HERE / "run_stage.py"), stage, "--base", str(BASE)]
        if stage == "atom":
            command += ["--candidate", str(candidate)]
        steps.append({
            "id": identifier,
            "command": command,
            "timeout_seconds": limit,
            "check": check,
        })
    plan = {"cwd": str(ROOT),
            "environment": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                            "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
                            "MPLCONFIGDIR": "/tmp/mochii_nei_fit_mpl"},
            "source_files": [str(path.resolve()) for path in source_files],
            "steps": steps}
    PLAN.write_text(json.dumps(plan, indent=2) + "\n")
    print(PLAN)


if __name__ == "__main__":
    main()
