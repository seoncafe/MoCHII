"""Freeze a reviewed Si I candidate and its numerical validation for products."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = HERE / "candidate"
PLAN = HERE / "plan.json"
STATE = HERE / "controller"
VALIDATION_DIR = HERE.parent / "si_i_radiative_completion"
SOURCE_DIR = HERE.parent / "si_i_electric_dipole_source_audit"
FITTING = ROOT / "tools/fitting"
COOLING_BUILD = HERE.parent / "bt_floor_cooling_products"
REGISTRY_BUILD = HERE.parent / "oxygen_valence_convergence"
APPROVED_CANDIDATE_SHA256 = "2c47d6e54cf4d1cfbc6340c5f01837889008a8fa8a0ac13afe0169f7f5d824df"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite_small(value, limit):
    return (isinstance(value, (int, float)) and math.isfinite(value)
            and 0 <= value < limit)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    args = parser.parse_args()
    if any(path.exists() for path in (PLAN, STATE, BASE)):
        raise FileExistsError("Existing Si I plan or products must be preserved")
    candidate, validation_path = args.candidate.resolve(), args.validation.resolve()
    if candidate != (SOURCE_DIR / "candidate_v2/nlevel_si_1.txt").resolve():
        raise ValueError("Only the reviewed corrected Si I candidate is accepted")
    if validation_path != (VALIDATION_DIR / "validation_manifest.json").resolve():
        raise ValueError("Only the complete Si I numerical manifest is accepted")
    record = json.loads((SOURCE_DIR / "candidate_v2/candidate_result.json").read_text())
    candidate_sha = digest(candidate)
    if (record.get("passed") is not True or record.get("candidate_sha256") != candidate_sha
        or candidate_sha != APPROVED_CANDIDATE_SHA256
        or [record.get(name) for name in ("nlev", "nrad", "nups", "nodes")]
        != [35, 250, 595, 8925] or record.get("bt_tmin_K") != 1003.02):
        raise ValueError("Si I corrected candidate identity or counts changed")
    validation = json.loads(validation_path.read_text())
    if (validation.get("passed") is not True
        or validation.get("candidate_sha256") != candidate_sha
        or Path(validation.get("candidate_path", "")).resolve() != candidate
        or validation.get("temperature_points") != 401
        or validation.get("densities_cm3") != [1e-6, 100.0, 1e6, 1e10]
        or [validation.get(name) for name in ("nlev", "nrad", "nups", "nodes")]
        != [35, 250, 595, 8925]
        or validation.get("bt_tmin_K") != 1003.02
        or validation.get("ledger_numerical_checks_passed") is not True
        or validation.get("source_accuracy_verified") is not False):
        raise ValueError("Complete Si I source/native/ledger checks are required")
    for key in ("maximum_relative_cooling_errors", "maximum_absolute_population_errors"):
        values = validation.get(key)
        if not isinstance(values, list) or len(values) != 4 or any(
            not finite_small(value, 1e-10) for value in values
        ):
            raise ValueError(f"Invalid or incomplete Si I {key}")
    if not finite_small(validation.get("maximum_relative_strength_error"), 1e-10):
        raise ValueError("Si I source collision strengths are not validated")

    source_files = [HERE / name for name in
                    ("prepare_plan.py", "run_stage.py", "run_controller_sync.sh", "README.md")]
    source_files += [candidate, validation_path,
                     SOURCE_DIR / "candidate_v2/candidate_result.json",
                     Path(record["source_module_path"]),
                     COOLING_BUILD / "build_result.json",
                     COOLING_BUILD / "build_mpi/check_cooling_fit.x",
                     REGISTRY_BUILD / "optimized_build_manifest.json",
                     REGISTRY_BUILD / "optimized_build/registered_check.x"]
    source_files += [Path(path) for path in record["source_sha256"]]
    source_files += [SOURCE_DIR / name for name in (
        "audit_si_i_electric_dipole.py", "Si_14.16.cff-mchf-lin.dat.mp")]
    source_files += [Path.home() / ".codex/scripts/calculation_workflow.py"]
    cooling_build_record = json.loads((COOLING_BUILD / "build_result.json").read_text())
    source_files += [Path(path) for path in cooling_build_record["source_sha256"]]
    source_files += [COOLING_BUILD / "build_mpi" / name for name in (
        "define.o", "utility.o", "hydrogenic_recombination_mod.o", "recomb_mod.o",
        "gaunt_vh14_mod.o", "nlevel_mod.o", "nlevel_cooling_mod.o",
        "cooling_mod.o", "gaunt.o")]
    source_files += [REGISTRY_BUILD / "optimized_build" / name for name in (
        "define.o", "nlevel_mod.o", "nlevel_cooling_mod.o")]
    source_files += [FITTING / name for name in (
        "fit_cooling.py", "chianti_cooling.py", "atomic_excitation_energies.py",
        "ionization_potentials.py", "validate_native_cooling_coefficients.py",
        "validate_registered_excitation_atom.py", "validate_stage_cooling.py",
        "atomic_output.py", "write_excitation_atom.py")]
    source_files += [ROOT / "src" / name for name in (
        "define.f90", "nlevel_mod.f90", "nlevel_cooling_mod.f90", "cooling_mod.f90")]
    source_files += [ROOT / "tests/ion_stages" / name for name in (
        "check_cooling_fit.f90", "check_registered_excitation_atom.f90")]
    for path in source_files:
        if not path.is_file():
            raise FileNotFoundError(f"Missing source/input before freeze: {path}")
    steps = []
    for identifier, stage, result, limit in (
        ("SiI_atom", "atom", "atom_result.json", 180),
        ("SiI_fit", "fit", "fit_result.json", 3600),
        ("SiI_registry", "registry", "result.json", 1800),
    ):
        command = ["/usr/bin/prlimit", "--as=8589934592", f"--cpu={limit - 60}",
                   "--", sys.executable, str(HERE / "run_stage.py"), stage,
                   "--base", str(BASE)]
        if stage == "atom":
            command += ["--candidate", str(candidate)]
        steps.append({"id": identifier, "command": command,
                      "timeout_seconds": limit,
                      "check": {"path": str(BASE / result),
                                "equals": ({"passed": True, "ion": "SiI"}
                                           if stage == "registry" else
                                           {"passed": True, "stage": stage})}})
    plan = {"cwd": str(ROOT),
            "environment": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                            "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
                            "MPLCONFIGDIR": "/tmp/mochii_sii_fit_mpl"},
            "source_files": [str(path.resolve()) for path in dict.fromkeys(source_files)],
            "steps": steps}
    PLAN.write_text(json.dumps(plan, indent=2) + "\n")
    print(PLAN)


if __name__ == "__main__":
    main()
