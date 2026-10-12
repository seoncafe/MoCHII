"""Prepare finite Ne I source and private cooling-product validation steps."""

import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OLD = HERE.parent / "source_validity_ne_i_mg_i/ne_i_candidate/nlevel_ne_1.txt"
OXYGEN = HERE.parent / "oxygen_valence_convergence/optimized_build"
PRODUCTS = HERE.parent / "ne_i_cooling_products"
COOLING = HERE.parent / "bt_floor_cooling_products/build_mpi"
PYTHON = Path(sys.executable)


def stage(identifier, action, result, seconds):
    return {
        "id": identifier,
        "command": ["/usr/bin/prlimit", "--as=8589934592",
                    f"--cpu={seconds - 30}", "--", str(PYTHON),
                    str(HERE / "run_profiles.py"), "--step", action],
        "timeout_seconds": seconds,
        "check": {"path": str(result), "equals": {"passed": True}},
    }


def main():
    target = HERE / "plan.json"
    if target.exists():
        raise FileExistsError(f"Preserve existing Ne I plan: {target}")
    raw = Path.home() / "RT_Codes/CHIANTI/dbase/ne/ne_1/ne_1"
    stout = Path.home() / "CLOUDY/c25.00/data/stout/ne/ne_1/ne_1"
    source_files = [
        HERE / "prepare_candidate.py", HERE / "candidate_result.json",
        HERE / "atom/nlevel_ne_1.txt", HERE / "run_profiles.py",
        HERE / "verify_independent_balance.py", HERE / "validate_complete_profiles.py",
        HERE / "launch_controller.sh", OLD,
        OXYGEN / "source_check.x", OXYGEN / "registered_check.x",
        OXYGEN / "nlevel_mod.o", OXYGEN / "nlevel_cooling_mod.o",
        OXYGEN / "define.o",
        HERE.parent / "oxygen_valence_convergence/optimized_build_manifest.json",
        COOLING / "check_cooling_fit.x",
        HERE.parent / "bt_floor_cooling_products/build_result.json",
        ROOT / "src/nlevel_mod.f90", ROOT / "src/nlevel_cooling_mod.f90",
        ROOT / "src/define.f90",
        ROOT / "tests/ion_stages/check_source_excitation_atom.f90",
        ROOT / "tests/ion_stages/check_nei_spontaneous_probabilities.py",
        ROOT / "tools/fitting/chianti_cooling.py",
        ROOT / "tools/fitting/write_excitation_atom.py",
        ROOT / "tools/fitting/neon_spontaneous_probabilities.py",
        ROOT / "tools/fitting/atomic_excitation_energies.py",
        ROOT / "tools/fitting/ionization_potentials.py",
        ROOT / "tools/fitting/validate_excitation_density_profiles.py",
        ROOT / "tools/fitting/fit_cooling.py",
        ROOT / "tools/fitting/validate_native_cooling_coefficients.py",
        ROOT / "tools/fitting/validate_registered_excitation_atom.py",
        PRODUCTS / "run_controller_sync.sh", PRODUCTS / "prepare_plan.py",
        PRODUCTS / "run_stage.py", Path.home() / ".codex/scripts/calculation_workflow.py",
        Path(str(raw) + ".elvlc"), Path(str(raw) + ".wgfa"),
        Path(str(raw) + ".scups"),
        Path(str(stout) + ".nrg"), Path(str(stout) + ".tp"),
        Path.home() / "CLOUDY/c25.00/data/ionization_potentials.dat",
    ]
    for path in source_files:
        if not path.is_file():
            raise FileNotFoundError(f"Missing Ne I input before freeze: {path}")
    steps = [
        stage("pilot", "pilot", HERE / "pilot_step.json", 180),
        stage("native", "native", HERE / "native_step.json", 1200),
        stage("direct", "direct", HERE / "direct_step.json", 1200),
        stage("compare", "compare", HERE / "comparison_step.json", 300),
        {"id": "balance",
         "command": ["/usr/bin/prlimit", "--as=8589934592", "--cpu=870",
                     "--", str(PYTHON), str(HERE / "verify_independent_balance.py")],
         "timeout_seconds": 900,
         "check": {"path": str(HERE / "independent_balance_result.json"),
                   "equals": {"passed": True}}},
        {"id": "manifest",
         "command": ["/usr/bin/prlimit", "--as=8589934592", "--cpu=270",
                     "--", str(PYTHON), str(HERE / "validate_complete_profiles.py")],
         "timeout_seconds": 300,
         "check": {"path": str(HERE / "validation_manifest.json"),
                   "equals": {"passed": True,
                              "stout_added_radiative_pairs_checked": True,
                              "ledger_numerical_checks_passed": True}}},
        {"id": "private_products",
         "command": ["/usr/bin/env", "NEI_PRODUCTS_APPROVED=yes",
                     f"NEI_CORRECTED_CANDIDATE={HERE / 'atom/nlevel_ne_1.txt'}",
                     f"NEI_VALIDATION_MANIFEST={HERE / 'validation_manifest.json'}",
                     "bash", str(PRODUCTS / "run_controller_sync.sh")],
         "timeout_seconds": 7200,
         "check": {"path": str(PRODUCTS / "candidate/result.json"),
                   "equals": {"passed": True, "ion": "NeI"}}},
    ]
    target.write_text(json.dumps({
        "cwd": str(ROOT),
        "environment": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                        "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
                        "MPLCONFIGDIR": "/tmp/mochii_nei_complete_mpl"},
        "source_files": [str(path.resolve()) for path in source_files],
        "steps": steps,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
