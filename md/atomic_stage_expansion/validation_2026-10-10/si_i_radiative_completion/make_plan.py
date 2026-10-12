"""Prepare a finite, source-hashed Si I numerical validation controller."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BUILD = HERE.parent / "oxygen_valence_convergence/optimized_build"
GAUGE = HERE.parent / "si_i_radiative_gauge_sensitivity"
PRODUCTS = HERE.parent / "si_i_cooling_products"
COOLING = HERE.parent / "bt_floor_cooling_products"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def step(identifier: str, command: list[str], result: str,
         timeout: int, expected: dict | None = None) -> dict:
    return {
        "id": identifier,
        "command": ["/usr/bin/prlimit", "--as=8589934592",
                    f"--cpu={timeout - 30}", "--", *command],
        "timeout_seconds": timeout,
        "check": {"path": str(HERE / result),
                  "equals": expected or {"passed": True}},
    }


def main() -> None:
    target = HERE / "plan.json"
    if target.exists():
        raise FileExistsError(f"Preserve existing Si I controller plan: {target}")
    candidate_path = (HERE.parent / "si_i_electric_dipole_source_audit/candidate_v2/"
                      "candidate_result.json")
    candidate = json.loads(candidate_path.read_text())
    if candidate.get("passed") is not True:
        raise ValueError("Si I candidate source checks have not passed")
    atom = Path(candidate["candidate_path"])
    module = Path(candidate["source_module_path"])
    if digest(atom) != candidate["candidate_sha256"] or digest(module) != candidate[
            "source_sha256"][str(module)]:
        raise ValueError("Si I candidate or source module changed")
    if (candidate["nlev"], candidate["nrad"], candidate["nups"],
            candidate["nodes"]) != (35, 250, 595, 8925):
        raise ValueError("Incomplete Si I candidate atom")
    source_files = [
        candidate_path, atom, module,
        HERE / "run_profiles.py", HERE / "check_energy_balance.py",
        HERE / "validate_profiles.py", HERE / "make_plan.py",
        HERE / "launch_controller.sh",
        HERE.parent / "si_i_electric_dipole_source_audit/result.json",
        HERE.parent / "si_i_electric_dipole_source_audit/audit_si_i_electric_dipole.py",
        HERE.parent / "si_i_electric_dipole_source_audit/Si_14.16.cff-mchf.lsj.db",
        HERE.parent / "si_i_electric_dipole_source_audit/Si_14.16.cff-mchf-lin.dat.mp",
        BUILD / "source_check.x", BUILD / "nlevel_mod.o", BUILD / "define.o",
        HERE.parent / "oxygen_valence_convergence/optimized_build_manifest.json",
        ROOT / "src/nlevel_mod.f90", ROOT / "src/define.f90",
        ROOT / "tests/ion_stages/check_source_excitation_atom.f90",
        ROOT / "tools/fitting/chianti_cooling.py",
        ROOT / "tools/fitting/validate_excitation_density_profiles.py",
        ROOT / "tools/fitting/write_excitation_atom.py",
        ROOT / "tools/fitting/atomic_output.py",
        ROOT / "tools/fitting/atomic_excitation_energies.py",
        ROOT / "tools/fitting/niv_level_energies.py",
        ROOT / "tools/fitting/ionization_potentials.py",
        ROOT / "tools/fitting/fit_cooling.py",
        ROOT / "tools/fitting/validate_native_cooling_coefficients.py",
        ROOT / "tools/fitting/validate_registered_excitation_atom.py",
        ROOT / "tools/fitting/validate_stage_cooling.py",
        GAUGE / "run_gauge_sensitivity.py", GAUGE / "README.md",
        PRODUCTS / "prepare_plan.py", PRODUCTS / "run_stage.py",
        PRODUCTS / "run_controller_sync.sh", PRODUCTS / "README.md",
        COOLING / "build_result.json", COOLING / "build_mpi/check_cooling_fit.x",
        BUILD / "registered_check.x", BUILD / "nlevel_cooling_mod.o",
        ROOT / "src/nlevel_cooling_mod.f90", ROOT / "src/cooling_mod.f90",
        Path.home() / ".codex/scripts/calculation_workflow.py",
        *(Path(name) for name in candidate["source_sha256"]),
    ]
    cooling_record = json.loads((COOLING / "build_result.json").read_text())
    if cooling_record.get("passed") is not True:
        raise ValueError("Isolated Si I cooling reader build is not validated")
    source_files.extend(Path(name) for name in cooling_record["source_sha256"])
    source_files.extend(COOLING / "build_mpi" / name for name in (
        "define.o", "utility.o", "hydrogenic_recombination_mod.o",
        "recomb_mod.o", "gaunt_vh14_mod.o", "nlevel_mod.o",
        "nlevel_cooling_mod.o", "cooling_mod.o", "gaunt.o"))
    source_files.extend(BUILD / name for name in ("define.o", "nlevel_mod.o",
                                                "nlevel_cooling_mod.o"))
    unique = []
    for path in source_files:
        resolved = path.resolve()
        if resolved not in unique:
            unique.append(resolved)
    for path in unique:
        if not path.is_file():
            raise FileNotFoundError(f"Missing Si I source before freeze: {path}")
    for name, expected in candidate["source_sha256"].items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Si I raw source changed before freeze: {name}")
    runner = [sys.executable, str(HERE / "run_profiles.py")]
    steps = [
        step("pilot", [*runner, "--step", "pilot"], "pilot_step.json", 180),
        step("native", [*runner, "--step", "native"], "native_step.json", 1200),
        step("direct", [*runner, "--step", "direct"], "direct_step.json", 1200),
        step("compare", [*runner, "--step", "compare"], "comparison_step.json", 300),
        step("energy_balance", [sys.executable, str(HERE / "check_energy_balance.py")],
             "energy_balance_step.json", 600),
        step("manifest", [sys.executable, str(HERE / "validate_profiles.py")],
             "validation_manifest.json", 180),
        step("gauge_sensitivity", [sys.executable,
                                   str(GAUGE / "run_gauge_sensitivity.py")],
             str(GAUGE / "result.json"), 1200),
        step("private_products", ["/usr/bin/env", "SII_PRODUCTS_APPROVED=yes",
                                  f"SII_CORRECTED_CANDIDATE={atom}",
                                  f"SII_VALIDATION_MANIFEST={HERE / 'validation_manifest.json'}",
                                  "bash", str(PRODUCTS / "run_controller_sync.sh")],
             str(PRODUCTS / "candidate/result.json"), 7200,
             {"passed": True, "ion": "SiI"}),
    ]
    target.write_text(json.dumps({
        "cwd": str(ROOT),
        "environment": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                        "MKL_NUM_THREADS": "1", "NUMEXPR_NUM_THREADS": "1",
                        "MPLCONFIGDIR": "/tmp/mochii_sii_radiative_mpl"},
        "source_files": [str(path) for path in unique],
        "steps": steps,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
