"""Generate and check private Si I cooling products from the corrected source."""

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
FITTING = ROOT / "tools/fitting"
sys.path.insert(0, str(FITTING))
import chianti_cooling as C  # noqa: E402
import fit_cooling as F  # noqa: E402
from atomic_output import atomic_write_text  # noqa: E402
from validate_stage_cooling import coefficients  # noqa: E402


SOURCE_DIR = HERE.parent / "si_i_electric_dipole_source_audit"
MODULE = SOURCE_DIR / "si_i_spontaneous_probabilities.py"
sys.path.insert(0, str(SOURCE_DIR))
COOLING_BUILD = HERE.parent / "bt_floor_cooling_products"
COOLING_DRIVER = COOLING_BUILD / "build_mpi/check_cooling_fit.x"
REGISTRY_BUILD = HERE.parent / "oxygen_valence_convergence"
REGISTRY_DRIVER = REGISTRY_BUILD / "optimized_build/registered_check.x"
COUNTS = (35, 250, 595, 8925)
FLOOR_K = 1003.02
FIT_TARGET = 0.003
REGISTRY_TARGET = 0.003


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def numerical_records(path):
    return b"\n".join(line for line in Path(path).read_bytes().splitlines()
                      if line.strip() and not line.lstrip().startswith(b"#")) + b"\n"


def result_file(path, data):
    if path.exists():
        raise FileExistsError(f"Preserve existing product: {path}")
    path.write_text(json.dumps(data, indent=2) + "\n")


def checked_source_atom():
    record = json.loads((SOURCE_DIR / "candidate_v2/candidate_result.json").read_text())
    if (record["candidate_sha256"] != digest(record["candidate_path"])
        or record["source_sha256"].get(str(MODULE)) != digest(MODULE)
        or tuple(record[name] for name in ("nlev", "nrad", "nups", "nodes")) != COUNTS
        or record["bt_tmin_K"] != FLOOR_K):
        raise ValueError("Si I corrected source identity changed")
    for name, expected in record["source_sha256"].items():
        if digest(name) != expected:
            raise ValueError(f"Si I source dependency changed: {name}")
    spec = importlib.util.spec_from_file_location("private_si_i_source", MODULE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    atom = module.source_atom()
    levels, radiation, collisions = atom
    if (len(levels), len(radiation), len(collisions), sum(len(row["xs"]) for row in collisions)) != COUNTS:
        raise ValueError("Private Si I source counts changed")
    if any(row.get("bt_tmin_K") != FLOOR_K for row in collisions):
        raise ValueError("Private Si I source temperature floor changed")
    return atom, record


def candidate_fit_file(name, elem, stage, amplitudes, temperatures,
                       error, relevant_error, output_dir):
    if (name, elem, stage) != ("SiI", "si", 1):
        raise ValueError("Only corrected Si I may use this source override")
    path = Path(output_dir) / "cooling_si_1.txt"
    lines = [
        "# Si I full 35-state low-density electron/spontaneous cooling fit",
        "# source: private corrected Si I candidate_v2 and si_i_spontaneous_probabilities.py",
        "# radiation: 229 existing CHIANTI pairs and 21 missing NIST-derived pairs; 250 total",
        "# collisions: 595 original CHIANTI pairs and 8925 original scaled source nodes",
        "# adopted numerical BT floor: 1003.02 K constant Upsilon below; Maxwellian rates use actual gas T",
        "# numerical fit domain: 1000--100000 K; source accuracy and 2005 final-paper methods unverified",
        "# optically thin stationary low-density limit; higher-density suppression is in the n-level atom",
        "# form: Lambda(T)=T^-1/2 sum_i A_i exp(-T_i/T) [erg cm^3 s^-1] per (n_e*n_ion)",
        "# signed coefficients represent one total positive cooling function",
        f"# maximum fitted relative difference on 201 fit temperatures: {error:.12g}",
        f"# maximum relevant relative difference: {relevant_error:.12g}",
        f"# candidate SHA-256: {digest(SOURCE_DIR / 'candidate_v2/nlevel_si_1.txt')}",
        f"# private source module SHA-256: {digest(MODULE)}",
        str(len(amplitudes)),
    ]
    lines += [f"  {amplitude:.16e}  {temperature:.16e}"
              for amplitude, temperature in zip(amplitudes, temperatures)]
    atomic_write_text(path, "\n".join(lines) + "\n")
    return str(path)


def atom_stage(base, candidate):
    _, record = checked_source_atom()
    if candidate.resolve() != Path(record["candidate_path"]).resolve():
        raise ValueError("Candidate path differs from corrected source record")
    output = base / "atom"
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    atom = output / "nlevel_si_1.txt"
    shutil.copyfile(candidate, atom)
    if digest(atom) != record["candidate_sha256"] or numerical_records(atom) != numerical_records(candidate):
        raise ValueError("Copied Si I atom differs from validated candidate")
    result_file(base / "atom_result.json", dict(
        passed=True, stage="atom", counts=COUNTS, bt_tmin_K=FLOOR_K,
        candidate_sha256=digest(candidate), atom_sha256=digest(atom),
        numerical_records_sha256=hashlib.sha256(numerical_records(atom)).hexdigest(),
        selected_ids=record["selected_ids"],
        density_ratio_cooling=record["density_ratio_cooling"]))


def fit_stage(base):
    atom_result = json.loads((base / "atom_result.json").read_text())
    atom, _ = checked_source_atom()
    if (atom_result.get("passed") is not True or tuple(atom_result["counts"]) != COUNTS
        or atom_result["atom_sha256"] != digest(base / "atom/nlevel_si_1.txt")):
        raise ValueError("Validated Si I source atom is required")
    output = base / "fit"
    if output.exists():
        raise FileExistsError(output)
    output.mkdir()
    original_bound = C.bound_excitation_atom
    original_fit_bound = F.bound_excitation_atom
    original_fit_cooling = F.cooling_low_density_limit
    original_fit_log = F.log_cooling_ground_state_limit
    original_writer = F.write_ion_file
    original_target = F.ERR_TARGET

    def selected_atom(elem, stage):
        if (elem, stage) != ("si", 1):
            raise ValueError("Private Si I fit cannot access another ion")
        return atom

    try:
        C.bound_excitation_atom = selected_atom
        F.bound_excitation_atom = selected_atom
        F.cooling_low_density_limit = C.cooling_low_density_limit
        F.log_cooling_ground_state_limit = C.log_cooling_ground_state_limit
        F.write_ion_file = candidate_fit_file
        F.ERR_TARGET = FIT_TARGET
        fitted_records = F.generate({"SiI": ("si", 1)}, str(output), plot_path=None)
        temperatures = np.geomspace(1000.0, 100000.0, 801)
        source = C.cooling_low_density_limit("si", 1, temperatures)
        if np.any(source == 0):
            log_source = C.log_cooling_ground_state_limit("si", 1, temperatures)
        else:
            log_source = np.log(source)
    finally:
        C.bound_excitation_atom = original_bound
        F.bound_excitation_atom = original_fit_bound
        F.cooling_low_density_limit = original_fit_cooling
        F.log_cooling_ground_state_limit = original_fit_log
        F.write_ion_file = original_writer
        F.ERR_TARGET = original_target
    fit_path = output / "cooling_si_1.txt"
    amplitude, excitation = coefficients(fit_path)
    fitted = F.evaluate(amplitude, excitation, temperatures)
    with np.errstate(divide="ignore", invalid="ignore"):
        difference = np.expm1(np.log(fitted) - log_source)
    if (not np.all(np.isfinite(log_source)) or not np.all(np.isfinite(fitted))
        or np.any(fitted <= 0) or not np.all(np.isfinite(difference))):
        raise ValueError("Nonfinite or nonpositive Si I cooling fit")
    error = float(np.max(np.abs(difference)))
    if error >= FIT_TARGET:
        raise ValueError(f"Si I 801-temperature source fit exceeds 0.3%: {error}")
    profile = output / "fit_profiles.npz"
    np.savez_compressed(profile, temperature_K=temperatures,
                        source_cooling=source, fitted_cooling=fitted,
                        relative_difference=difference)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(6.5, 6.5),
                             gridspec_kw={"height_ratios": [3, 1]})
    axes[0].loglog(temperatures, source, label="Corrected full 35-state source")
    axes[0].loglog(temperatures, fitted, "--", label="Serialized cooling fit")
    axes[0].set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
    axes[0].legend()
    axes[1].semilogx(temperatures, 100.0 * difference)
    axes[1].set_xlabel(r"$T$ [K]")
    axes[1].set_ylabel(r"Fit difference [\%]")
    fig.tight_layout()
    plot = output / "SiI_fit.pdf"
    fig.savefig(plot)
    plt.close(fig)
    build = json.loads((COOLING_BUILD / "build_result.json").read_text())
    if (build.get("passed") is not True
        or build["executables"].get("check_cooling_fit.x") != digest(COOLING_DRIVER)
        or any(build["source_sha256"].get(str((ROOT / "src" / name).resolve()))
               != digest(ROOT / "src" / name)
               for name in ("define.f90", "cooling_mod.f90"))
        or build["source_sha256"].get(
            str((ROOT / "tests/ion_stages/check_cooling_fit.f90").resolve()))
            != digest(ROOT / "tests/ion_stages/check_cooling_fit.f90")):
        raise ValueError("Compiled cooling coefficient reader changed")
    native_dir = base / "native_reader"
    command = [sys.executable, str(FITTING / "validate_native_cooling_coefficients.py"),
               "--coefficients", str(fit_path), "--driver", str(COOLING_DRIVER),
               "--work-dir", str(native_dir)]
    with (output / "native_reader.log").open("w") as stream:
        subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                       check=True, timeout=180)
    native = json.loads((native_dir / "result.json").read_text())
    if (native.get("passed") is not True or native.get("points") != 801
        or not math.isfinite(native["maximum_resolved_relative_error"])
        or native["maximum_resolved_relative_error"] >= 1e-9):
        raise ValueError("Compiled coefficient reader failed")
    result_file(base / "fit_result.json", dict(
        passed=True, stage="fit", temperature_points=801,
        temperature_K=[1000.0, 100000.0], maximum_relative_fit_error=error,
        coefficients_sha256=digest(fit_path), fit_plot_sha256=digest(plot),
        fit_profiles_sha256=digest(profile),
        source_atom_sha256=atom_result["atom_sha256"], native_reader=native,
        compiled_reader_sha256=digest(COOLING_DRIVER),
        compiled_reader_build_manifest_sha256=digest(COOLING_BUILD / "build_result.json"),
        fit_term_count=len(amplitude), source_fit_temperatures=len(F.T)))


def registry_stage(base):
    atom = json.loads((base / "atom_result.json").read_text())
    fit = json.loads((base / "fit_result.json").read_text())
    if (atom.get("passed") is not True or fit.get("passed") is not True
        or atom["atom_sha256"] != digest(base / "atom/nlevel_si_1.txt")
        or fit["coefficients_sha256"] != digest(base / "fit/cooling_si_1.txt")):
        raise ValueError("Si I atom and coefficients must pass first")
    build_manifest = REGISTRY_BUILD / "optimized_build_manifest.json"
    build = json.loads(build_manifest.read_text())
    linked = [REGISTRY_DRIVER] + [REGISTRY_DRIVER.parent / name for name in (
        "define.o", "nlevel_mod.o", "nlevel_cooling_mod.o")]
    for path in linked:
        relative = str(path.relative_to(ROOT))
        if build["sha256"].get(relative) != digest(path):
            raise ValueError(f"Compiled registry object changed: {path}")
    output = base / "registered"
    command = [sys.executable, str(FITTING / "validate_registered_excitation_atom.py"),
               "--source-dir", str(base / "atom"), "--element", "si", "--stage", "1",
               "--driver", str(REGISTRY_DRIVER), "--work-dir", str(output),
               "--timeout", "1500"]
    with (base / "registered.log").open("w") as stream:
        subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                       check=True, timeout=1560)
    data = json.loads((output / "result.json").read_text())
    errors = [data.get(name) for name in (
        "worst_grid_relative_error", "worst_middle_relative_error", "worst_apply_relative_error")]
    if (data.get("passed") is not True or data.get("samples") != 40
        or data.get("table_cells") != 64561
        or tuple(data.get(name) for name in ("nlev", "nrad", "nups")) != COUNTS[:3]
        or data.get("source_sha256") != atom["atom_sha256"]
        or any(not isinstance(value, (float, int)) or not math.isfinite(value)
               or value < 0 for value in errors)
        or errors[0] >= 1e-10 or max(errors[1:]) >= REGISTRY_TARGET):
        raise ValueError("Si I density table or application exceeds 0.3%")
    result_file(base / "result.json", dict(
        passed=True, ion="SiI", source_atom=atom, fit=fit, registered=data,
        registry_build_manifest_sha256=digest(build_manifest),
        registry_linked_sha256={str(path): digest(path) for path in linked},
        source_accuracy_verified=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("atom", "fit", "registry"))
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--candidate", type=Path)
    args = parser.parse_args()
    base = args.base.resolve()
    base.mkdir(parents=True, exist_ok=True)
    if args.stage == "atom":
        if args.candidate is None:
            parser.error("Corrected source atom path is required")
        atom_stage(base, args.candidate.resolve())
    else:
        {"fit": fit_stage, "registry": registry_stage}[args.stage](base)


if __name__ == "__main__":
    main()
