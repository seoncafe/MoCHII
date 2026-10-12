"""Stage one Mg III product check after the shared source is approved."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
FIT_CODE = ROOT / "tools/fitting"
sys.path.insert(0, str(FIT_CODE))
import chianti_cooling as C  # noqa: E402
import fit_cooling as F  # noqa: E402
from validate_stage_cooling import coefficients  # noqa: E402
from write_excitation_atom import write_excitation_atom  # noqa: E402


COOLING_BUILD = HERE.parent / "bt_floor_cooling_products"
COOLING_DRIVER = COOLING_BUILD / "build_mpi/check_cooling_fit.x"
REGISTRY_DRIVER = (HERE.parent / "oxygen_valence_convergence" /
                   "optimized_build/registered_check.x")
EXPECTED = (167, 3364, 13861, 207915)
STALE_CANDIDATE_SHA256 = "46e71f59c1c885a46ec4304f4a26bc1549e716a5490c6118d4cd242b16711873"
THEORETICAL_LEVEL_167_CM1 = 634430.0
MODEL_NOTES = (
    "Mg III full 167-state source; first 27 observed-state spontaneous A values are unchanged.",
    "Higher E1 rates preserve source dipole strength: A scales as photon frequency cubed times source upper-g / physical upper-g.",
    "The 112 higher gf=0 rates retain source A as an unidentified-multipole effective-rate approximation.",
    "The adopted Upsilon is constant below the published 1800 K source minimum; Maxwellian rates use actual gas T.",
    "Higher-level radiative classification and source accuracy remain unverified; numerical validation does not establish physical acceptance.",
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def numerical_records(path):
    """Return exact serialized physical records while ignoring provenance text."""
    return b"\n".join(line for line in Path(path).read_bytes().splitlines()
                      if line.strip() and not line.lstrip().startswith(b"#")) + b"\n"


def level_167_energy(path):
    lines = numerical_records(path).decode("ascii").splitlines()
    if lines[0] != "NLEVEL_VERSION 2" or lines[1] != "NLEV 167":
        raise ValueError("Unexpected Mg III level header")
    tokens = lines[2 + 166].split()
    if len(tokens) != 3 or tokens[0] != "167":
        raise ValueError("Missing source level 167")
    return float(tokens[2])


def write_result(path, result):
    if path.exists():
        raise FileExistsError(f"Preserve existing product: {path}")
    path.write_text(json.dumps(result, indent=2) + "\n")


def annotate_product(path):
    body = path.read_bytes()
    notes = b"".join(("# " + note + "\n").encode("ascii") for note in MODEL_NOTES)
    path.write_bytes(notes + body)


def generate_atom(base, candidate):
    if digest(candidate) == STALE_CANDIDATE_SHA256:
        raise ValueError("Former level-167 Mg III candidate is physically stale")
    if level_167_energy(candidate) != THEORETICAL_LEVEL_167_CM1:
        raise ValueError("Corrected source level 167 must use its theoretical energy")
    output = base / "atom"
    if output.exists():
        raise FileExistsError(f"Preserve existing atom: {output}")
    output.mkdir(parents=True)
    metadata = write_excitation_atom("mg", 3, output_dir=output)
    atom = Path(metadata["path"])
    annotate_product(atom)
    counts = tuple(metadata[key] for key in ("nlev", "nrad", "nups", "nodes"))
    if counts != EXPECTED or metadata["bt_tmin_K"] != 1800.0:
        raise ValueError(f"Unexpected Mg III source atom: {counts}, {metadata['bt_tmin_K']}")
    actual = numerical_records(atom)
    expected = numerical_records(candidate)
    if actual != expected:
        raise ValueError("Generated atom's numerical records differ from validated candidate")
    result = dict(passed=True, stage="atom", counts=counts,
                  bt_tmin_K=1800.0, numerical_records_sha256=hashlib.sha256(actual).hexdigest(),
                  atom_sha256=digest(atom), candidate_sha256=digest(candidate),
                  level_167_energy_cm1=level_167_energy(atom),
                  selected_ids=metadata["selected_ids"],
                  density_ratio_cooling=metadata["density_ratio_cooling"])
    if not result["density_ratio_cooling"] or len(result["selected_ids"]) != 167:
        raise ValueError("Generated atom is not the complete bound source")
    write_result(base / "atom_result.json", result)


def fit_and_validate(base):
    atom = json.loads((base / "atom_result.json").read_text())
    if not atom.get("passed") or atom["counts"] != list(EXPECTED):
        raise ValueError("Source atom was not validated")
    output = base / "fit"
    if output.exists():
        raise FileExistsError(f"Preserve existing fit: {output}")
    output.mkdir()
    command = [sys.executable, str(FIT_CODE / "fit_cooling.py"),
               "--ion-spec", "mg:3", "--output-dir", str(output), "--no-plot"]
    with (output / "fit.log").open("w") as stream:
        subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                       check=True, timeout=3000)
    coefficients_path = output / "cooling_mg_3.txt"
    annotate_product(coefficients_path)
    amplitude, energy = coefficients(coefficients_path)
    temperatures = np.geomspace(1000.0, 100000.0, 801)
    source = C.cooling_low_density_limit("mg", 3, temperatures)
    fitted = F.evaluate(amplitude, energy, temperatures)
    if (not np.all(np.isfinite(source)) or np.any(source < 0)
        or not np.all(np.isfinite(fitted)) or np.any(fitted <= 0)):
        raise ValueError("Nonfinite or nonpositive dense cooling profile")
    if np.any(source == 0):
        log_source = C.log_cooling_ground_state_limit("mg", 3, temperatures)
        log_fit = F.log_cooling_fit(amplitude, energy, temperatures)
        difference = np.expm1(log_fit - log_source)
    else:
        difference = fitted / source - 1.0
    error = float(np.max(np.abs(difference)))
    if not np.all(np.isfinite(difference)) or error >= 0.01:
        raise ValueError(f"Dense source fit exceeds 1%: {error}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(6.5, 6.5),
                             gridspec_kw={"height_ratios": [3, 1]})
    axes[0].loglog(temperatures, source, label="Full 167-state source")
    axes[0].loglog(temperatures, fitted, "--", label="Serialized fit")
    axes[0].set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
    axes[0].legend()
    axes[1].semilogx(temperatures, 100.0 * difference)
    axes[1].set_xlabel(r"$T$ [K]")
    axes[1].set_ylabel(r"Fit difference [\%]")
    fig.tight_layout()
    fig.savefig(output / "MgIII_fit.pdf")
    plt.close(fig)
    build = json.loads((COOLING_BUILD / "build_result.json").read_text())
    if (build.get("passed") is not True
        or build["executables"].get("check_cooling_fit.x") != digest(COOLING_DRIVER)
        or any(build["source_sha256"].get(str((ROOT / "src" / name).resolve()))
               != digest(ROOT / "src" / name)
               for name in ("define.f90", "cooling_mod.f90"))):
        raise ValueError("Cooling reader build does not match its linked source record")
    linked_objects = [COOLING_BUILD / "build_mpi" / name for name in (
        "define.o", "utility.o", "hydrogenic_recombination_mod.o",
        "recomb_mod.o", "gaunt_vh14_mod.o", "nlevel_mod.o",
        "nlevel_cooling_mod.o", "cooling_mod.o", "gaunt.o")]
    linked_object_hashes = {str(path): digest(path) for path in linked_objects}
    native_dir = base / "native_reader"
    native_command = [sys.executable, str(FIT_CODE / "validate_native_cooling_coefficients.py"),
                      "--coefficients", str(coefficients_path), "--driver",
                      str(COOLING_DRIVER), "--work-dir", str(native_dir)]
    with (output / "native_reader.log").open("w") as stream:
        subprocess.run(native_command, cwd=ROOT, stdout=stream,
                       stderr=subprocess.STDOUT, check=True, timeout=180)
    native = json.loads((native_dir / "result.json").read_text())
    if (not native.get("passed") or native["points"] != 801
        or not math.isfinite(native["maximum_resolved_relative_error"])
        or native["maximum_resolved_relative_error"] >= 1e-9):
        raise ValueError("Actual compiled cooling reader failed")
    result = dict(passed=True, stage="fit", temperature_points=801,
                  temperature_K=[1000.0, 100000.0], max_relative_fit_error=error,
                  coefficients_sha256=digest(coefficients_path),
                  fit_figure_sha256=digest(output / "MgIII_fit.pdf"),
                  native_reader=native,
                  cooling_reader_build_record_sha256=digest(COOLING_BUILD / "build_result.json"),
                  cooling_reader_linked_object_sha256=linked_object_hashes,
                  cooling_reader_source_sha256=build["source_sha256"])
    write_result(base / "fit_result.json", result)


def validate_registry(base):
    atom = json.loads((base / "atom_result.json").read_text())
    fit = json.loads((base / "fit_result.json").read_text())
    if not atom.get("passed") or not fit.get("passed"):
        raise ValueError("Earlier product checks are incomplete")
    registry_build_path = HERE.parent / "oxygen_valence_convergence/optimized_build_manifest.json"
    registry_build = json.loads(registry_build_path.read_text())
    registry_objects = [REGISTRY_DRIVER.parent / name for name in (
        "define.o", "nlevel_mod.o", "nlevel_cooling_mod.o")]
    registry_hashes = {str(path.relative_to(ROOT)): digest(path)
                       for path in [REGISTRY_DRIVER] + registry_objects}
    if any(registry_build["sha256"].get(path) != sha
           for path, sha in registry_hashes.items()):
        raise ValueError("Optimized registry driver differs from isolated build record")
    output = base / "registered"
    if output.exists():
        raise FileExistsError(f"Preserve existing registry: {output}")
    command = [sys.executable, str(FIT_CODE / "validate_registered_excitation_atom.py"),
               "--source-dir", str(base / "atom"), "--element", "mg", "--stage", "3",
               "--driver", str(REGISTRY_DRIVER), "--work-dir", str(output),
               "--timeout", "1500"]
    with (base / "registered.log").open("w") as stream:
        subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                       check=True, timeout=1560)
    data = json.loads((output / "result.json").read_text())
    if (not data.get("passed") or data.get("samples") != 40
        or data.get("table_cells") != 64561
        or tuple(int(data[key]) for key in ("nlev", "nrad", "nups")) != EXPECTED[:3]
        or data["source_sha256"] != atom["atom_sha256"]
        or any(not math.isfinite(data[name]) for name in (
            "worst_grid_relative_error", "worst_middle_relative_error",
            "worst_apply_relative_error"))
        or data["worst_grid_relative_error"] >= 1e-10
        or data["worst_middle_relative_error"] >= 0.005
        or data["worst_apply_relative_error"] >= 0.005):
        raise ValueError("Registered density table or application fails acceptance")
    write_result(base / "result.json", dict(
        passed=True, ion="MgIII", source_atom=atom, fit=fit, registered=data,
        registry_build_manifest_sha256=digest(registry_build_path),
        registry_linked_object_and_driver_sha256=registry_hashes,
        note=("40 native checks include temperatures outside the published "
              "1,800 K minimum; those low-temperature checks establish only "
              "the adopted constant-Upsilon numerical continuation.")))


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
            parser.error("The validated corrected candidate path is required")
        generate_atom(base, args.candidate.resolve())
    else:
        {"fit": fit_and_validate, "registry": validate_registry}[args.stage](base)


if __name__ == "__main__":
    main()
