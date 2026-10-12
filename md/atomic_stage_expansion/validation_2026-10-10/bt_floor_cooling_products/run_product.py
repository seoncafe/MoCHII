"""Stage one cold-floor cooling fit and validate its native consumers."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PREVIOUS = HERE.parent / "bt_source_temperature_floor_corrected"
IONS = {"CIV": ("c", 4, "CIV"), "SiV": ("si", 5, "SiV")}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ion", choices=IONS, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    element, stage, previous_name = IONS[args.ion]
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError("Preserve previous products; use a new directory")
    output.mkdir(parents=True)
    python = sys.executable
    fit_command = [python, str(ROOT / "tools/fitting/validate_stage_cooling.py"),
                   "--ion", args.ion, "--element", element, "--stage", str(stage),
                   "--output-dir", str(output)]
    with (output / "fit.log").open("w") as log:
        subprocess.run(fit_command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=1500)
    fitted = json.loads((output / f"{args.ion}.json").read_text())
    coefficients = output / f"cooling_{element}_{stage}.txt"
    if not fitted.get("passed") or fitted["max_relative_fit_error"] >= .01:
        raise ValueError("Dense rounded fit failed")
    native_dir = output / "native_reader"
    native_command = [python, str(ROOT / "tools/fitting/validate_native_cooling_coefficients.py"),
                      "--coefficients", str(coefficients), "--driver",
                      str(HERE / "build_mpi/check_cooling_fit.x"),
                      "--work-dir", str(native_dir)]
    with (output / "native_reader.log").open("w") as log:
        subprocess.run(native_command, cwd=ROOT, stdout=log,
                       stderr=subprocess.STDOUT, check=True, timeout=180)
    native = json.loads((native_dir / "result.json").read_text())
    if not native.get("passed") or native["points"] != 801:
        raise ValueError("Native reader failed")
    earlier = json.loads((PREVIOUS / previous_name / "result.json").read_text())
    atom_source = Path(earlier["path"])
    atom_dir = output / "atom"
    atom_dir.mkdir()
    atom = atom_dir / f"nlevel_{element}_{stage}.txt"
    shutil.copyfile(atom_source, atom)
    if digest(atom) != digest(atom_source):
        raise ValueError("Copied source atom differs")
    registry_dir = output / "registered"
    registry_command = [python, str(ROOT / "tools/fitting/validate_registered_excitation_atom.py"),
                        "--source-dir", str(atom_dir), "--element", element,
                        "--stage", str(stage), "--driver",
                        str(HERE / "build_mpi/check_registered.x"),
                        "--work-dir", str(registry_dir), "--timeout", "1200"]
    with (output / "registered_validation.log").open("w") as log:
        subprocess.run(registry_command, cwd=ROOT, stdout=log,
                       stderr=subprocess.STDOUT, check=True, timeout=1260)
    registry = json.loads((registry_dir / "result.json").read_text())
    if not registry.get("passed") or registry["samples"] != 40:
        raise ValueError("Registered table failed")
    result = {
        "passed": True, "ion": args.ion, "physical_scope": (
            "Source atom and cold-constant Upsilon are numerically represented. "
            "Source accuracy below Tmin and above the last source node is not certified. "
            "Registered 100 to 1e6 K samples outside source validity check only software interpolation."),
        "fit": fitted, "native_reader": native, "registered": registry,
        "atom_sha256": digest(atom),
        "outputs_sha256": {str(path): digest(path) for path in (
            coefficients, output / f"{args.ion}.json", output / f"{args.ion}.pdf",
            native_dir / "result.json", native_dir / "native_cooling.txt",
            registry_dir / "result.json", registry_dir / "registered.txt")},
    }
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(output / "result.json")


if __name__ == "__main__":
    main()
