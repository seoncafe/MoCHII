"""Build isolated cooling readers from current sources."""

import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BUILD = HERE / "build_mpi"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    BUILD.mkdir(exist_ok=False)
    sources = [ROOT / name for name in (
        "src/define.f90", "src/utility.f90",
        "src/hydrogenic_recombination_mod.f90", "src/recomb_mod.f90",
        "src/gaunt_vh14_mod.f90", "src/nlevel_mod.f90",
        "src/nlevel_cooling_mod.f90", "src/cooling_mod.f90", "src/gaunt.f90",
        "tests/ion_stages/check_cooling_fit.f90",
        "tests/ion_stages/check_registered_excitation_atom.f90")]
    before = {str(path): digest(path) for path in sources}
    flags = ["mpifort", "-O2", "-cpp", "-DMPI", "-fcheck=all",
             "-ffpe-trap=invalid,zero,overflow", "-ffree-line-length-none",
             "-J" + str(BUILD), "-I" + str(BUILD),
             "-I" + str(ROOT / "src")]
    commands = []
    for path in sources[:9]:
        command = flags + ["-c", str(path)]
        subprocess.run(command, cwd=BUILD, check=True, timeout=120)
        commands.append(command)
    for driver, objects, output in (
        (sources[9], ["define.o", "utility.o", "hydrogenic_recombination_mod.o",
                      "recomb_mod.o", "gaunt_vh14_mod.o", "nlevel_mod.o",
                      "nlevel_cooling_mod.o", "cooling_mod.o", "gaunt.o"],
         "check_cooling_fit.x"),
        (sources[10], ["define.o", "nlevel_mod.o", "nlevel_cooling_mod.o"],
         "check_registered.x")):
        command = flags + [str(driver), *objects, "-o", output]
        subprocess.run(command, cwd=BUILD, check=True, timeout=120)
        commands.append(command)
    after = {str(path): digest(path) for path in sources}
    if before != after:
        raise RuntimeError("Source changed during isolated build")
    result = {"passed": True, "source_sha256": before, "commands": commands,
              "executables": {name: digest(BUILD / name) for name in
                              ("check_cooling_fit.x", "check_registered.x")}}
    (HERE / "build_result.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
