#!/usr/bin/env python3
"""Build the current source-node driver without touching shared objects."""

import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[4]
BUILD = Path(__file__).resolve().parent / "build"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    BUILD.mkdir(exist_ok=False)
    flags = ["mpifort", "-O2", "-fcheck=all",
             "-ffpe-trap=invalid,zero,overflow", "-ffree-line-length-none",
             "-J" + str(BUILD), "-I" + str(BUILD)]
    sources = [ROOT / "src/define.f90", ROOT / "src/nlevel_mod.f90",
               ROOT / "tests/ion_stages/check_source_excitation_atom.f90"]
    before = {str(path): digest(path) for path in sources}
    commands = []
    for path in sources[:2]:
        command = flags + ["-c", str(path)]
        subprocess.run(command, cwd=BUILD, check=True, timeout=120)
        commands.append(command)
    executable = BUILD / "check_source_excitation_atom.x"
    command = flags + [str(sources[2]), "define.o", "nlevel_mod.o",
                       "-o", str(executable)]
    subprocess.run(command, cwd=BUILD, check=True, timeout=120)
    commands.append(command)
    after = {str(path): digest(path) for path in sources}
    result = {"passed": before == after and executable.exists(),
              "source_sha256": before, "executable": str(executable),
              "executable_sha256": digest(executable), "commands": commands}
    (BUILD.parent / "build_result.json").write_text(json.dumps(result, indent=2) + "\n")
    if not result["passed"]:
        raise RuntimeError("Build source changed")


if __name__ == "__main__":
    main()
