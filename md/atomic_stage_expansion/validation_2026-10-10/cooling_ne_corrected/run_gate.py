#!/usr/bin/env python3
"""Run the isolated cooling gate with finite resources and preserve its output."""
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def limit_memory():
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))


def main():
    root = Path(__file__).resolve().parents[4]
    work = Path(__file__).resolve().parent
    exe = work / "isolated_build/test_cool_table"
    atomic = root / "data/atomic"
    names = (("o", 3), ("o", 2), ("n", 2), ("s", 2), ("ne", 3),
             ("c", 2), ("fe", 3), ("c", 1), ("n", 1))
    inputs = [exe, root / "src/define.f90", root / "src/nlevel_mod.f90",
              root / "src/nlevel_cooling_mod.f90",
              root / "tests/cooling_ne/test_cool_table.f90",
              root / "tests/cooling_ne/build_test.sh"]
    for name, stage in names:
        inputs.extend((atomic / f"nlevel_{name}_{stage}.txt",
                       atomic / f"cooling_{name}_{stage}.txt"))
    before = {str(path): digest(path) for path in inputs}
    try:
        run = subprocess.run([str(exe), str(atomic)], cwd=work, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             timeout=1600, preexec_fn=limit_memory, check=False)
        output = run.stdout
        code = run.returncode
    except subprocess.TimeoutExpired as error:
        output = error.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        code = 124
    (work / "test_cool_table.log").write_text(output)
    after = {str(path): digest(path) for path in inputs}
    gates = {gate: f"{gate} PASS" in output for gate in ("G-a", "G-c", "G-e")}
    result = {"passed": code == 0 and all(gates.values()) and before == after,
              "return_code": code, "gates": gates,
              "sources_unchanged": before == after, "sha256": before,
              "log": str(work / "test_cool_table.log")}
    (work / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "sha256"}))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
