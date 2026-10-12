#!/usr/bin/env python3
"""Measure complete N II and S II source cooling against installed atoms."""

import json
from pathlib import Path
import resource

from check_existing_density import BASE, ROOT, DRIVER, digest, inspect_atom
import chianti_cooling as C


def main():
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))
    ions = (("NII", "n", 2), ("SII", "s", 2))
    files = [Path(__file__), ROOT / "tools/fitting/chianti_cooling.py",
             ROOT / "tools/fitting/validate_excitation_density_profiles.py",
             ROOT / "tools/fitting/atomic_excitation_energies.py",
             ROOT / "tools/fitting/ionization_potentials.py", DRIVER]
    from ionization_potentials import SOURCE
    files.append(SOURCE)
    for _, element, stage in ions:
        stem = Path(C.ion_dir(element, stage)) / f"{element}_{stage}"
        files.extend(Path(f"{stem}.{extension}")
                     for extension in ("elvlc", "wgfa", "scups"))
        files.append(ROOT / f"data/atomic/nlevel_{element}_{stage}.txt")
    before = {str(path): digest(path) for path in files}
    measured = {label: inspect_atom(label, element, stage)
                for label, element, stage in ions}
    after = {str(path): digest(path) for path in files}
    result = {"passed": before == after, "sources_unchanged": before == after,
              "source_sha256": before, "ions": measured}
    (BASE / "n_s_result.json").write_text(json.dumps(result, indent=2) + "\n")
    if before != after:
        raise RuntimeError("Source changed during calculation")
    print(json.dumps({"passed": result["passed"], "ions": list(measured)}))


if __name__ == "__main__":
    main()
