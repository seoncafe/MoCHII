"""Run the linked Fortran cooling reader on the staged S I coefficients."""
from pathlib import Path
import json
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import chianti_cooling as source
from validate_added_charge_transfer import compile_physics_check

HERE = Path(__file__).resolve().parent


def main():
    work = HERE / "native_check"
    work.mkdir(exist_ok=True)
    executable = compile_physics_check("tests/ion_stages/check_cooling_fit.f90", work)
    output = work / "si_native_801.txt"
    subprocess.run([str(executable), str(HERE / "cooling_s_1.txt"), str(output)],
                   check=True)
    rows = np.loadtxt(output)
    if rows.shape != (801, 2):
        raise ValueError("Incomplete native S I cooling profile")
    temperature, measured = rows[:,0], rows[:,1]
    reference = source.cooling_low_density_limit("s", 1, temperature)
    difference = np.abs(measured/reference-1)
    result = {"points":len(temperature),
              "temperature_K":[float(temperature[0]),float(temperature[-1])],
              "maximum_absolute_relative_error":float(difference.max()),
              "worst_temperature_K":float(temperature[difference.argmax()]),
              "all_finite_positive":bool(np.all(np.isfinite(measured)) and np.all(measured>0))}
    (HERE/"native_result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
