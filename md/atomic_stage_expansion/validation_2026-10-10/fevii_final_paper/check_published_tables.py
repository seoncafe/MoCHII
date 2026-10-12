"""Compare final journal Tables 3--5 with the unchanged Stout Fe VII files."""
from pathlib import Path
import hashlib
import json
import math
import re
import subprocess

HERE = Path(__file__).resolve().parent
PDF = Path.home() / "MoCHII/references/Witthoeft_Badnell_2008_AA_481_543.pdf"
STEM = Path.home() / "CLOUDY/c25.00/data/stout/fe/fe_7/fe_7"
RY_CM1 = 109737.31568539


def scientific_number(value):
    match = re.fullmatch(r"([0-9.]+)([+\u2212]\d+)?", value)
    if not match:
        raise ValueError(value)
    return float(match[1]) * 10 ** int((match[2] or "0").replace("\u2212", "-"))


def main():
    raw = subprocess.check_output(["pdftotext", "-raw", str(PDF), "-"], text=True)
    records = re.findall(r"^(\d+) 3p[56]\n[^\n]+\n[A-Z]\d+ ([0-9.]+)(?: ([0-9.]+))?", raw, re.M)
    if len(records) != 189 or {int(r[0]) for r in records} != set(range(1, 190)):
        raise ValueError("Incomplete final Table 3")
    levels = {int(f[0]): float(f[1]) for line in Path(str(STEM)+".nrg").read_text().splitlines()[1:]
              if line and line[0].isdigit() for f in [line.split()]}
    rates = {(int(f[1]), int(f[2])): float(f[3])
             for line in Path(str(STEM)+".tp").read_text().splitlines()
             if line.startswith("A\t") for f in [line.split()]}
    rows = []
    for index, present, observed in records:
        index = int(index)
        rows.append(dict(level=index, present_Ry=float(present),
                         NIST_Ry=float(observed) if observed else None,
                         Stout_cm1=levels[index], Stout_Ry=levels[index]/RY_CM1,
                         difference_from_present_Ry=levels[index]/RY_CM1-float(present)))
    radiation = []
    table4 = raw.split("i j Present N&S Berr.", 1)[1].split("In Fe6+", 1)[0]
    for lo, hi, value in re.findall(r"^(\d+) (\d+) ([0-9.]+(?:[+\u2212]\d+)?)", table4, re.M):
        pair = int(lo), int(hi)
        radiation.append(dict(table=4, lower=pair[0], upper=pair[1],
                              published=scientific_number(value), Stout=rates[pair]))
    table5 = raw.split("i j 36-config. 40-config. W&K", 1)[1].split("(10%) with Warner", 1)[0]
    for lo, hi, old, value in re.findall(r"^(\d+) (\d+) ([0-9.]+[+\u2212]\d+)\n([0-9.]+[+\u2212]\d+)", table5, re.M):
        pair = int(lo), int(hi)
        radiation.append(dict(table=5, lower=pair[0], upper=pair[1],
                              published=scientific_number(value), Stout=rates[pair]))
    if (sum(r["table"] == 4 for r in radiation) != 27
            or sum(r["table"] == 5 for r in radiation) != 30
            or any(not math.isclose(r["published"], r["Stout"], rel_tol=1e-13)
                   for r in radiation)):
        raise ValueError("Incomplete or mismatched radiative tables")
    result = dict(pdf_sha256=hashlib.sha256(PDF.read_bytes()).hexdigest(),
                  conversion_cm1_per_Ry=RY_CM1, paper_rows=len(rows), rows=rows,
                  maximum_difference_levels_10_189_Ry=max(abs(r["difference_from_present_Ry"]) for r in rows if r["level"]>9),
                  maximum_difference_levels_1_9_NIST_Ry=max(abs(r["Stout_Ry"]-r["NIST_Ry"]) for r in rows if r["level"]<=9),
                  radiative_rows=radiation)
    (HERE/"level_correspondence.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("rows", "radiative_rows")}))
    print("Published radiative rows:", len(radiation))


if __name__ == "__main__":
    main()
