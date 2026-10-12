"""Run the public S I cooling generator and record a complete staged product."""
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "tools/fitting"))
import fit_cooling as F


def main():
    output = HERE / "normal_generator"
    output.mkdir(exist_ok=True)
    F.generate({"SI": ("s", 1)}, output_dir=output,
               plot_path=output / "cooling_fits.png")
    product = output / "cooling_s_1.txt"
    lines = [line for line in product.read_text().splitlines()
             if line.strip() and not line.startswith("#")]
    count = int(lines[0])
    if len(lines) != count + 1 or count < 1:
        raise ValueError("S I cooling product is incomplete")
    result = {"generated": True, "term_count": count,
              "product_sha256": hashlib.sha256(product.read_bytes()).hexdigest()}
    (output / "generated.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
