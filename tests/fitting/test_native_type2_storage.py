"""Native-coordinate fits retain the adopted collision-strength source."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools/fitting"))
from chianti_cooling import bound_excitation_atom, upsilon
from fit_nlevel import build_ion, eval_ups_fit, fit_transition


class NativeType2StorageTest(unittest.TestCase):
    def test_carbon_positive_interpolation_survives_native_rounding(self):
        levels, _, transitions = bound_excitation_atom("c", 1)
        self.assertGreaterEqual(len(levels), 5)
        temperature = np.logspace(3, 5, 160001)
        for pair in ((1, 2), (1, 3), (2, 3)):
            with self.subTest(pair=pair):
                transition = next(row for row in transitions
                                  if (row["ll"], row["ul"]) == pair)
                original_constant = transition["cups"]
                original_type = transition["ttype"]
                source_before = upsilon(transition, temperature)
                lower, upper, coefficients, logarithmic, error = fit_transition(transition)
                source_after = upsilon(transition, temperature)
                np.testing.assert_array_equal(source_before, source_after)
                self.assertEqual(transition["cups"], original_constant)
                self.assertEqual(transition["ttype"], original_type)
                self.assertEqual(transition["native_ttype"], 2)
                self.assertNotEqual(transition["native_cups"], original_constant)
                self.assertEqual(len(coefficients), 16)
                fitted = eval_ups_fit(transition, lower, upper, coefficients,
                                      logarithmic, temperature)
                self.assertLess(np.max(np.abs(fitted/source_before-1)), 0.002)
                self.assertLess(error, 0.002)

    def test_staged_rows_preserve_source_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            path, levels, _, _, worst = build_ion(
                "CI", "c", 1, 5, output_dir=directory)
            self.assertEqual(levels, 5)
            self.assertLess(worst, 0.002)
            lines = Path(path).read_text().splitlines()
            metadata = [line for line in lines if line.startswith("# native ")]
            self.assertEqual(len(metadata), 3)
            self.assertTrue(any("source type 2 C 6.247000e+02" in line
                                for line in metadata))
            position = next(i for i, line in enumerate(lines)
                            if line.startswith("NUPS "))
            count = int(lines[position].split()[1])
            rows = {tuple(map(int, fields[:2])): fields
                    for fields in (line.split() for line in lines[position+1:position+1+count])}
            for pair in ((1, 2), (1, 3), (2, 3)):
                fields = rows[pair]
                self.assertEqual(int(fields[2]), 2)
                self.assertEqual(int(fields[8]), 16)
                self.assertNotEqual(float(fields[3]),
                                    {(1,2):624.7,(1,3):208.2,(2,3):312.4}[pair])


if __name__ == "__main__":
    unittest.main()
