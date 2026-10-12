"""Focused checks for selective atomic-data generation safety."""

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
FITTING = ROOT / "tools" / "fitting"
sys.path.insert(0, str(FITTING))

import atomic_output
import make_element_data


class GenerationSafety(unittest.TestCase):
    def test_imports_do_not_write(self):
        with tempfile.TemporaryDirectory() as directory:
            code = (
                "import sys; sys.path.insert(0, sys.argv[1]); "
                "import make_element_data, fit_cooling, fit_nlevel"
            )
            subprocess.run([sys.executable, "-c", code, str(FITTING)],
                           cwd=directory, check=True, capture_output=True)
            self.assertEqual(os.listdir(directory), [])

    def test_missing_threshold_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "element_c.txt"
            output.write_text("unchanged\n")
            with patch.dict(make_element_data.ETH, {}, clear=True):
                with self.assertRaisesRegex(ValueError, "Missing selected atomic sources.*ETH"):
                    make_element_data.generate([("c", 6, 4)], directory)
            self.assertEqual(output.read_text(), "unchanged\n")
            self.assertEqual(os.listdir(directory), ["element_c.txt"])

    def test_failed_replace_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "rate.txt"
            output.write_text("unchanged\n")
            with patch.object(atomic_output.os, "replace", side_effect=OSError("test failure")):
                with self.assertRaises(OSError):
                    atomic_output.atomic_write_text(output, "new\n")
            self.assertEqual(output.read_text(), "unchanged\n")
            self.assertEqual(os.listdir(directory), ["rate.txt"])

    def test_default_scope_is_stages_one_to_seven(self):
        # Stages I-VII for every element (make_element_data.ELEMENTS; 77
        # stages, 66 transitions).
        self.assertEqual(len(make_element_data.ELEMENTS), 11)
        self.assertTrue(all(n == 7 for _, _, n in make_element_data.ELEMENTS))

    def test_default_nitrogen_charge_exchange_sources(self):
        # The default nitrogen scope generates without refusal: N4+ + H from
        # the updated KF96 database (Ferland et al. 2013, Cloudy
        # ctrecombdata.dat), N5+ and N6+ + H from the Langevin capture limit.
        with tempfile.TemporaryDirectory() as directory:
            make_element_data.generate([("n", 7, 7)], directory)
            self.assertEqual(os.listdir(directory), ["element_n.txt"])
            blocks = (Path(directory) / "element_n.txt").read_text().split("TRANSITION ")[1:]
        self.assertEqual([int(b.split()[0]) for b in blocks], [1, 2, 3, 4, 5, 6])
        self.assertIn("# CXR: updated KF96 database (Ferland et al. 2013", blocks[3])
        self.assertIn("CXR 8 ", blocks[3])
        for block in blocks[4:]:
            self.assertIn("Ferland 1997 Langevin capture upper limit", block)

    def test_default_silicon_charge_exchange_source(self):
        # Si2+ + H -> Si+ from the updated KF96 database (Ferland et al. 2013).
        with tempfile.TemporaryDirectory() as directory:
            make_element_data.generate([("si", 14, 7)], directory)
            blocks = (Path(directory) / "element_si.txt").read_text().split("TRANSITION ")[1:]
        self.assertIn("# CXR: updated KF96 database (Ferland et al. 2013", blocks[1])

    def test_stage_beyond_default_cannot_silently_omit_charge_exchange(self):
        # A stage above the default scope still needs resolved charge
        # exchange and product energies; nothing is written without them.
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "n stage 7: unresolved added-boundary"):
                make_element_data.generate([("n", 7, 8)], directory)
            self.assertEqual(os.listdir(directory), [])

if __name__ == "__main__":
    unittest.main()
