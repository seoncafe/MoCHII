"""Check final-journal Fe recombination identities and source limitations."""

import importlib.util
import tempfile
import unittest
import zipfile
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "tools/fitting/zhang_fe_recombination.py"
SPEC = importlib.util.spec_from_file_location("zhang_fe_recombination", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


class ZhangFeRecombinationTest(unittest.TestCase):
    def test_author_grid_identity_and_known_rates(self):
        expected = {
            3: (1, 90., 1.81e-10, 9000., 1.43e-11),
            4: (4, 160., 3.06e-10, 8000., 2.93e-11),
            5: (4, 250., 2.43e-10, 12500., 7.06e-11),
            6: (3, 360., 3.86e-10, 7200., 7.32e-11),
        }
        for charge, (parents, first_t, first_dr, middle_t, middle_dr) in expected.items():
            with self.subTest(charge=charge):
                grid = AUDIT.read_ground_dr_grid(AUDIT.ARCHIVE, charge)
                self.assertEqual(len(grid), 19)
                self.assertEqual(grid[0], (first_t, first_dr))
                self.assertIn((middle_t, middle_dr), grid)
                self.assertEqual(AUDIT.SEQUENCES[charge][2], parents)

    def test_final_source_and_physical_domain(self):
        AUDIT.verify_final_recombination_source()
        for charge in range(3,7):
            for temperature in (20*charge**2,1e4,1e7*charge**2):
                self.assertGreater(AUDIT.rr_fit(charge,temperature),0)
                self.assertGreater(AUDIT.dr_fit(charge,temperature),0)
            for temperature in (0,-1,float('nan'),float('inf')):
                with self.assertRaises(ValueError):
                    AUDIT.rr_fit(charge,temperature)
                with self.assertRaises(ValueError):
                    AUDIT.dr_fit(charge,temperature)

    def test_fit_disagreement_is_visible(self):
        measured = AUDIT.compare(AUDIT.ARCHIVE)["charge_states"]
        self.assertEqual(measured["4"]["maximum_valid_grid_difference_temperature_K"], 320.)
        self.assertGreater(measured["4"]["maximum_valid_grid_relative_difference"], 0.04)
        self.assertEqual(measured["5"]["maximum_valid_grid_difference_temperature_K"], 500.)
        self.assertGreater(measured["5"]["maximum_valid_grid_relative_difference"], 0.05)
        self.assertEqual(measured["6"]["maximum_valid_grid_difference_temperature_K"], 1800.)
        self.assertGreater(measured["6"]["maximum_valid_grid_relative_difference"], 0.03)

    def test_wrong_parent_identity_is_rejected(self):
        member = "supplement/adf09_Fe4+_S1"
        with zipfile.ZipFile(AUDIT.ARCHIVE) as archive:
            wrong = archive.read(member).replace(b"SEQ='V '", b"SEQ='SC'", 1)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "wrong.zip"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr(member, wrong)
            with self.assertRaisesRegex(ValueError, "identity"):
                AUDIT.read_ground_dr_grid(path, 4)


if __name__ == "__main__":
    unittest.main()
