"""Direct checks of selected published charge-transfer rows and indexing."""

import importlib.util
import math
from pathlib import Path
import sys
import tempfile
import unittest


SOURCE = Path(__file__).resolve().parents[2] / "tools/fitting/charge_exchange_expansion.py"
spec = importlib.util.spec_from_file_location("charge_exchange_expansion", SOURCE)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


class ChargeExchangeExpansionTests(unittest.TestCase):
    def test_published_rates_and_charge_index(self):
        rows = module.cloudy_recombination_rows()
        verified = module.published_added_rates()
        self.assertEqual(len(verified), 10)
        magnesium = verified['mg', 4]
        self.assertEqual(magnesium.coefficients, (6.36, .55, 3.86, -5.19))
        self.assertEqual((magnesium.min_temperature, magnesium.max_temperature), (1000., 30000.))
        self.assertEqual(magnesium.cloudy_energy_defect_ev, 8.60)
        self.assertNotEqual(rows["n", 4].coefficients,
                            (3.28, 0.52, -0.52, -0.19))
        self.assertNotIn(("n", 4), verified)
        self.assertEqual(verified["c", 4].coefficients[0], 332.46)
        self.assertEqual(verified["c", 4].min_temperature, 10)
        self.assertEqual(verified["c", 4].max_temperature, 100000)
        self.assertTrue(math.isclose(
            verified["c", 4].rate(10000),
            332.46e-9 * (1 - .995 * math.exp(-.00158)), rel_tol=1e-12))
        self.assertEqual(verified["c", 4].rate(1), verified["c", 4].rate(10))
        self.assertEqual(verified["c", 4].rate(1e6),
                         verified["c", 4].rate(1e5))
        self.assertEqual(verified["n", 3].cloudy_energy_defect_ev, 6.4)


    def test_reject_invalid_temperature_and_source(self):
        row = module.published_added_rates()["o", 3]
        for value in (0, -1, float("nan"), float("inf")):
            with self.subTest(temperature=value):
                with self.assertRaises(ValueError):
                    row.rate(value)
        with tempfile.TemporaryDirectory() as directory:
            bad = Path(directory) / "ctrecombdata.dat"
            bad.write_text("201903042\n1 2 3\n")
            with self.assertRaises(ValueError):
                module.cloudy_recombination_rows(bad)

    def test_existing_kf96_intervals_and_silicon_revision(self):
        published = module.published_existing_rates()
        cloudy = module.cloudy_recombination_rows()
        self.assertEqual(len(published), 13)
        # Si2+ + H is the updated KF96 database (Ferland et al. 2013, Cloudy
        # ctrecombdata.dat), 500-1e5 K, in place of the printed KF96 fit
        # (1.23, 0.24, 3.17, 0.00418; 10-1e6 K) (author's decision,
        # 2026-10-11).
        self.assertEqual((published["si", 2].min_temperature,
                          published["si", 2].max_temperature), (500., 1e5))
        self.assertEqual(published["si", 2].coefficients,
                         (6.77, 0.0736, -0.43, -0.11))
        self.assertEqual(published["si", 2].coefficients,
                         cloudy["si", 2].coefficients)
        self.assertEqual((published["c", 2].min_temperature,
                          published["c", 2].max_temperature),
                         (5000., 50000.))
        self.assertEqual((published["fe", 3].min_temperature,
                          published["fe", 3].max_temperature),
                         (1000., 100000.))


if __name__ == "__main__":
    unittest.main()
