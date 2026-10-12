"""Focused checks of candidate Fe recombination data and unit conversions."""

import importlib.util
import math
from pathlib import Path
import unittest


SCRIPT = Path(__file__).resolve().parents[2] / "tools/fitting/fe_recombination_expansion.py"
SPEC = importlib.util.spec_from_file_location("fe_recombination_expansion", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FeRecombinationExpansionTest(unittest.TestCase):
    def test_stage_mapping_and_chianti_fit_coefficients(self):
        expected = {
            5: (7.8e-12, 0.682, (4.3285e5, 7.8214e5), (3.8023e-2, 1.6010e-2)),
            6: (1.51e-11, 0.699, (6.2896e5, 1.1605e6), (8.0047e-2, 2.4014e-2)),
            7: (2.62e-11, 0.728, (5.2800e5, 4.1776e6), (9.2054e-2, 4.1024e-2)),
        }
        for stage, (a, eta, energy, coef) in expected.items():
            with self.subTest(stage=stage):
                row = MODULE.read_chianti_fe(MODULE.DEFAULT_CHIANTI, stage)
                self.assertEqual(row["transition"], f"Fe {stage} -> Fe {stage - 1}")
                self.assertEqual((row["rr_A"], row["rr_eta"]), (a, eta))
                self.assertEqual(row["dr_E_K"][:2], energy)
                self.assertEqual(row["dr_c_cm3_s_K3_2"][:2], coef)
                self.assertTrue(all(v == 0 for v in row["dr_E_K"][2:]))
                self.assertTrue(all(v == 0 for v in row["dr_c_cm3_s_K3_2"][2:]))
                for temperature in (1e3, 1e4, 1e5, 1e6):
                    rr, dr = MODULE.chianti_rates(row, temperature)
                    expected_rr = a * (temperature / 1e4) ** -eta
                    expected_dr = sum(c * math.exp(-e / temperature)
                                      for c, e in zip(coef, energy)) / temperature ** 1.5
                    self.assertAlmostEqual(rr, expected_rr, delta=expected_rr * 1e-14)
                    self.assertAlmostEqual(dr, expected_dr, delta=max(expected_dr * 1e-14, 1e-99))

    def test_cloudy_gu_starts_at_fe_vii_and_uses_ev(self):
        row = MODULE.read_cloudy_gu_fe_vii(MODULE.DEFAULT_CLOUDY)
        self.assertEqual(row["stage"], 7)
        self.assertEqual(row["c_cm3_s_eV3_2"][:2], (2.50507e-11, 5.60226e-11))
        self.assertEqual(row["E_eV"][:2], (8.30501e-2, 8.52897e-1))
        self.assertIn("private communication", row["provenance"])
        temperature = 1e4
        temperature_eV = temperature / 11604.518121550082
        expected = sum(c * math.exp(-e / temperature_eV)
                       for c, e in zip(row["c_cm3_s_eV3_2"], row["E_eV"])) / temperature_eV ** 1.5
        self.assertAlmostEqual(MODULE.cloudy_gu_fe_vii_dr(row, temperature), expected,
                               delta=expected * 1e-14)

    def test_invalid_temperature_rejected(self):
        chianti = MODULE.read_chianti_fe(MODULE.DEFAULT_CHIANTI, 5)
        cloudy = MODULE.read_cloudy_gu_fe_vii(MODULE.DEFAULT_CLOUDY)
        for value in (-1.0, 0.0, float("nan"), float("inf")):
            with self.subTest(temperature=value):
                with self.assertRaises(ValueError):
                    MODULE.chianti_rates(chianti, value)
                with self.assertRaises(ValueError):
                    MODULE.cloudy_gu_fe_vii_dr(cloudy, value)


if __name__ == "__main__":
    unittest.main()
