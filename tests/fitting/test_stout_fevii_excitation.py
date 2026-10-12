"""Fe VII Stout source identity, direction, and cold-rate checks."""

from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools/fitting"))

import chianti_cooling as chianti
from stout_excitation_atoms import fevii_excitation_atom


class FeviiExcitation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.levels, cls.radiation, cls.collisions = fevii_excitation_atom()

    def test_complete_bound_atom_and_radiative_directions(self):
        self.assertEqual(len(self.levels), 189)
        self.assertEqual(len(self.collisions), 17766)
        self.assertEqual(len(self.radiation), 17531)
        self.assertEqual(self.levels[1], (5.0, 0.0, 0.0))
        self.assertTrue(all(self.levels[u][1] > self.levels[l][1]
                            for l, u, _, _, rate in self.radiation if rate > 0))
        self.assertTrue(all(self.levels[row["ul"]][1] > self.levels[row["ll"]][1]
                            for row in self.collisions))

    def test_source_temperatures_and_actual_cold_boltzmann_factor(self):
        row = next(row for row in self.collisions
                   if row["ll"] == 1 and row["ul"] == 10)
        grid = row["temperature_grid"]
        self.assertEqual((len(grid), grid[0], grid[-1]), (13, 9800.0, 9.8e7))
        temperature = np.array([1000.0, 3000.0, 9800.0])
        strength = chianti.upsilon(row, temperature)
        np.testing.assert_allclose(strength, row["upsilon_grid"][0], rtol=1e-14)
        gap = (self.levels[10][1] - self.levels[1][1]) * chianti.HC_ERG_CM
        rate = (chianti.COLL_PREF * strength / (self.levels[1][0] * np.sqrt(temperature))
                * np.exp(-gap / (chianti.K_B_ERG * temperature)))
        self.assertTrue(np.all(np.diff(rate) > 0))
        self.assertTrue(np.all(np.isfinite(rate)))


if __name__ == "__main__":
    unittest.main()
