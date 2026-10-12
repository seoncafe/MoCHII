"""Focused checks for CHIANTI energy-policy diagnostics."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "fitting"))

from atomic_level_energies import audit_ion, audit_transitions, select_level_energies


class AtomicLevelEnergiesTest(unittest.TestCase):
    def test_missing_observation_does_not_become_zero_energy(self):
        levels = {1: (1, 0, 0), 2: (3, -1, 100), 3: (1, 80, 120)}
        self.assertEqual(select_level_energies(levels, "mixed"),
                         {1: 0, 2: 100, 3: 80})
        self.assertEqual(select_level_energies(levels, "observed"),
                         {1: 0, 2: None, 3: 80})

    def test_reversed_and_equal_rows_remain_visible(self):
        levels = {1: (1, 0, 0), 2: (3, 100, 100), 3: (1, 80, 100)}
        radiation = [(2, 3, 5000, 0, 1)]
        collisions = [{"ll": 2, "ul": 3, "de": 0.01}]
        mixed = audit_transitions(levels, radiation, collisions, "mixed")
        theoretical = audit_transitions(levels, radiation, collisions, "theoretical")
        self.assertEqual(mixed["radiative_reversed"][0]["gap_cm1"], -20)
        self.assertEqual(mixed["collision_reversed"][0]["pair"], [2, 3])
        self.assertEqual(theoretical["radiative_equal"][0]["gap_cm1"], 0)
        self.assertEqual(theoretical["collision_equal"][0]["pair"], [2, 3])

    def test_n_iv_observation_reverses_real_radiative_row(self):
        record = audit_ion("n", 4)
        mixed = record["policies"]["mixed"]
        theoretical = record["policies"]["theoretical"]
        self.assertIn([107, 116], [row["pair"] for row in mixed["radiative_reversed"]])
        self.assertNotIn([107, 116], [row["pair"] for row in theoretical["radiative_reversed"]])


if __name__ == "__main__":
    unittest.main()
