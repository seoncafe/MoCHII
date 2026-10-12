"""Check the identified Ar II and Ca IV fine-structure source records."""

from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools' / 'fitting'))
from chianti_cooling import upsilon
from stout_fine_structure_atoms import stout_fine_structure_atom


class FineStructureAtomTest(unittest.TestCase):
    def test_published_ground_pair_and_cold_continuation(self):
        for element, stage, gap, rate, first, last in (
            ('ar', 2, 1431.583, .0530, 2.48, 2.71),
            ('ca', 4, 3118.200, .543, 1.00, 1.53),
        ):
            with self.subTest(element=element, stage=stage):
                levels, radiation, collisions = stout_fine_structure_atom(element, stage)
                self.assertEqual(levels[1][0], 4.)
                self.assertEqual(levels[2][0], 2.)
                self.assertEqual(levels[2][1], gap)
                self.assertEqual(radiation[0][4], rate)
                self.assertEqual(len(collisions), 1)
                transition = collisions[0]
                self.assertEqual(len(transition['temperature_grid']), 11)
                self.assertEqual(transition['upsilon_grid'][0], first)
                self.assertEqual(transition['upsilon_grid'][-1], last)
                self.assertEqual(float(upsilon(transition, 1000.)[0]), first)
                with self.assertRaises(ValueError):
                    upsilon(transition, transition['temperature_grid'][-1] * 1.01)

    def test_changed_raw_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stem = root / 'ar' / 'ar_2'
            stem.mkdir(parents=True)
            (stem / 'ar_2.nrg').write_text('17 09 05\n1 0 4\n', encoding='ascii')
            with self.assertRaisesRegex(ValueError, 'Changed Stout source'):
                stout_fine_structure_atom('ar', 2, root)

    def test_table_node_coordinates_are_positive_and_ordered(self):
        for element, stage in (('ar', 2), ('ca', 4)):
            _, _, collisions = stout_fine_structure_atom(element, stage)
            row = collisions[0]
            self.assertTrue(np.all(np.diff(row['temperature_grid']) > 0))
            self.assertTrue(np.all(row['upsilon_grid'] > 0))
            self.assertAlmostEqual(float(upsilon(row, row['temperature_grid'][5])[0]),
                                   float(row['upsilon_grid'][5]), places=12)


if __name__ == '__main__':
    unittest.main()
