"""Focused checks for version-2 excitation source serialization."""

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools' / 'fitting'))
import write_excitation_atom as writer


class SourceStorageTest(unittest.TestCase):
    def setUp(self):
        self.levels = {1: (1., 0.), 4: (3., 100.), 8: (5., 200.)}
        # Two source rows of one pair (e.g. E1 and M2) are summed. Rows are
        # given lower-upper by energy: the writer refuses a reversed row
        # (chianti_cooling.bound_excitation_atom reorients E1 rows first).
        self.radiation = [(1, 4, 0., 0., 2.), (1, 4, 0., 0., 3.)]
        self.collision = [dict(ll=4, ul=1, de=.01, cups=1., ttype=6,
                               xs=np.array([0., .5, 1.]),
                               ys=np.array([-2., -1., 0.]))]

    def test_compact_indices_raw_nodes_and_repeatability(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'source.scups'
            source.write_text('identified source bytes', encoding='ascii')
            with (patch.object(writer, 'bound_excitation_atom',
                               return_value=(self.levels, self.radiation, self.collision)),
                  patch.object(writer, '_source_files', return_value=[source])):
                result = writer.write_excitation_atom('X', 2, 2, directory)
                contents = Path(result['path']).read_bytes()
                writer.write_excitation_atom('X', 2, 2, directory)
                self.assertEqual(contents, Path(result['path']).read_bytes())
            records = [line for line in contents.decode().splitlines()
                       if line and not line.startswith('#')]
            self.assertEqual(records[:4], ['NLEVEL_VERSION 2', 'NLEV 2',
                                            '1 1 0', '2 3 100'])
            self.assertEqual(records[4:8], ['NRAD 1', '1 2 5', 'NODES 3', 'NUPS 1'])
            self.assertEqual(records[8], '1 2 BT_SOURCE 6 1 0.01 CUBIC 3')
            self.assertEqual(records[9:], ['0 -2', '0.5 -1', '1 0'])
            self.assertEqual(result['selected_ids'], [1, 4])
            self.assertEqual(result['nodes'], 3)
            self.assertIn('SHA-256', contents.decode())

    def test_complete_atom_records_density_ratio(self):
        with tempfile.TemporaryDirectory() as directory:
            with (patch.object(writer, 'bound_excitation_atom',
                               return_value=(self.levels, self.radiation, self.collision)),
                  patch.object(writer, '_source_files', return_value=[])):
                result = writer.write_excitation_atom('X', 2, None, directory)
            self.assertTrue(result['density_ratio_cooling'])
            self.assertTrue(Path(result['path']).read_text().endswith('DENSITY_RATIO 1\n'))

    def test_negative_nonlogarithmic_source_is_rejected(self):
        bad = dict(self.collision[0], ttype=2)
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(writer, 'bound_excitation_atom',
                              return_value=(self.levels, self.radiation, [bad])):
                with self.assertRaisesRegex(ValueError, 'Negative nonlogarithmic'):
                    writer.write_excitation_atom('x', 2, 2, directory)
            self.assertFalse((Path(directory) / 'nlevel_x_2.txt').exists())

    def test_temperature_table_preserves_kelvin_and_strength(self):
        table = dict(self.collision[0], ttype=2,
                     temperature_grid=np.array([9800., 1.e5, 9.8e7]),
                     upsilon_grid=np.array([.1, .2, .3]))
        lower, upper, model, ttype, cups, de, interpolation, x, y = \
            writer._collision_row(self.levels, table, 1)
        self.assertEqual((lower, upper, model, ttype, cups, de, interpolation),
                         (1, 4, 'TEMPERATURE_TABLE', 0, 0., 0., 'LOG_PCHIP_T'))
        np.testing.assert_array_equal(x, table['temperature_grid'])
        np.testing.assert_array_equal(y, table['upsilon_grid'])

    def test_zero_interval_selects_pchip(self):
        transition = dict(self.collision[0], ttype=2,
                          ys=np.array([0., 0., 1.]))
        row = writer._collision_row(self.levels, transition, 1)
        self.assertEqual(row[6], 'PCHIP')
        np.testing.assert_array_equal(row[-1], transition['ys'])

    def test_degenerate_collision_preserves_source_scaling(self):
        levels = {1: (1., 0.), 4: (3., 100.), 8: (5., 100.)}
        transition = dict(self.collision[0], ll=8, ul=4)
        row = writer._collision_row(levels, transition, 1)
        self.assertEqual(row[:2], (4, 8))
        self.assertEqual(row[5], transition['de'])
        np.testing.assert_array_equal(row[-1], transition['ys'])
        # A positive radiative probability cannot describe a zero-energy line.
        with self.assertRaisesRegex(ValueError, 'Zero-energy transition'):
            writer._energy_pair(levels, 4, 8)


if __name__ == '__main__':
    unittest.main()
