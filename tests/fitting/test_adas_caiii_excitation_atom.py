"""Private Ca III ADF04 candidate: source identity, physics fields, and guards."""

from pathlib import Path
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools/fitting'))

import chianti_cooling as C
from adas_caiii_excitation_atom import (
    BOUND_LEVEL_COUNT, IONIZATION_CM1, SOURCE_SHA256,
    _COMPACT_REAL, _parse_caiii_adf04, adas_caiii_excitation_atom)


class CaIIISource(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = Path(os.environ.get('ADAS_CAIII_ADF04',
                                       '/tmp/adas_ca2_dw_ic.dat'))
        if not cls.path.is_file():
            raise unittest.SkipTest('Private OPEN-ADAS Ca III file is absent')
        cls.raw = cls.path.read_bytes()
        cls.levels, cls.radiation, cls.collisions = adas_caiii_excitation_atom(
            cls.path)

    def test_complete_bound_atom_and_source_quantities(self):
        from hashlib import sha256
        self.assertEqual(sha256(self.raw).hexdigest(), SOURCE_SHA256)
        self.assertEqual(len(self.levels), BOUND_LEVEL_COUNT)
        self.assertEqual(len(self.collisions), 741)
        self.assertEqual(len(self.radiation), 663)
        self.assertEqual(self.levels[1], (1., 0., 0.))
        self.assertEqual(self.levels[2][0], 1.)
        self.assertEqual(self.levels[3][0], 3.)
        self.assertTrue(all(self.levels[i][1] < IONIZATION_CM1
                            for i in self.levels))
        self.assertTrue(all(u > l and self.levels[u][1] > self.levels[l][1]
                            and a >= 0 for l, u, _, _, a in self.radiation))
        self.assertEqual(len({(row['ll'], row['ul']) for row in self.collisions}),
                         741)
        self.assertTrue(all(np.all(row['upsilon_grid'] > 0)
                            for row in self.collisions))
        self.assertEqual(next(a for l, u, _, _, a in self.radiation
                              if (l, u) == (1, 3)), 2.14e6)
        self.assertNotIn((1, 2), {(l, u) for l, u, _, _, _ in self.radiation})
        closed, _ = C.radiative_cascade_probabilities(
            self.levels, self.collisions, self.radiation)
        self.assertEqual(closed, [1, 2])
        marked_decay = self.radiation + [(1, 2, 1., 0., 1e-30)]
        closed_with_false_decay, _ = C.radiative_cascade_probabilities(
            self.levels, self.collisions, marked_decay)
        self.assertEqual(closed_with_false_decay, [1])

    def test_temperature_nodes_are_upsilon_not_rates(self):
        row = self.collisions[0]
        grid = row['temperature_grid']
        self.assertEqual((len(grid), grid[0], grid[-1]),
                         (14, 900., 18_000_000.))
        np.testing.assert_allclose(C.upsilon(row, grid),
                                   row['upsilon_grid'], rtol=5e-14)
        self.assertEqual(row['upsilon_grid'][0], 4.71e-2)
        self.assertEqual(row['lim'], 0.)  # Optional Bethe limit is not a node.
        lower, upper = row['ll'], row['ul']
        gap = (self.levels[upper][1]-self.levels[lower][1])*C.HC_ERG_CM
        requested = np.array([1000., 1800., 10000.])
        rate = (C.COLL_PREF*C.upsilon(row, requested)
                * np.exp(-gap/(C.K_B_ERG*requested))
                / (self.levels[lower][0]*np.sqrt(requested)))
        self.assertTrue(np.all(np.isfinite(rate)))
        self.assertTrue(np.all(rate > 0))
        self.assertTrue(np.all(np.diff(rate) > 0))
        with self.assertRaises(ValueError):
            C.upsilon(row, 18_000_001.)

    def test_source_hash_and_duplicate_pair_guards(self):
        lines = self.raw.decode('ascii').splitlines(keepends=True)
        first = next(i for i, line in enumerate(lines)
                     if line.lstrip().startswith('2    1 1.00-30'))
        lines[first+1] = lines[first]
        changed = ''.join(lines).encode('ascii')
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            _parse_caiii_adf04(changed)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'changed.dat'
            path.write_bytes(changed)
            with self.assertRaisesRegex(ValueError, 'Changed OPEN-ADAS'):
                adas_caiii_excitation_atom(path)

    def test_invalid_effective_strength_is_rejected(self):
        lines = self.raw.decode('ascii').splitlines(keepends=True)
        first = next(i for i, line in enumerate(lines)
                     if line.lstrip().startswith('2    1 1.00-30'))
        fields = list(_COMPACT_REAL.finditer(lines[first]))
        source = fields[1].group()
        lines[first] = lines[first].replace(source, '-' + source, 1)
        with self.assertRaisesRegex(ValueError, 'Nonphysical'):
            _parse_caiii_adf04(''.join(lines).encode('ascii'))


if __name__ == '__main__':
    unittest.main()
