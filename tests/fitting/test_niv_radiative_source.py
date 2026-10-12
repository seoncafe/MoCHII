"""N IV energies must preserve radiative source indices and output integrity."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools/fitting'))
import fit_nlevel as F
import niv_level_energies as N
from atomic_excitation_energies import selected_energies
import chianti_cooling as C


class NivRadiativeSource(unittest.TestCase):
    def test_radiative_indices_and_energies(self):
        stem = Path(C.ion_dir('n', 4))/'n_4'
        raw = C.read_elvlc(str(stem)+'.elvlc')
        radiation = C.read_wgfa(str(stem)+'.wgfa')
        collisions = C.read_scups(str(stem)+'.scups')
        levels = selected_energies('n', 4, raw, radiation, collisions)
        self.assertEqual(len(levels), 238)
        self.assertAlmostEqual(levels[107][1], 588134.6)
        self.assertTrue(all(levels[u][1] > levels[l][1] for l,u,*_ in radiation))
        self.assertTrue(all(levels[tr['ul']][1] > levels[tr['ll']][1]
                            for tr in collisions))
        with tempfile.TemporaryDirectory() as directory:
            # Fitting spline coefficients is independent of the energy choice.
            with patch.object(F, 'fit_transition', return_value=(0., 1., [1.], 0, 0.)):
                path, nl, nr, nu, _ = F.build_ion('NIV', 'n', 4, 136,
                                                   output_dir=directory)
            self.assertEqual((nl, nr, nu), (136, 3272, 9180))
            written = Path(path).read_text()
            self.assertIn('107 116', written)
            self.assertIn('588134.6000', written)
            self.assertIn(N.STOUT_SHA256, written)

    def test_changed_source_preserves_existing_output(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'nlevel_n_4.txt'
            path.write_text('preserved\n')
            with patch.object(N, 'STOUT_SHA256', '0'*64):
                with self.assertRaisesRegex(ValueError, 'Stout energy source hash changed'):
                    F.build_ion('NIV', 'n', 4, 136, output_dir=directory)
            self.assertEqual(path.read_text(), 'preserved\n')


if __name__ == '__main__':
    unittest.main()
