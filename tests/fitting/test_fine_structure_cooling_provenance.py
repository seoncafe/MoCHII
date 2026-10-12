"""Protect atomic source files while writing staged fine-structure fits."""

from hashlib import sha256
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools' / 'fitting'))
import fit_cooling
from write_excitation_atom import _source_files


class FineStructureCoolingProvenanceTest(unittest.TestCase):
    def test_staging_output_does_not_replace_a_source(self):
        for name, element, stage in (('ArII', 'ar', 2), ('CaIV', 'ca', 4)):
            with self.subTest(ion=name), tempfile.TemporaryDirectory() as directory:
                source_paths = _source_files(element, stage)
                before = {path: sha256(path.read_bytes()).hexdigest()
                          for path in source_paths}
                expected = Path(directory) / f'cooling_{element}_{stage}.txt'
                actual = fit_cooling.write_ion_file(
                    name, element, stage, np.array([1.e-20]),
                    np.array([1000.]), 0., 0., output_dir=directory)
                self.assertEqual(Path(actual), expected)
                self.assertTrue(expected.is_file())
                self.assertIn('Pelan and Berrington', expected.read_text())
                after = {path: sha256(path.read_bytes()).hexdigest()
                         for path in source_paths}
                self.assertEqual(after, before)


if __name__ == '__main__':
    unittest.main()
