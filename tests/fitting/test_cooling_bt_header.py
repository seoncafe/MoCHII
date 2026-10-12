"""Cooling coefficient headers retain identified BT temperature limits."""

from hashlib import sha256
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools/fitting'))
import fit_cooling as F
from chianti_cooling import ion_dir


class CoolingBtHeaderTest(unittest.TestCase):
    def test_identified_cold_source_metadata_and_source_preservation(self):
        for name, element, stage, minimum, maximum in (
                ('CIV', 'c', 4, 3200, 32000000),
                ('MgIII', 'mg', 3, 1800, 18000000),
                ('SiV', 'si', 5, 5000, 50000000)):
            with self.subTest(ion=name), tempfile.TemporaryDirectory() as directory:
                raw = Path(ion_dir(element, stage)) / f'{element}_{stage}.scups'
                before = sha256(raw.read_bytes()).hexdigest()
                amplitude = np.array([1e-19, -1e-21])
                excitation = np.array([10000., 20000.])
                product = Path(F.write_ion_file(name, element, stage,
                                                amplitude, excitation,
                                                .002, .001, directory))
                content = product.read_text()
                self.assertEqual(product.parent, Path(directory))
                self.assertIn(f'published BT collision temperature interval: '
                              f'{minimum}--{maximum} K', content)
                self.assertIn(f'constant source Upsilon below {minimum} K; '
                              'actual T enters Maxwellian rates', content)
                self.assertIn('cold collision-strength accuracy is unverified',
                              content)
                numeric = [line.split() for line in content.splitlines()
                           if line.strip() and not line.startswith('#')]
                self.assertEqual(int(numeric[0][0]), 2)
                self.assertEqual(float(numeric[2][0]), amplitude[1])
                self.assertEqual(sha256(raw.read_bytes()).hexdigest(), before)

    def test_no_floor_source_does_not_claim_cold_continuation(self):
        with tempfile.TemporaryDirectory() as directory:
            product = Path(F.write_ion_file('CIII', 'c', 3,
                                            np.array([1e-20]),
                                            np.array([10000.]),
                                            .001, .001, directory))
            content = product.read_text()
            self.assertNotIn('published BT collision temperature interval', content)
            self.assertNotIn('adopted cold approximation', content)
            numeric = [line.split() for line in content.splitlines()
                       if line.strip() and not line.startswith('#')]
            self.assertEqual(int(numeric[0][0]), 1)
            self.assertEqual(float(numeric[1][0]), 1e-20)


if __name__ == '__main__':
    unittest.main()
