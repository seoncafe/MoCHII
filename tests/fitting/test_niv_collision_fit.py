"""N IV low-level strength peaks must fit inside native coefficient capacity.

Historical (2026-10-11): this checked the Chebyshev N IV atom of
fit_nlevel.py, which no longer produces any installed N IV atom. Both
excitation models install source-node atoms (write_excitation_atom.py,
data/atomic/excitation/<model>/nlevel_n_4.txt) that store the original
collision-strength nodes, so no coefficient-capacity fit is made. Measured
before it was skipped: 0.39% against the 0.2% fit_nlevel.ERR_TARGET.
fit_nlevel.py remains the generator of the H I atom only.
"""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/fitting'))
import fit_nlevel as N


@unittest.skip('historical: the fit_nlevel Chebyshev N IV atom is not installed; see the module docstring')
class NivCollisionFit(unittest.TestCase):
    def test_metastable_strength_peaks_on_dense_temperature_grid(self):
        source=Path(N.ion_dir('n',4))/'n_4.scups'
        selected=[tr for tr in N.read_scups(source) if tr['ll']==1 and tr['ul'] in (2,3,4)]
        self.assertEqual(len(selected),3)
        temperatures=np.geomspace(1e3,1e5,1601)
        for tr in selected:
            lo,hi,coefficients,logarithmic,error=N.fit_transition(tr)
            self.assertLessEqual(len(coefficients),16)
            fitted=N.eval_ups_fit(tr,lo,hi,coefficients,logarithmic,temperatures)
            reference=N.upsilon(tr,temperatures)
            self.assertTrue(np.all(fitted>0))
            self.assertLess(np.max(abs(fitted/reference-1)),N.ERR_TARGET)


if __name__=='__main__':
    unittest.main()
