"""Physical excitation tails remain finite in log space without a rate floor."""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/fitting'))
import chianti_cooling as C
import fit_cooling as F


class ExcitationUnderflow(unittest.TestCase):
    def test_ground_excitation_matches_exact_logarithmic_coefficient(self):
        energy=2e6/C.HC_OVER_K
        levels={1:(1.,0.,0.),2:(3.,energy,energy)}
        radiation=[(1,2,1.,0.,1.)]
        collision=dict(ll=1,ul=2,de=energy/109737.,ttype=2,cups=1.,
                       xs=np.array([0.,1.]),ys=np.array([1.,1.]))
        with patch.object(C,'read_elvlc',return_value=levels), \
             patch.object(C,'read_wgfa',return_value=radiation), \
             patch.object(C,'read_scups',return_value=[collision]):
            source=C.cooling_low_density_limit('synthetic',1,F.T)
            logarithmic=C.log_cooling_ground_state_limit('synthetic',1,F.T)
        amplitude=C.COLL_PREF*energy*C.HC_ERG_CM
        exact=np.log(amplitude)-.5*np.log(F.T)-2e6/F.T
        np.testing.assert_allclose(logarithmic,exact,rtol=1e-14)
        self.assertEqual(source[0],0)
        self.assertTrue(np.all(np.isfinite(logarithmic)))
        a,t=F.fit_ground_state_log_cooling(logarithmic,1,[2.1e6])
        self.assertLess(F.log_relative_cooling_error(a,t,logarithmic),1e-10)

    def test_real_high_ions_have_logarithmic_excitation_tails(self):
        for element,stage in [('c',5),('n',6),('s',7)]:
            temperatures=np.array([1000.,10000.,100000.])
            source=C.cooling_low_density_limit(element,stage,temperatures)
            logs=C.log_cooling_ground_state_limit(element,stage,temperatures)
            self.assertEqual(source[0],0)
            self.assertTrue(np.all(np.isfinite(logs)))
            np.testing.assert_allclose(np.exp(logs[1:]),source[1:],rtol=1e-12)


if __name__=='__main__':unittest.main()
