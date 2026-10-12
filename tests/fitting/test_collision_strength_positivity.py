"""Physical nonnegativity and source-node preservation in collision rates."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/fitting'))
import chianti_cooling as C


class CollisionStrengthPositivity(unittest.TestCase):
    def test_nev_channel_remains_present_below_old_artificial_crossing(self):
        rows=C.read_scups(C.ion_dir('ne',5)+'/ne_5.scups')
        tr=next(t for t in rows if (t['ll'],t['ul'])==(6,15))
        temperature=np.geomspace(1000,100000,1601)
        self.assertTrue(np.all(C.upsilon(tr,temperature)>0))
        # Type 2 stores Upsilon itself: every original node is retained.
        interpolator=C.scaled_collision_strength_interpolator(tr)
        np.testing.assert_allclose(interpolator(tr['xs']),tr['ys'],rtol=2e-14,atol=1e-18)

    def test_logarithmic_ordinates_are_not_negative_physical_strengths(self):
        tr=dict(ll=1,ul=2,de=1.,cups=1.,ttype=6,
                xs=np.array([0.,.5,1.]),ys=np.array([-3.,-2.,-1.]))
        values=C.upsilon(tr,[1000,10000,100000])
        self.assertTrue(np.all(values>0))

    def test_invalid_negative_source_region_is_rejected(self):
        tr=dict(ll=1,ul=2,de=1.,cups=1.,ttype=2,
                xs=np.array([0.,.5,1.]),ys=np.array([1.,.5,-1.]))
        with self.assertRaisesRegex(ValueError,'negative source nodes'):
            C.upsilon(tr,[1e12])


if __name__=='__main__':unittest.main()
