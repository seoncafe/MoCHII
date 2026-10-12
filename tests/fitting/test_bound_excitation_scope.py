"""Ordinary bound excitation excludes autoionizing and unreachable states."""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/fitting'))
import chianti_cooling as C


class BoundExcitationScope(unittest.TestCase):
    def test_autoionizing_and_unpopulated_states_are_excluded(self):
        # C II has a physical edge near 24 eV. State 3 is well above it;
        # state 4 is bound but has no upward route from ground-connected states.
        raw={1:(1.,0.,0.),2:(3.,10000.,10000.),3:(5.,500000.,500000.),4:(3.,15000.,15000.)}
        rad=[(1,2,1.,0.,1.),(1,3,1.,0.,1.),(2,4,1.,0.,1.)]
        coll=[dict(ll=1,ul=u,de=(raw[u][1])/109737.,ttype=2,cups=1.,xs=np.array([0.,1.]),ys=np.ones(2)) for u in [2,3]]
        with patch.object(C,'read_excitation_atom',return_value=(raw,rad,coll)):
            levels,radiation,transitions=C.bound_excitation_atom('c',2)
            measured=C.cooling_low_density_limit('c',2,[10000.])[0]
        self.assertEqual(set(levels),{1,2});self.assertEqual(len(radiation),1);self.assertEqual(len(transitions),1)
        gap=10000*C.HC_ERG_CM
        expected=C.COLL_PREF*gap/100*np.exp(-gap/(C.K_B_ERG*10000))
        self.assertAlmostEqual(measured/expected,1.,places=14)


if __name__=='__main__':unittest.main()
