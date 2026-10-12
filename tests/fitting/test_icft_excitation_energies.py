"""Observed/ICFT identities and independent stationary cooling checks."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/fitting'))
import chianti_cooling as C
from atomic_excitation_energies import selected_energies
from test_cooling_stationary_limit import finite_density_cooling


class ICFTExcitation(unittest.TestCase):
    def atom(self,el,stage):
        p=Path(C.ion_dir(el,stage))/f'{el}_{stage}'
        raw=C.read_elvlc(str(p)+'.elvlc');rad=C.read_wgfa(str(p)+'.wgfa');coll=C.read_scups(str(p)+'.scups')
        return raw,selected_energies(el,stage,raw,rad,coll),rad,coll

    def test_ov_bound_mapping_restores_rounded_gaps(self):
        raw,levels,rad,coll=self.atom('o',5)
        self.assertEqual(len(levels),166)
        for i in range(1,41):self.assertEqual(levels[i][1],raw[i][1])
        self.assertGreater(levels[73][1]-levels[72][1],0)
        self.assertGreater(levels[143][1]-levels[142][1],0)
        self.assertAlmostEqual(levels[87][1],832807.7824)

    def test_clvi_preserves_low_observations_and_high_icft(self):
        raw,levels,rad,coll=self.atom('cl',6)
        for i in range(1,45):self.assertEqual(levels[i][1],raw[i][1])
        self.assertEqual(levels[79][1],raw[79][2])
        self.assertGreater(levels[80][1],levels[79][1])

    def test_full_finite_density_balance_matches_fit_target(self):
        # The cooling fits are made to the stationary state at n_e = 1e-6
        # cm^-3 (chianti_cooling.reference_log_cooling; decision of
        # 2026-10-11), not to the strict n_e -> 0 limit, which a placeholder
        # spontaneous rate can drain. An independently assembled full balance
        # at that density must reproduce the target.
        for el,stage in [('o',5),('cl',6)]:
            lev,rad,coll=C.bound_excitation_atom(el,stage)
            temperatures=np.array([2000.,10000.,100000.])
            target=np.exp(C.reference_log_cooling(lev,rad,coll,temperatures)[0])
            actual=np.array([finite_density_cooling(lev,rad,coll,t,C.REFERENCE_DENSITY)
                             for t in temperatures])
            np.testing.assert_allclose(actual,target,rtol=1e-9,atol=0)


if __name__=='__main__':unittest.main()
