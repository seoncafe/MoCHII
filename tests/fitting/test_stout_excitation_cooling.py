"""Cl VII bound-state scope and the published collision-temperature interval."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/fitting'))
import chianti_cooling as C
import make_element_data as G
from stout_excitation_cooling import clvii_outer_atom,clvii_outer_cooling,effective_collision_strength_temperature_table


class StoutExcitation(unittest.TestCase):
    def test_bound_model_excludes_core_states_above_threshold(self):
        levels,radiation,grid,collisions=clvii_outer_atom()
        threshold=G.ETH[17,11]
        self.assertEqual(len(levels),32)
        self.assertTrue(all(e*C.HC_ERG_CM/C.EV_ERG<threshold for g,e in levels.values()))
        self.assertTrue(all(l<=32 and u<=32 for l,u,a in radiation))

    def test_invalid_temperatures_are_rejected(self):
        for t in [0.,-1.,98000001.,float('nan')]:
            with self.assertRaises(ValueError):clvii_outer_cooling([t])
        value=clvii_outer_cooling([9800.,24500.,49000.,98000.])
        self.assertTrue(np.all(np.isfinite(value)))
        self.assertTrue(np.all(value>0))

    def test_cold_strength_is_constant_but_rates_use_actual_temperature(self):
        levels,rad,grid,coll=clvii_outer_atom()
        temperature=np.array([1000.,2000.,9799.,9800.])
        expected=np.zeros_like(temperature)
        for lower,upper,strength in coll:
            actual=effective_collision_strength_temperature_table(grid,strength,temperature)
            np.testing.assert_allclose(actual,strength[0],rtol=2e-15)
            if lower==1:
                gap=(levels[upper][1]-levels[1][1])*C.HC_ERG_CM
                expected+=C.COLL_PREF*strength[0]*gap/(levels[1][0]*np.sqrt(temperature))*np.exp(-gap/(C.K_B_ERG*temperature))
        np.testing.assert_allclose(clvii_outer_cooling(temperature),expected,rtol=3e-15)
        np.testing.assert_allclose(C.cooling_low_density_limit('cl',7,temperature),expected,rtol=4e-15)
        self.assertLess(expected[0],expected[-1]*1e-20)

    def test_published_temperature_nodes_are_preserved(self):
        levels,rad,grid,coll=clvii_outer_atom()
        for lower,upper,strength in coll:
            np.testing.assert_allclose(effective_collision_strength_temperature_table(grid,strength,grid),strength,rtol=3e-15)

    def test_generic_source_ceiling_accepts_only_rounding(self):
        grid=np.array([1000.,10000.,100000.])
        strength=np.array([1.,2.,3.])
        rounded=np.nextafter(grid[-1],np.inf)
        result=effective_collision_strength_temperature_table(grid,strength,[grid[-1],rounded])
        np.testing.assert_array_equal(result,[result[0],result[0]])
        with self.assertRaises(ValueError):
            effective_collision_strength_temperature_table(grid,strength,[100000.1])


if __name__=='__main__':unittest.main()
