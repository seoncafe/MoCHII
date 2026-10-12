"""Published Fe limit and reduced-mass capture physics, without plasma activation."""
from pathlib import Path
import sys
import unittest
import math
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/fitting'))
import charge_exchange_capture as C
import make_element_data as G
import tempfile

class PolarizationCapture(unittest.TestCase):
    def test_published_rounded_iron_limit(self):
        for charge in (5,6):
            self.assertLess(abs(C.langevin_hydrogen_capture_bound(charge,55.845)/(1.92e-9*charge)-1),.002)

    def test_reduced_mass_and_charge_dependence(self):
        iron=C.langevin_hydrogen_capture_bound(5,55.845)
        nitrogen=C.langevin_hydrogen_capture_bound(5,14.0067)
        self.assertAlmostEqual(nitrogen/iron,math.sqrt(C.hydrogen_ion_reduced_mass(55.845)/C.hydrogen_ion_reduced_mass(14.0067)),places=13)
        self.assertGreater(nitrogen,iron)
        self.assertAlmostEqual(C.langevin_hydrogen_capture_bound(6,55.845)/iron,6/5,places=13)

    def test_characteristic_barrier_scaling(self):
        a=C.polarization_capture_conditions(5,55.845,1000.)
        b=C.polarization_capture_conditions(5,55.845,16000.)
        self.assertAlmostEqual(a['barrier_radius_a0']/b['barrier_radius_a0'],2.)
        self.assertAlmostEqual(b['critical_angular_momentum_hbar']/a['critical_angular_momentum_hbar'],2.)

    def test_invalid_collision_conditions(self):
        for q in (0,-1,True,5.5):
            with self.assertRaises(ValueError):C.langevin_hydrogen_capture_bound(q,55.845)
        for mass in (0,-1,float('nan'),float('inf')):
            with self.assertRaises(ValueError):C.langevin_hydrogen_capture_bound(5,mass)
        for t in (0,-1,float('nan'),float('inf')):
            with self.assertRaises(ValueError):C.polarization_capture_conditions(5,55.845,t)

    def test_selected_rates_and_separate_mean_heat(self):
        rates=C.langevin_recombination_rates()
        self.assertEqual(len(rates),22)
        for (el,q),rate in rates.items():
            self.assertEqual(G.CX_HIGH[el,q],(1,[rate,0,0,0,0,0]))
            atomic_number=next(n for name,n,_ in G.ELEMENTS if name==el)
            edge=G.ETH[atomic_number,atomic_number-q+1]
            ion,heat,split,note,table=G.charge_exchange_heat_row(el,q,edge,False,True)
            self.assertEqual((ion,split),(0,0))
            self.assertAlmostEqual(heat,2.86*(q-1),places=12)
            exit_coulomb_eV=(q-1)*C.ELEMENTARY_CHARGE_C/(4*math.pi*C.EPSILON_0*9.5*C.BOHR_M)
            self.assertLess(abs(heat/exit_coulomb_eV-1),0.002)
            self.assertLess(heat,edge-G.IP_H)
            self.assertIn('upper limit',note)
        for q in (4,True,5.5):
            with self.assertRaises(ValueError):C.avoided_crossing_mean_heat(q)

    def test_generated_high_charge_reactions(self):
        with tempfile.TemporaryDirectory() as directory:
            G.generate([('n',7,6),('fe',26,7)],directory)
            for el,q in [('n',5),('fe',5),('fe',6)]:
                text=(Path(directory)/f'element_{el}.txt').read_text()
                block=text.split(f'TRANSITION {q}\n')[1].split('TRANSITION ')[0]
                self.assertIn('CXR 1 ',block)
                self.assertIn('unit-transfer probability',block)
                self.assertIn(f'CXHEAT 0.000000 {2.86*(q-1):.6f} 0.000000',block)

if __name__=='__main__':unittest.main()
