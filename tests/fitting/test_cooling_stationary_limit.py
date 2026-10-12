"""Steady low-density cooling against finite-density statistical equilibrium."""
from pathlib import Path
import sys
import unittest
import tempfile
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools/fitting'))
import chianti_cooling as C
from atomic_excitation_energies import selected_energies


def finite_density_cooling(levels, radiation, collisions, temperature, density,
                           return_population=False):
    # Independently assemble the full level balance, including all excited
    # populations and their spontaneous photon emission, at finite n_e.
    active = sorted({i for tr in collisions for i in (tr['ll'],tr['ul'])} |
                    {i for l,u,w,g,a in radiation for i in (l,u)})
    pos = {i:k for k,i in enumerate(active)}
    matrix = np.zeros((len(active),len(active)))
    power = np.zeros(len(active))

    def rate(source, destination, value):
        matrix[pos[destination],pos[source]] += value
        matrix[pos[source],pos[source]] -= value

    for l,u,w,g,a in radiation:
        if a > 0:
            rate(u,l,a)
            power[pos[u]] += a*(levels[u][1]-levels[l][1])*C.HC_ERG_CM
    for tr in collisions:
        l,u = sorted((tr['ll'],tr['ul']),key=lambda k:(levels[k][1],k))
        delta = C.transition_energy_erg(levels,tr)
        common = C.COLL_PREF*C.upsilon(tr,temperature)[0]/np.sqrt(temperature)
        rate(l,u,density*common/levels[l][0]*np.exp(-delta/(C.K_B_ERG*temperature)))
        rate(u,l,density*common/levels[u][0])
    scale = np.max(np.abs(matrix),axis=1)
    matrix /= np.where(scale > 0,scale,1)[:,None]
    matrix[0,:] = 1
    rhs = np.zeros(len(active)); rhs[0] = 1
    fractions = np.linalg.solve(matrix,rhs)
    cooling = fractions@power/density
    if return_population:
        return cooling, fractions
    return cooling


class StationaryCooling(unittest.TestCase):
    def test_radiative_direction_is_independent_of_index_order(self):
        levels={1:(1.,0.),2:(5.,10000.),3:(3.,1000.)}
        radiation=[(1,2,1.,0.,1.),(3,2,1.,0.,2.)]
        transitions=[{'ll':1,'ul':2},{'ll':1,'ul':3}]
        closed,probability=C.radiative_cascade_probabilities(levels,transitions,radiation)
        self.assertEqual(closed,[1,3])
        np.testing.assert_allclose(probability[2],[1/3,2/3])

    def test_collision_pair_direction_is_independent_of_index_order(self):
        levels={1:(1.,0.,0.),2:(5.,10000.,10000.),3:(3.,1000.,1000.)}
        energy=C._level_energy_cm1(levels)
        radiation=[(1,2,1.,0.,1.),(3,2,1.,0.,2.)]
        pairs=[(1,2),(1,3),(2,3)]
        coll=[dict(ll=l,ul=u,de=abs(energy[u][1]-energy[l][1])/109737.,
            ttype=2,cups=1.,xs=np.array([0.,1.]),ys=np.ones(2)) for l,u in pairs]
        temperature=np.array([1e3,1e4,1e5])
        with patch.object(C,'read_elvlc',return_value=levels),patch.object(C,'read_wgfa',return_value=radiation),patch.object(C,'read_scups',return_value=coll):
            measured=C.cooling_low_density_limit('synthetic',1,temperature)
        expected=np.array([finite_density_cooling(energy,radiation,coll,t,1e-6) for t in temperature])
        np.testing.assert_allclose(measured,expected,rtol=1e-9)
        self.assertGreater(C.transition_energy_erg(energy,coll[-1]),0)

    def test_two_closed_states(self):
        levels = {1:(1.,0.,0.),2:(3.,1000.,1000.),3:(5.,10000.,10000.)}
        energy = C._level_energy_cm1(levels)
        radiation = [(1,3,1.,0.,1.),(2,3,1.,0.,2.)]
        transitions = []
        for l,u in [(1,2),(1,3),(2,3)]:
            transitions.append(dict(ll=l,ul=u,de=(energy[u][1]-energy[l][1])/109737.,
                                    ttype=2,cups=1.,xs=np.array([0.,1.]),ys=np.array([1.,1.])))
        temperatures = np.array([1e3,1e4,1e5])
        with patch.object(C,'read_elvlc',return_value=levels), \
             patch.object(C,'read_wgfa',return_value=radiation), \
             patch.object(C,'read_scups',return_value=transitions):
            limit = C.cooling_low_density_limit('synthetic',1,temperatures)
            imposed_ground = C.cooling_effective('synthetic',1,temperatures,pop='coronal')
        reference = np.array([finite_density_cooling(energy,radiation,transitions,t,1e-6)
                              for t in temperatures])
        np.testing.assert_allclose(limit,reference,rtol=1e-9)
        self.assertGreater(np.max(abs(limit/imposed_ground-1)),.01)

    def test_current_ciii_and_siiii(self):
        for element,stage in [('c',3),('si',3)]:
            levels,radiation,transitions = C.bound_excitation_atom(element,stage)
            temperatures = np.array([1e3,1e4,1e5])
            limit = C.cooling_low_density_limit(element,stage,temperatures)
            reference = np.array([finite_density_cooling(levels,radiation,transitions,t,1e-12)
                                  for t in temperatures])
            np.testing.assert_allclose(limit,reference,rtol=1e-8)
            self.assertTrue(np.all(limit > 0))

    def test_elvlc_optional_columns(self):
        main = f"{3:7d}{'2s 2p':30s}{'':5s}{3:5d}{'P':5s}{1.:5.1f}{67272.:15.3f}{67458.:15.3f}"
        for extra in ['', ', 123.0, 456.0', '       3    2    1']:
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory)/'atom.elvlc'
                path.write_text(main+extra+'\n -1\n')
                self.assertEqual(C.read_elvlc(path), {3:(3.,67272.,67458.)})

    def test_niv_energy_columns_and_consistent_source(self):
        stem = Path(C.ion_dir('n',4))/'n_4'
        raw = C.read_elvlc(str(stem)+'.elvlc')
        self.assertEqual(raw[3], (3.,67272.,67458.))
        radiation = C.read_wgfa(str(stem)+'.wgfa')
        transitions = C.read_scups(str(stem)+'.scups')
        levels = selected_energies('n',4,raw,radiation,transitions)
        temperatures = np.array([1e3,1e4,1e5])
        limit = C.cooling_low_density_limit('n',4,temperatures)
        self.assertTrue(np.all(np.isfinite(limit)))
        self.assertTrue(np.all(limit > 0))
        finite = [finite_density_cooling(levels,radiation,transitions,t,1e-12,
                                         return_population=True)
                  for t in temperatures]
        reference = np.array([cooling for cooling, _ in finite])
        np.testing.assert_allclose(limit,reference,rtol=1e-8)
        for _, populations in finite:
            self.assertTrue(np.all(np.isfinite(populations)))
            self.assertGreaterEqual(np.min(populations), -64*np.finfo(float).eps)
            self.assertAlmostEqual(np.sum(populations), 1.0)
        effective = C.cooling_effective('n',4,temperatures,pop='coronal')
        self.assertTrue(np.all(np.isfinite(effective)))
        self.assertTrue(np.all(effective > 0))

    def test_single_closed_state_reduces_to_ground_excitation(self):
        temperature = np.logspace(3,5,51)
        expected = C.cooling_effective('c',2,temperature,pop='coronal')
        actual = C.cooling_low_density_limit('c',2,temperature)
        np.testing.assert_allclose(actual,expected,rtol=1e-14)


if __name__ == '__main__':
    unittest.main()
