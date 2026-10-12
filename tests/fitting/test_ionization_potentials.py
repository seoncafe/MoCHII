"""Physical valence thresholds and retained independent fit energies."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/fitting'))
import ionization_potentials as I
import make_element_data as G


class IonizationPotentials(unittest.TestCase):
    def test_nist_identity_and_unit_conversion(self):
        values = I.read_ionization_potentials()
        self.assertAlmostEqual(values[6,6],90820.348*I.HC_EV_CM,places=12)
        self.assertAlmostEqual(values[18,15],480560.0*I.HC_EV_CM,places=10)
        for element,z in I.ELEMENTS:
            for stage in range(1,min(z,7)+1):
                self.assertGreater(values[z,z-stage+1],0)

    def test_occupied_valence_edges_and_independent_ci(self):
        ph,ninn,ntot = G.parse_ph1_phfit2()
        ci = G.parse_cf()
        for el,z,n in G.ELEMENTS:
            for stage in range(1,n):
                ne = z-stage+1
                inner,rows = G.phfit2_shells(ph,ninn,ntot,z,ne)
                self.assertEqual(rows[-1][1],G.ETH[z,ne])
                for row in rows[:-1]:
                    # Corresponding coefficients preserve the source inner edge.
                    self.assertTrue(any(row[1:]==v for key,v in ph.items()
                                        if key[:2]==(z,ne)))
        self.assertEqual(ci[18,15][0],52.3)
        self.assertGreater(G.ETH[18,15],59)

    def test_charge_exchange_term_reader_uses_fixed_columns(self):
        self.assertEqual(G.term_J('n_4',[3]),[1.0])
        self.assertAlmostEqual(G.chianti_term_energy('n_4',[3],'2s 2p 3P'),
                               67272.*G.HC_EV_CM,places=12)

    def test_source_version_duplicate_and_footer(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'potentials.dat'
            path.write_text('20221115\n 6 0 90820.348 0.009\n***\nnot numeric bibliography\n')
            self.assertEqual(len(I.read_ionization_potentials(path)),1)
            path.write_text('20221115\n6 0 1 0\n6 0 2 0\n')
            with self.assertRaisesRegex(ValueError,'duplicate'):
                I.read_ionization_potentials(path)
            path.write_text('20220101\n6 0 1 0\n')
            with self.assertRaisesRegex(ValueError,'snapshot'):
                I.read_ionization_potentials(path)


if __name__=='__main__':
    unittest.main()
