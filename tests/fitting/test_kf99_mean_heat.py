"""Validate the separately published mean heat and selective atomic generation."""
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/fitting'))
import charge_exchange_expansion as R
import make_element_data as G
from kf96_rate_table import kf96_charge_transfer_table, charge_transfer_rate_on_grid


class MeanChargeTransferHeat(unittest.TestCase):
    def test_published_heat_and_energy_budget(self):
        # Independently transcribed from the visually checked KF99 Table 1.
        expected={('c',4):11.30,('n',3):6.40,('o',3):5.00,('o',4):8.47,
                  ('ne',3):5.82,('ne',4):8.60,('s',4):8.60,
                  ('ar',4):8.60,('mg',3):5.73,('mg',4):8.60}
        self.assertEqual(R.KF99_ADDED_HEAT,expected)
        for key,heat in expected.items():
            el,q=key; z=R.ATOMIC_NUMBER[el]
            dip=G.ETH[z,z-q+1]-G.IP_H
            ion,rec,split,note,table=G.charge_exchange_heat_row(el,q,G.ETH[z,z-q+1],False,True)
            self.assertAlmostEqual(rec,heat,places=12)
            self.assertEqual((ion,split,table),(0,0,None))
            self.assertGreater(dip-rec,0)  # Mean excitation, not an assigned level.
            self.assertIn('mean heat',note)

    def test_added_tables_over_full_published_intervals(self):
        for key,reaction in R.kf96_added_rates().items():
            ts,ks=kf96_charge_transfer_table(reaction)
            ts=[float(f'{v:.6e}') for v in ts];ks=[float(f'{v:.6e}') for v in ks]
            self.assertLessEqual(len(ts),1024)
            for j in range(2001):
                t=reaction.min_temperature*(reaction.max_temperature/reaction.min_temperature)**(j/2000)
                self.assertLess(abs(charge_transfer_rate_on_grid(t,ts,ks)/reaction.rate(t)-1),1e-5)

    def test_all_ten_selected_and_missing_heat_still_rejected(self):
        roster=[('c',6,5),('n',7,4),('o',8,5),('ne',10,5),
                ('s',16,5),('ar',18,5),('mg',12,5)]
        with tempfile.TemporaryDirectory() as directory:
            G.generate(roster,directory)
            for (el,q),heat in R.KF99_ADDED_HEAT.items():
                lines=(Path(directory)/f'element_{el}.txt').read_text().splitlines()
                start=lines.index(f'TRANSITION {q}')
                row=next(v for v in lines[start:] if v.startswith('CXHEAT ')).split()
                self.assertEqual(float(row[2]),heat)
                self.assertTrue(any(v.startswith('CXRTAB ') for v in lines[start:]))
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'element_c.txt';path.write_text('preserved\n')
            absent=dict(G.CX_HEAT);absent.pop(('c',4))
            with patch.dict(G.CX_HEAT,absent,clear=True):
                with self.assertRaisesRegex(ValueError,'c transition 4: carried rec charge exchange has no CX_HEAT entry'):
                    G.generate([('c',6,5)],directory)
            self.assertEqual(path.read_text(),'preserved\n')


if __name__=='__main__': unittest.main()
