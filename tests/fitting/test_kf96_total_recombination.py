"""Further published total rates, uncertain channels, and existing omissions."""
import math
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/fitting'))
import charge_exchange_expansion as R
import make_element_data as G
from kf96_rate_table import kf96_charge_transfer_table, charge_transfer_rate_on_grid

class PublishedTotalRecombination(unittest.TestCase):
    def test_independent_published_totals_and_means(self):
        # Direct KF96 Tables 1/2 and KF99 Table 1 transcriptions, except N4+ +
        # H: the updated KF96 database (Ferland et al. 2013, Cloudy
        # ctrecombdata.dat, 1e3-1e6 K) in place of the printed KF96 Table 1
        # fit (3.28, 0.52, -0.52, -0.19; 1e3-3e5 K) (author's decision,
        # 2026-10-11).
        expected={('n',4):(2.95,.55,-.39,-1.07),('cl',3):(1.88,.32,1.77,-5.70),
                  ('cl',4):(7.27,.29,1.04,-10.14),('ca',3):(.0317,2.12,12.06,-.40),
                  ('ca',4):(2.68,.69,-.68,-4.47),('fe',4):(14.60,.0357,-.92,-.37)}
        heats={('n',4):11.,('cl',3):8.,('cl',4):9.,('ca',3):6.6,('ca',4):9.9,('fe',4):10.}
        rows=R.kf96_recombination_total_rates()
        self.assertEqual(R.KF99_RECOMBINATION_MEAN_HEAT,heats)
        for key,fit in expected.items():
            self.assertEqual(rows[key].coefficients,fit)
            self.assertEqual(rows[key].min_temperature,1000.)
            self.assertEqual(rows[key].max_temperature,1e6 if key==('n',4) else 3e4)
            t4=1.;a,b,c,d=fit
            self.assertAlmostEqual(rows[key].rate(10000)/(a*1e-9*(1+c*math.exp(d))),1.,places=13)
            z=R.ATOMIC_NUMBER[key[0]];threshold=G.ETH[z,z-key[1]+1]
            _,heat,split,_,_=G.charge_exchange_heat_row(*key,threshold,False,True)
            self.assertAlmostEqual(heat,heats[key],places=12)
            self.assertEqual(split,0.)
        self.assertEqual(rows['n',4].coefficients,R.cloudy_recombination_rows()['n',4].coefficients)

    def test_full_domain_serialized_tables(self):
        for key,row in R.kf96_recombination_total_rates().items():
            ts,ks=kf96_charge_transfer_table(row)
            ts=[float(f'{x:.6e}') for x in ts];ks=[float(f'{x:.6e}') for x in ks]
            self.assertLessEqual(len(ts),1024)
            for i in range(2001):
                t=row.min_temperature*(row.max_temperature/row.min_temperature)**(i/2000)
                self.assertLess(abs(charge_transfer_rate_on_grid(t,ts,ks)/row.rate(t)-1),1e-5)

    def test_existing_omissions_and_selected_nitrogen(self):
        with tempfile.TemporaryDirectory() as directory:
            G.generate([('cl',17,5),('ca',20,5),('n',7,5)],directory)
            for el,q in [('cl',3),('cl',4),('ca',3),('ca',4),('n',4)]:
                lines=(Path(directory)/f'element_{el}.txt').read_text().splitlines()
                start=lines.index(f'TRANSITION {q}')
                stop=next((i for i in range(start+1,len(lines)) if lines[i].startswith('TRANSITION ')),len(lines))
                block=lines[start:stop]
                self.assertIn('CXR 8 0 0 0 0 0 0',block)
                self.assertTrue(any(v.startswith('CXRTAB ') for v in block))
                heat=next(v for v in block if v.startswith('CXHEAT ')).split()
                self.assertEqual(float(heat[2]),R.KF99_RECOMBINATION_MEAN_HEAT[el,q])
                if el=='cl': self.assertTrue(any('selected inside sampled' in v for v in block))
        with tempfile.TemporaryDirectory() as directory:
            G.generate([('fe',26,5)],directory)
            self.assertIn('Zhang et al. 2026',(Path(directory)/'element_fe.txt').read_text())
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError,'unresolved added-boundary'):
                G.generate([('fe',26,8)],directory)
            self.assertFalse((Path(directory)/'element_fe.txt').exists())

if __name__=='__main__': unittest.main()
