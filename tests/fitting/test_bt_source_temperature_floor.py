"""Physical source-temperature continuation for identified BT collisions."""

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools/fitting"))
import chianti_cooling as C
import write_excitation_atom as writer


class BtSourceTemperatureFloorTest(unittest.TestCase):
    def test_identified_source_records_only(self):
        for element, stage, minimum in (("c", 4, 3200.),
                                         ("mg", 3, 1800.),
                                         ("si", 5, 5000.)):
            with self.subTest(element=element, stage=stage):
                _, _, transitions = C.read_excitation_atom(element, stage)
                self.assertGreater(len(transitions), 0)
                self.assertTrue(all(tr.get("bt_tmin_K") == minimum
                                    for tr in transitions))
        _, _, other = C.read_excitation_atom("c", 3)
        self.assertTrue(all("bt_tmin_K" not in tr for tr in other))

    def test_six_bt_types_keep_cold_upsilon_constant(self):
        temperature = np.array([500., 1000., 3199.9, 3200., 3200.1, 6000.])
        x = np.array([0., .2, .5, .8, 1.])
        for kind in range(1, 7):
            with self.subTest(kind=kind):
                y = np.array([-1., -.8, -1.1, -.5, -.4]) if kind == 6 else \
                    np.array([.1, .3, .2, .5, .4])
                record = dict(ll=1, ul=2, de=.02, cups=2., ttype=kind,
                              xs=x, ys=y, bt_tmin_K=3200.)
                measured = C.upsilon(record, temperature)
                self.assertTrue(np.all(np.isfinite(measured)))
                self.assertTrue(np.all(measured >= 0))
                np.testing.assert_array_equal(measured[:4],
                                              np.repeat(measured[3], 4))
                original = dict(record)
                del original["bt_tmin_K"]
                self.assertEqual(measured[3], C.upsilon(original, [3200.])[0])
                self.assertEqual(measured[-1], C.upsilon(original, [6000.])[0])

    def test_writer_records_only_identified_floor(self):
        levels = {1: (1., 0.), 2: (3., 100.)}
        radiation = [(1, 2, 0., 0., 1.)]
        transition = dict(ll=1, ul=2, de=.02, cups=2., ttype=5,
                          xs=np.array([0., .5, 1.]),
                          ys=np.array([.1, .2, .3]), bt_tmin_K=3200.)
        with tempfile.TemporaryDirectory() as directory:
            with (patch.object(writer, "bound_excitation_atom",
                               return_value=(levels, radiation, [transition])),
                  patch.object(writer, "_source_files", return_value=[])):
                result = writer.write_excitation_atom("c", 4, None, directory)
            contents = Path(result["path"]).read_text()
            self.assertIn("BT_TMIN_K 3200\n", contents)
            self.assertLess(contents.index("BT_TMIN_K 3200"),
                            contents.index("DENSITY_RATIO 1"))
            self.assertEqual(result["bt_tmin_K"], 3200.)
            self.assertIn("3200..32000000 K", contents)

    def test_writer_rejects_missing_identified_record(self):
        levels = {1: (1., 0.), 2: (3., 100.)}
        transition = dict(ll=1, ul=2, de=.02, cups=2., ttype=2,
                          xs=np.array([0., .5, 1.]),
                          ys=np.array([.1, .2, .3]))
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(writer, "bound_excitation_atom",
                              return_value=(levels, [(1, 2, 0., 0., 1.)], [transition])):
                with self.assertRaisesRegex(ValueError, "Inconsistent BT source"):
                    writer.write_excitation_atom("c", 4, None, directory)

    def test_rejects_nonphysical_actual_temperature(self):
        bt = dict(ll=1, ul=2, de=.02, cups=2., ttype=2,
                  xs=np.array([0., .5, 1.]), ys=np.array([.1, .2, .3]),
                  bt_tmin_K=3200.)
        table = dict(temperature_grid=np.array([1000., 2000.]),
                     upsilon_grid=np.array([.1, .2]))
        for record in (bt, table):
            for temperature in (0., -1., np.nan, np.inf, -np.inf,
                                np.array([1000., np.nan])):
                with self.subTest(table=record is table,
                                  temperature=str(temperature)):
                    with self.assertRaisesRegex(ValueError,
                                                'finite and positive'):
                        C.upsilon(record, temperature)


if __name__ == "__main__":
    unittest.main()
