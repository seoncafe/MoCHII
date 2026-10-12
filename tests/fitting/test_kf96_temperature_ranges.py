"""Check KF96 rate tables against independent fit evaluations."""

import importlib.util
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]


def load(name, relative_path):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


rates = load("charge_exchange_expansion_for_grid",
             "tools/fitting/charge_exchange_expansion.py")
table = load("kf96_rate_table", "tools/fitting/kf96_rate_table.py")


class Kf96TemperatureRanges(unittest.TestCase):
    def test_generator_writes_published_tables(self):
        expected = rates.kf96_existing_rates()
        seen = set()
        with tempfile.TemporaryDirectory() as directory:
            subprocess.run([sys.executable,
                            str(ROOT / "tools/fitting/make_element_data.py"),
                            "--output-dir", directory], check=True,
                           capture_output=True, text=True)
            for element in {key[0] for key in expected}:
                path = Path(directory) / f"element_{element}.txt"
                transition = None
                reaction_form = None
                for line in path.read_text().splitlines():
                    fields = line.split()
                    if not fields or fields[0].startswith("#"):
                        continue
                    if fields[0] == "TRANSITION":
                        transition = int(fields[1])
                    elif fields[0] == "CXR":
                        reaction_form = int(fields[1])
                    elif fields[0] == "CXRTAB" and (element, transition) in expected:
                        key = (element, transition)
                        self.assertEqual(reaction_form, 8)
                        count = int(fields[1])
                        self.assertEqual(len(fields), 2 + 2 * count)
                        temperatures = [float(x) for x in fields[2:2 + count]]
                        coefficients = [float(x) for x in fields[2 + count:]]
                        reaction = expected[key]
                        self.assertEqual(temperatures[0],
                                         reaction.min_temperature)
                        self.assertEqual(temperatures[-1],
                                         reaction.max_temperature)
                        for fraction in (0., .1, .5, .9, 1.):
                            temperature = reaction.min_temperature * (
                                reaction.max_temperature /
                                reaction.min_temperature) ** fraction
                            actual = table.charge_transfer_rate_on_grid(
                                temperature, temperatures, coefficients)
                            self.assertLessEqual(
                                abs(actual / reaction.rate(temperature) - 1),
                                1e-5)
                        seen.add(key)
        self.assertEqual(seen, set(expected))

    def test_all_existing_rates_and_clamps(self):
        for key, reaction in rates.published_existing_rates().items():
            with self.subTest(reaction=key):
                temperatures, coefficients = table.kf96_charge_transfer_table(
                    reaction)
                self.assertLessEqual(len(temperatures), 1024)
                self.assertEqual(temperatures[0], reaction.min_temperature)
                self.assertEqual(temperatures[-1], reaction.max_temperature)
                self.assertTrue(all(a < b for a, b in zip(
                    temperatures, temperatures[1:])))
                # The element writer uses six digits after the decimal point.
                printed_t = [float(f"{x:.6e}") for x in temperatures]
                printed_k = [float(f"{x:.6e}") for x in coefficients]
                largest_error = 0.0
                for j in range(2001):
                    temperature = reaction.min_temperature * (
                        reaction.max_temperature / reaction.min_temperature
                    ) ** (j / 2000)
                    calculated = table.charge_transfer_rate_on_grid(
                        temperature, printed_t, printed_k)
                    largest_error = max(largest_error,
                                        abs(calculated / reaction.rate(
                                            temperature) - 1))
                self.assertLessEqual(largest_error, 1e-5)
                for temperature, endpoint in (
                    (reaction.min_temperature / 10, printed_k[0]),
                    (reaction.max_temperature * 10, printed_k[-1]),
                ):
                    self.assertEqual(table.charge_transfer_rate_on_grid(
                        temperature, printed_t, printed_k), endpoint)

    def test_capacity_and_invalid_temperature(self):
        reaction = rates.published_existing_rates()["c", 2]
        with self.assertRaisesRegex(ValueError, "MAX_CX_TAB"):
            table.kf96_charge_transfer_table(reaction, max_nodes=40)
        temperatures, coefficients = table.kf96_charge_transfer_table(
            reaction)
        for bad in (0, -1, math.nan, math.inf):
            with self.subTest(temperature=bad):
                with self.assertRaises(ValueError):
                    table.charge_transfer_rate_on_grid(
                        bad, temperatures, coefficients)


if __name__ == "__main__":
    unittest.main()
