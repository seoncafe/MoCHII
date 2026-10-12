"""Run focused atomic-generation checks and record a controller result."""
import argparse
import json
from pathlib import Path
import unittest

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--result', type=Path, required=True)
args = parser.parse_args()
suite = unittest.defaultTestLoader.discover(str(Path(__file__).parent))
result = unittest.TextTestRunner(verbosity=2).run(suite)
args.result.parent.mkdir(parents=True, exist_ok=True)
args.result.write_text(json.dumps(dict(passed=result.wasSuccessful(), tests=result.testsRun,
                                      failures=len(result.failures), errors=len(result.errors)), indent=2)+'\n')
raise SystemExit(0 if result.wasSuccessful() else 1)
