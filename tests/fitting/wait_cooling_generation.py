"""Wait for a declared cooling controller and verify complete staged products."""
import argparse
import hashlib
import json
from pathlib import Path
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--state', type=Path, required=True)
parser.add_argument('--products', type=Path, required=True)
parser.add_argument('--result', type=Path, required=True)
args = parser.parse_args()
deadline = time.monotonic()+900
while True:
    state = json.loads(args.state.read_text())
    if state['status'] == 'failed':
        raise RuntimeError('Source cooling controller failed; preserve its products')
    if state['status'] == 'completed':
        break
    if time.monotonic() > deadline:
        raise TimeoutError('Cooling completion deadline exceeded')
    time.sleep(2)
products = {}
for name,el,stage in [('CIII','c',3),('SiIII','si',3),('SV','s',5)]:
    check = json.loads((args.products/name/f'{name}.json').read_text())
    path = args.products/name/f'cooling_{el}_{stage}.txt'
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if not check['passed'] or check['product_sha256'] != digest:
        raise RuntimeError('Changed or invalid staged cooling file: '+str(path))
    if not (args.products/name/f'{name}.pdf').is_file():
        raise RuntimeError('Missing complete temperature curve: '+name)
    products[name] = digest
args.result.parent.mkdir(parents=True,exist_ok=True)
args.result.write_text(json.dumps(dict(passed=True,products=products),indent=2)+'\n')
