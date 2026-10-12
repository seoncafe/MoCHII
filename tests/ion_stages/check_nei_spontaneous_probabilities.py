"""Check identified Ne I decay channels and the adopted finite-node floor."""

from pathlib import Path
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'tools/fitting'))
import chianti_cooling as C  # noqa: E402
import neon_spontaneous_probabilities as Ne  # noqa: E402
import write_excitation_atom as W  # noqa: E402

stem = Path(C.ion_dir('ne', 1)) / 'ne_1'
levels = C.read_elvlc(str(stem) + '.elvlc')
source_a = C.read_wgfa(str(stem) + '.wgfa')
source_collisions = C.read_scups(str(stem) + '.scups')
assert (len(levels), len(source_a), len(source_collisions)) == (31, 142, 465)
assert sum(len(row['xs']) for row in source_collisions) == 6975
assert not any(upper == 15 for lower, upper, _, _, rate in source_a if rate > 0)

selected, radiation, collisions = C.bound_excitation_atom('ne', 1)
assert (len(selected), len(radiation), len(collisions)) == (31, 144, 465)
assert selected == {index: (g, observed) for index, (g, observed, theory)
                    in levels.items()}
assert radiation[:142] == source_a
assert radiation[142:] == Ne.nei_identified_spontaneous_probabilities(
    levels, source_a, source_collisions, stem)[142:]
assert {(lower, upper): rate for lower, upper, _, _, rate in radiation[142:]} == {
    (3, 15): 9.00e5, (5, 15): 6.82e7}
assert C.radiative_cascade_probabilities(selected, collisions, radiation)[0] == [1]
assert all(row['bt_tmin_K'] == 1000.36 for row in collisions)
for raw, adopted in zip(source_collisions, collisions):
    assert raw['ll'] == adopted['ll'] and raw['ul'] == adopted['ul']
    assert (raw['xs'] == adopted['xs']).all()
    assert (raw['ys'] == adopted['ys']).all()

saved = Ne.STOUT_SHA256
try:
    Ne.STOUT_SHA256 = ('0' * 64, '0' * 64)
    try:
        Ne.nei_identified_spontaneous_probabilities(
            levels, source_a, source_collisions, stem)
    except ValueError as error:
        assert 'Stout source identity' in str(error)
    else:
        raise AssertionError('Changed Ne I radiative source was accepted')
finally:
    Ne.STOUT_SHA256 = saved

with tempfile.TemporaryDirectory() as directory:
    result = W.write_excitation_atom('ne', 1, output_dir=directory)
    assert (result['nlev'], result['nrad'], result['nups'], result['nodes']) == (
        31, 144, 465, 6975)
    assert result['bt_tmin_K'] == 1000.36
    text = Path(result['path']).read_text()
    assert 'not a published physical temperature range' in text
    assert 'A(3,15)=900000 and A(5,15)=68200000' in text
    assert 'BT_TMIN_K 1000.3599999999999' in text or 'BT_TMIN_K 1000.36' in text

print('Ne I source identity, two restored A values, unique radiative sink, '
      'original collision nodes, and 1000.36 K floor passed')
