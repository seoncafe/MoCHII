from pathlib import Path
import hashlib,json,subprocess
ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
cmd=['make','DEBUG=1','-j4']
subprocess.run(cmd,cwd=ROOT,check=True,timeout=600)
BUILD=HERE/'build';BUILD.mkdir(exist_ok=False)
cmd2=['mpiifort','-O0','-no-ipo','-no-ftz','-check','all,noarg_temp_created','-fpe0','-traceback','-I'+str(ROOT/'src'),'-module',str(BUILD),'-o',str(BUILD/'check_source_excitation_atom.x'),str(ROOT/'tests/ion_stages/check_source_excitation_atom.f90'),str(ROOT/'src/define.o'),str(ROOT/'src/nlevel_mod.o')]
subprocess.run(cmd2,cwd=ROOT,check=True,timeout=60)
files=[ROOT/'MoCHII.x',ROOT/'src/nlevel_mod.o',ROOT/'src/nlevel_cooling_mod.o',ROOT/'src/species_mod.o',BUILD/'check_source_excitation_atom.x']
assert (ROOT/'src/nlevel_mod.o').stat().st_mtime_ns>=(ROOT/'src/nlevel_mod.f90').stat().st_mtime_ns
r={'passed':True,'command':cmd,'compile_command':cmd2,'gradual_underflow':'Focused native tests linked with -no-ftz; production defaults unchanged','products':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}}
(HERE/'build_result.json').write_text(json.dumps(r,indent=2)+'\n')
