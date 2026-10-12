from pathlib import Path
import hashlib,json,shutil,subprocess
ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
BUILD=HERE/'build'
BUILD.mkdir(exist_ok=False)
sources=['src/define.f90','src/nlevel_mod.f90','src/nlevel_cooling_mod.f90','tests/ion_stages/check_registered_excitation_atom.f90']
for name in sources:
 shutil.copyfile(ROOT/name,BUILD/Path(name).name)
flags=['mpifort','-O2','-fcheck=all','-ffpe-trap=invalid,zero,overflow','-ffree-line-length-none','-J'+str(BUILD),'-I'+str(BUILD)]
commands=[]
for name in sources[:3]:
 cmd=flags+['-c',Path(name).name];subprocess.run(cmd,cwd=BUILD,check=True);commands.append(cmd)
cmd=flags+['check_registered_excitation_atom.f90','define.o','nlevel_mod.o','nlevel_cooling_mod.o','-o','check_registered.x']
subprocess.run(cmd,cwd=BUILD,check=True);commands.append(cmd)
r={'passed':True,'commands':commands,'sources':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sources},'executable_sha256':hashlib.sha256((BUILD/'check_registered.x').read_bytes()).hexdigest()}
(HERE/'build_result.json').write_text(json.dumps(r,indent=2)+'\n')
