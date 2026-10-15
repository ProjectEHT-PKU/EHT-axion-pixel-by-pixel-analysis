from pathlib import Path
import ast,subprocess,tempfile,os,shutil,re,json
BASE=Path(tempfile.mkdtemp(prefix='raptor-build-run-probe-'))
print('PROBE_DIR',BASE)
core=Path('/tmp/eht-github-raptor-GckNId/raptor')
run=Path('/tmp/eht-github-raptor-GckNId/run/raptor_run')
# Reconstruct ONLY model.in f.write expressions, without importing/executing Run.py.
tree=ast.parse((run/'Run.py').read_text())
fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='Run')
params={'MBH':6.2e9,'DISTANCE':16800,'M_UNIT':3e25,'R_LOW':1,'R_HIGH':1,'INCLINATION':163,'IMG_WIDTH':100,'IMG_HEIGHT':100,'CAM_SIZE_X':46.2261,'CAM_SIZE_Y':46.2261,'FREQS_PER_DEC':1,'FREQ_MIN':230e9,'STEPSIZE':0.005,'AXION_NORM':0.1,'AXION_OMEGA':0.3,'AXION_PHASE':15,'MAX_LEVEL':1}
text=''
for n in fn.body:
 if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Attribute) and isinstance(n.value.func.value,ast.Name) and n.value.func.value.id=='f' and n.value.func.attr=='write':
  text+=eval(compile(ast.Expression(n.value.args[0]),'<extracted-write-expression>','eval'),{'__builtins__':{}},{'params_json':params})
(BASE/'generated-model.in').write_text(text)
print('model.in newline count',text.count('\n'),'bytes',len(text))
# Use actual read_model function, but no model init/ray tracing.
s=(core/'src/core.c').read_text();a=s.index('void read_model(');b=s.index('// For a single block',a)
body=s[a:b]
c='''#include <stdio.h>
#include <stdlib.h>
#define MSUN 1.989e33
#define KPCTOCM 3.086e21
#define MKSBHAC 6
#define metric MKSBHAC
char GRMHD_FILE[256];
double MBH,M_UNIT,TIME_INIT,INCLINATION,axion_norm,axion_omega,phase_index,R_HIGH,R_LOW,FREQS_PER_DEC,FREQ_MIN,SOURCE_DIST,CAM_SIZE_X,CAM_SIZE_Y,STEPSIZE;
int IMG_WIDTH,IMG_HEIGHT,max_level;
'''+body+'''int main(int argc,char **argv){read_model(argv);printf("RESULT norm=%g omega=%g phase=%g max_level=%d width=%d output_index=%g\\n",axion_norm,axion_omega,phase_index,max_level,IMG_WIDTH,TIME_INIT);return 0;}
'''
(BASE/'parse-only.c').write_text(c)
subprocess.run(['gcc','-O0',str(BASE/'parse-only.c'),'-o',str(BASE/'parse-only')],check=True)
for f in [BASE/'generated-model.in',run/'model.in']:
 p=subprocess.run([str(BASE/'parse-only'),str(f),'snapshot.dat','16'],text=True,capture_output=True)
 print('PARSE',str(f),p.stdout.strip());print(p.stderr.strip())
# GCC preprocessor only: provide an empty hdf5.h stub, NOT a scientific compile.
(BASE/'stubs').mkdir();(BASE/'stubs/hdf5.h').write_text('/* audit preprocessing stub only */\n')
local=Path('/home/canoziia/MyProject/PHYS/vmdk/home/yuxin/Work/Work/EHT/EHT-axion-pixel-by-pixel')
for file in [local/'raptor/src/core.c',local/'raptor_run/model.c']:
 args=['gcc','-E','-H','-dM','-I'+str(local/'raptor/src'),'-I'+str(local/'raptor_run'),'-I'+str(BASE/'stubs'),str(file)]
 p=subprocess.run(args,text=True,capture_output=True)
 print('PREPROCESS',file,'rc',p.returncode)
 for l in p.stderr.splitlines():
  if str(local) in l: print(l)
 for l in p.stdout.splitlines():
  if re.search(r'^#define (SIGMA_CUT|metric |POL |DEBUG |LIGHT_TRANSPORT|RAD_TRANS)',l): print(l)
# Reuse original makefile in isolated toy source tree to test VPATH and stale headers.
source=BASE/'toy-core/src';source.mkdir(parents=True)
work=BASE/'toy-run';work.mkdir();(work/'build').mkdir()
shutil.copyfile(run/'makefile',work/'makefile')
names=re.search(r'^SOURCES=(.*)$',(run/'makefile').read_text(),re.M)[1].split()
for name in names:
 (source/name).write_text('/* remote source */\n')
 (work/'build'/Path(name).with_suffix('.o')).touch()
(work/'model.c').write_text('/* local override */\n')
(work/'definitions.h').write_text('#define SIGMA_CUT 1000\n')
(work/'RAPTOR').touch()
# timestamps older source, object/binary fixed; header later: make must still report no rebuild.
for p in source.iterdir():os.utime(p,(100,100))
os.utime(work/'model.c',(100,100))
for p in (work/'build').iterdir():os.utime(p,(200,200))
os.utime(work/'RAPTOR',(300,300));os.utime(work/'definitions.h',(400,400))
env=os.environ.copy();env.update(RAPTOR=str(source.parent),PWD=str(work))
for args in [['make','-n','all'],['make','-n','-B','all'],['make','-n','-B','all','RAPTOR='+str(BASE/'different-core')]]:
 p=subprocess.run(args,cwd=work,env=env,text=True,capture_output=True)
 print('MAKE_PROBE',args,'rc',p.returncode); print(p.stdout.strip());print(p.stderr.strip())
print('ALL_PROBES_ONLY_IN',BASE)
