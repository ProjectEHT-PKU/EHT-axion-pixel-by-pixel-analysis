import sys,json,tempfile
from pathlib import Path
import numpy as np
sys.path.insert(0,'/home/canoziia/MyProject/PHYS/eht-script-pipeline')
from eht_axion.pipeline.common import stage_directory,verify_stage,manifest_path
from eht_axion.pipeline.migrate import migrate_run,stage_files
r=Path(tempfile.mkdtemp(prefix='migration-',dir='/tmp/phys-independent-review'));s=r/'old';s.mkdir();raw=r/'raw.txt';raw.write_text('original')
def create(rel,inputs,writer):
    p=s/rel
    with stage_directory(p,'test',{},inputs) as t:writer(t)
    (p/'manifest.json').rename(p/'stage.json')
    return p
p=create('01_raptor/alpha-0.1',[raw],lambda t:(np.save(t/'stokes.npy',np.ones((2,17,2,2,4))), (t/'preprocess.json').write_text('{}')))
f=create('02_response/alpha-0.1',stage_files(p),lambda t:(np.save(t/'fits_1d.npy',np.ones((2,4,3))),np.save(t/'fits_2d.npy',np.ones((2,2,2,3))),(t/'fit.json').write_text('{}')))
for dim in ['1d','2d']:
    create('03_templates/'+dim,[f/f'fits_{dim}.npy'],lambda t:(np.savez(t/'templates.npz',data=[1]),(t/'summarize.json').write_text('{}')))
def eht(t):
    np.save(t/'mask.npy',np.ones((2,2),bool));np.savez(t/'selection.npz',data=[1]);(t/'report.json').write_text('{}')
    for dim,shape in [('1d',(4,4)),('2d',(4,2,2))]:
        fields={key:np.ones(shape) for key in ['mean','std','center','unwrapped_mean','concentration','circular_std','count']}
        if dim=='1d':fields['covariance']=np.ones((4,4,4))
        else:fields['intensity_mean']=np.ones(shape)
        np.savez(t/f'summary_{dim}.npz',**fields)
    for day in [95,96,100,101]:
        np.save(t/f'stokes_{day}.npy',np.ones((2,2,2,4)));np.save(t/f'evpa_1d_{day}.npy',np.ones((2,4)))
e=create('04_eht',[raw],eht)
bp=create('01_raptor_background',[raw],lambda t:(np.save(t/'stokes.npy',np.ones((2,1,2,2,4))),(t/'preprocess.json').write_text('{}')))
b=create('05_background',stage_files(bp)+[e/'mask.npy'],lambda t:((t/'background.json').write_text('{}'),np.savez(t/'background_diagnostics.npz',data=[1])))
for dim in ['1d','2d']:
    create('06_observations/'+dim,stage_files(e)+[s/'03_templates'/dim/'templates.npz'],lambda t:(np.savez(t/'observations.npz',data=[1]),(t/'observations.json').write_text('{}')))
inp=[b/'stage.json',b/'background.json']
for dim in ['1d','2d']:inp.extend([s/'06_observations'/dim/'stage.json',s/'03_templates'/dim/'stage.json',s/'06_observations'/dim/'observations.npz',s/'03_templates'/dim/'templates.npz'])
def inf(t):
    for dim in ['1d','2d']:
        for mode in ['std','cor','cor_grmhd']:
            np.savez(t/f'{mode}_{dim}.npz',alphas=[.1],cuts=[0],variant='old',observation='null',prior='log-uniform');(t/f'{mode}_{dim}.json').write_text('{}')
    (t/'run.json').write_text(json.dumps({'options':{},'results':{'std_1d':{'alphas':[.1]}}}))
i=create('07_inference',inp,inf)
# Deliberately no figures: entire migration fails late, leaves completed per-stage directories.
d=r/'partial'
try:migrate_run(s,d)
except Exception as exc:print('MIGRATION_LATE_FAILURE',type(exc).__name__,str(exc))
print('PARTIAL_DEST_EXISTS',d.exists(), 'STAGES',len(list(d.rglob('manifest.json'))))
try:migrate_run(s,d)
except Exception as exc:print('MIGRATION_RETRY',type(exc).__name__,str(exc))
create('08_figures',stage_files(i),lambda t:(t/'limits_by_dimension.png').write_bytes(b'synthetic'))
d2=r/'complete';migrate_run(s,d2)
orig=s/'02_response/alpha-0.1/stage.json';orig.write_text('invalidated original source manifest')
print('ACTUAL_MIGRATED_RESPONSE_ORIGIN_NOT_VERIFIED',verify_stage(d2/'raptor/response/alpha-0.10',recursive=True)['complete'])
print('MIGRATED_RESPONSE_INPUTS',[e['path'] for e in verify_stage(d2/'raptor/response/alpha-0.10')['inputs']])
print('ROOT',r)
