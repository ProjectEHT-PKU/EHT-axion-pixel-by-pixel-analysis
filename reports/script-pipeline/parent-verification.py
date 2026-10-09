from pathlib import Path
import os,sys,tempfile,json,contextlib,io,importlib.util,shutil
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MPLBACKEND='Agg',PYTHONDONTWRITEBYTECODE='1')
repo=Path('/home/canoziia/MyProject/PHYS/eht-script-pipeline');sys.path.insert(0,str(repo))
root=Path(tempfile.mkdtemp(prefix='phys-parent-check-')); print('ROOT',root,flush=True)
import numpy as np
from eht_axion.pipeline import common,workflow,layout,raptor,eht,convergence,inference
from eht_axion.polarization import integrate_stokes
from eht_axion.limits import LimitOptions,calculate_from_products
from eht_axion.observations import ObservationData
from eht_axion.model import omega,evpa_templates,cloud_amplitude_scaling
from scipy.special import i0e
from scipy.optimize import brentq
from eht_axion.statistics import posterior_upper_limit
results={}
def record(k,v):results[k]=v;print(k,json.dumps(v,default=str),flush=True)
# Independent boundary probes, no actual input files.
cfg=root/'bad.yaml';cfg.write_text((repo/'configs/research.yaml').read_text().replace('historical_mask: true','historical_mask: "false"').replace('inherit_2d_mask_in_1d: false','inherit_2d_mask_in_1d: "false"'))
c=workflow.read_config(cfg);record('B3_yaml',{k:[type(c['eht'][k]).__name__,bool(c['eht'][k])] for k in ('historical_mask','inherit_2d_mask_in_1d')})
x=np.array([[0,-.1],[np.pi/2,.2],[0,.1],[np.pi/2,-.05]])
r=[]
for eps in (0,1e-12):
 y=x.copy();y[:,0]+=eps;s=eht.summarize_evpa(y);r.append(dict(eps=eps,center=s['center'].tolist(),cross=float(eht.ensemble_covariance(y,s)[0,1]),resultant=s['concentration'].tolist()))
record('B4_circle',r)
f=np.ones((2,2,4));f[...,1]=np.inf;record('B9_inf_evpa',raptor._evpa(f).tolist())
f=np.ones((2,2,4));f[...,1]=np.nan;record('B9_nan_integral',integrate_stokes(f,np.ones((1,4))).tolist())
f=np.zeros((2,2,4),dtype=np.uint8);f[...,0]=200;f[...,1]=2;f[...,2]=1
record('B9_integer_integral',[float(integrate_stokes(f,np.ones((1,4)))[0]),float(integrate_stokes(f.astype(float),np.ones((1,4)))[0])])
raw=root/'raw.txt';raw.write_text('v1');p=root/'proc';r=root/'resp';t=root/'template'
with common.stage_directory(p,'test',{},[raw]) as tmp:(tmp/'x.txt').write_text('p')
with common.stage_directory(r,'test',{},[p/'manifest.json',p/'x.txt'],imported_from={'path':'old'}) as tmp:np.savez(tmp/'response.npz',x=[1])
with common.stage_directory(t,'test',{},[r/'response.npz']) as tmp:(tmp/'template.txt').write_text('t')
raw.write_text('v2')
record('B1_recursive_gap',common.verify_stage(t,recursive=True,allow_imported=False)['complete'])
try:common.verify_stage(r,recursive=True)
except Exception as ex:record('B1_direct_rejection',str(ex))
external=root/'external';external.mkdir();alias=root/'stage-link';alias.symlink_to(external/'missing',target_is_directory=True)
with common.stage_directory(alias,'test',{},[]) as tmp:(tmp/'x.txt').write_text('x')
record('B2_dangling_stage',(external/'missing/x.txt').is_file())
out=root/'lowlevel';out.mkdir();(out/'stokes.npy').symlink_to(external/'missing.npy')
raptor.preprocess_stokes(np.ones((1,2,2,4)),out,sigma=0,init_fov=2,final_fov=2,rotation_angle=0,final_shape=2)
record('B2_dangling_file',(external/'missing.npy').is_file())
shardroot=root/'shards';d=shardroot/'raptor/res-400/alpha-0.30';d.mkdir(parents=True)
for name in ('stokes-t900-t990.npy','stokes-t1000-t1090.npy'):np.save(d/name,np.ones((1,1,2,2,4)))
record('B5_order',[p.name for p in layout.input_raptor_files(shardroot,.3,400)])
for p in d.glob('stokes-t*'):p.unlink()
np.save(d/'stokes.npy',np.ones((1,1,2,2,4)));legacy=shardroot/'raptor-stokes/res-400/alpha-0.3';legacy.mkdir(parents=True);np.save(legacy/'stokes-t1600-t1690.npy',np.ones((1,1,2,2,4)))
record('B5_priority',[str(p.relative_to(shardroot)) for p in layout.input_raptor_files(shardroot,.3,400)])
obs=ObservationData(tuple(np.zeros(1) for _ in range(4)),tuple(np.ones(1) for _ in range(4)),tuple(np.eye(1)*.5 for _ in range(4)),(np.zeros(1),np.zeros(1)),(np.eye(1),np.eye(1)),np.array([0]))
a={.3:np.array([[0.],[0.]])};ph={.3:np.array([[0.],[0.]])}
record('B10_legacy_prior',{prior:float(calculate_from_products(obs,a,ph,[.3],options=LimitOptions(legacy=True,prior=prior,coupling_min=1,coupling_max=101,coupling_points=9))[1][0]) for prior in ('uniform','log-uniform')})
# Create a fresh tiny chain independently of subagent's mutated artifacts.
spec=importlib.util.spec_from_file_location('fixture',repo/'tests/test_pipeline_end_to_end.py');fixture=importlib.util.module_from_spec(spec);spec.loader.exec_module(fixture)
inputs=root/'input';fixture.make_inputs(inputs)
text=(repo/'configs/research.yaml').read_text()
for before,after in [('alphas: [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.44]','alphas: [0.1, 0.3]'),('final_shape: 64','final_shape: 8'),('topset_size: 1572','topset_size: 40'),('min_radius_rg: 2.0','min_radius_rg: 0.0'),('historical_mask: true','historical_mask: false'),('covariance_ridge: 0.0','covariance_ridge: 0.001'),('source_resolution: 400','source_resolution: 100'),('phase_count: 64','phase_count: 16'),('coupling_points: 120','coupling_points: 30')]:text=text.replace(before,after)
cfg=root/'config.yaml';cfg.write_text(text);out=root/'output';workflow.run_stage(cfg,out,input_root=inputs,trust_pickle=True);config=workflow.read_config(cfg)
record('fresh_chain',str(out))
raw=inputs/'raptor-stokes/res-100/alpha-0.1.npy';before=raw.read_bytes();v=np.load(raw);v[...,1]*=.91;np.save(raw,v)
try:inference.run_inference(out,config,rebuild=True);record('B1_actual_stale_raw_accepted',True)
except Exception as ex:record('B1_actual_stale_raw_accepted',str(ex))
finally:raw.write_bytes(before)
r=out/'raptor/response/alpha-0.10/manifest.json';before=r.read_bytes();v=json.loads(before);v['imported_from']={'path':'external-origin'};r.write_text(json.dumps(v))
try:inference.run_inference(out,config,rebuild=True);record('B1_actual_imported_accepted',True)
except Exception as ex:record('B1_actual_imported_accepted',str(ex))
finally:r.write_bytes(before)
# Ensure original fine reference re-generated before convergence.
inference.run_inference(out,config,rebuild=True)
t=out/'raptor/templates/1d/templates.npz';before=t.read_bytes()
with np.load(t) as z:v={k:z[k] for k in z.files}
v['amplitude_mean']*=2;np.savez(t,**v)
try:
 report=convergence.check_convergence(out,root/'convergence.json',dimensions=('1d',),modes=('std',));record('B6_stale_convergence_accepted',report['checks'])
finally:t.write_bytes(before)
alias=root/'run-alias';alias.symlink_to(out,target_is_directory=True);before=(out/'figures/manifest.json').stat().st_mtime_ns
workflow.run_stage(cfg,alias,stage='figures',rebuild=True)
record('B2_root_symlink_rebuilt',(out/'figures/manifest.json').stat().st_mtime_ns!=before)
# Own numerical checks.
alpha=brentq(lambda x:omega(x)*.112411*120-np.pi,.2,.3);a={alpha:np.array([[1.],[0.]])};ph={alpha:np.array([[.037],[0.]])}
grid={}
for n in (64,128):
 _,cut,diag=calculate_from_products(obs,a,ph,[alpha],options=LimitOptions(coupling_points=120,phase_count=n));grid[str(n)]=float(cut[0])
s,_=evpa_templates(alpha,a,ph,phase_count=4);s=s.reshape(4,-1)*cloud_amplitude_scaling(alpha);H=np.array([[s[0]@s[0],s[0]@s[1]],[s[0]@s[1],s[1]@s[1]]]);w=np.maximum(np.linalg.eigvalsh(H),0);g=np.array(diag['couplings']);exact=-g*g*w[0]/2+np.log(i0e(g*g*(w[1]-w[0])/4));grid['exact']=float(np.log10(posterior_upper_limit(g,exact)));record('S8_phasegrid',dict(alpha=alpha,**grid))
from eht_axion.pipeline.transport import legacy_rms_amplitude
psi=np.arange(16)*2*np.pi/16
record('transport_rms_bias',[float(legacy_rms_amplitude(.1*np.cos(psi+x))/.1) for x in (0,np.pi/2)])
record('superradiance',{'OmegaH':.9375/(2*(1+np.sqrt(1-.9375**2))),'omega':[omega(x) for x in (.35,.4,.44)]})
(root/'results.json').write_text(json.dumps(results,indent=2,default=str)+'\n')
print('DONE',root,flush=True)
