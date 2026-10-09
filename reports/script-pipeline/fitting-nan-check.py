import os,sys,json,tempfile
from pathlib import Path
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MPLBACKEND='Agg')
sys.path.insert(0,'/home/canoziia/MyProject/PHYS/eht-script-pipeline')
import numpy as np
from eht_axion.pipeline.raptor import fit_phase_dft,summarize_fits
from eht_axion.pipeline.products import load_templates
root=Path(tempfile.mkdtemp(prefix='eht-fit-nan-'));print('ROOT',root,flush=True)
psi=np.arange(16)*2*np.pi/16;X=np.column_stack([np.ones(16),np.cos(psi),np.sin(psi)])
rng=np.random.default_rng(393);rows=[]
for A in [.0,.02,.3,.75,.9,1.4]:
 y=.6+A*np.cos(psi+.4);wrapped=(y+np.pi/2)%np.pi-np.pi/2
 fitted=fit_phase_dft(wrapped[:,None])[0]
 coeff=np.linalg.lstsq(X,y,rcond=None)[0];expected=np.array([np.hypot(*coeff[1:]),np.arctan2(-coeff[2],coeff[1]),coeff[0]])
 rows.append({'trueA':A,'dft':fitted.tolist(),'unwrapped_lstsq':expected.tolist(),'residual_using_fitted':float(np.sqrt(np.mean((((wrapped-fitted[2]-fitted[0]*np.cos(psi+fitted[1])+np.pi/2)%np.pi)-np.pi/2)**2))) if A else None})
print('DFT_VS_LS',json.dumps(rows),flush=True)
for kind in ['sparse','mixed_zero','all_zero','disjoint','negative']:
 fits=np.full((141,1,3),np.nan)
 if kind=='sparse':fits[:2,0]=[[.1,.3,.2],[.2,.4,.2]]
 if kind=='mixed_zero':fits[:,0,0]=0;fits[:2,0]=[[.2,.3,.2],[.2,.3,.2]]
 if kind=='all_zero':fits[:,0,0]=0
 if kind=='disjoint':fits[:2,0,0]=[.1,.2];fits[-2:,0,1]=[.3,.4]
 if kind=='negative':fits[:2,0]=[[-.1,.3,.2],[-.2,.4,.2]]
 result=summarize_fits({.3:fits},root/kind)
 with np.load(result['templates']) as z: values={k:z[k].tolist() for k in z.files}
 _,_,_,valid=load_templates(result['templates'],'1d')
 print('SUMMARY',kind,json.dumps(values),'valid',valid.tolist(),flush=True)
# Read only saved compact products. Do not validate them as current-code caches or rerun.
run=Path('/home/canoziia/MyProject/PHYS/eht-pipeline-runs');report={}
for dim in ['1d','2d']:
 p=run/'raptor/templates'/dim/'templates.npz'
 if not p.is_file():continue
 with np.load(p,allow_pickle=False) as z:
  arrays={k:z[k] for k in z.files};alphas=arrays['alphas'];rows=[]
  for i,alpha in enumerate(alphas):
   a=arrays['amplitude_count'][i];ph=arrays['phase_count'][i];finite=np.logical_and.reduce([np.isfinite(arrays[k][i]) for k in ['amplitude_mean','amplitude_std','phase_mean','phase_std']]);row=dict(alpha=float(alpha),spatial=int(a.size),finite=int(finite.sum()),amplitude_count_range=[int(a.min()),int(a.max())],phase_count_range=[int(ph.min()),int(ph.max())],counts_differ=int((a!=ph).sum()),finite_incomplete=int((finite&((a<141)|(ph<141))).sum()))
   if finite.any():row['phase_std_range_rad']=[float(arrays['phase_std'][i][finite].min()),float(arrays['phase_std'][i][finite].max())]
   response=run/'raptor/response'/f'alpha-{alpha:.2f}'/f'response-{dim}.npz'
   if response.is_file():
    with np.load(response,allow_pickle=False) as r:
     A=r['amplitude'];P=r['phase'];joint=np.isfinite(A)&np.isfinite(P);zero=np.isfinite(A)&(A==0);row['snapshots']=len(A);row['zero_A_samples']=int(zero.sum());row['zero_A_phi_nan']=int((zero&~np.isfinite(P)).sum());row['joint_count_range']=[int(joint.sum(axis=0).min()),int(joint.sum(axis=0).max())];row['finite_with_zero_snapshots']=int((finite&zero.any(axis=0)).sum());row['finite_A_missing_phi_nonzero']=int((np.isfinite(A)&(A!=0)&~np.isfinite(P)).sum())
     # z expectation can be formed without defining phase for A=0.
     validcoef=np.isfinite(A)&((A==0)|np.isfinite(P));C=np.where(validcoef,A*np.exp(1j*np.where(np.isfinite(P),P,0)),0);n=validcoef.sum(axis=0);zmean=C.sum(axis=0)/np.maximum(n,1);zapprox=arrays['amplitude_mean'][i]*np.exp(1j*arrays['phase_mean'][i]);sel=finite&(np.abs(zapprox)>1e-15)
     if sel.any():row['complex_mean_vs_Amean_phase_relative_quantiles']=np.quantile(np.abs(zmean[sel]-zapprox[sel])/np.abs(zapprox[sel]),[.5,.9,.99,1]).tolist()
   rows.append(row);print('REAL_COMPACT',dim,json.dumps(row),flush=True)
  report[dim]=rows
(root/'real-compact-counts.json').write_text(json.dumps(report,indent=2)+'\n')
print('DONE',root,flush=True)
