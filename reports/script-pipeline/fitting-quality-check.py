import os,sys,json,tempfile
from pathlib import Path
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
sys.path.insert(0,'/home/canoziia/MyProject/PHYS/eht-script-pipeline')
import numpy as np
from eht_axion.pipeline.raptor import fit_phase_dft,_evpa
run=Path('/home/canoziia/MyProject/PHYS/eht-pipeline-runs');root=Path(tempfile.mkdtemp(prefix='eht-fit-quality-'));psi=np.arange(16)*2*np.pi/16;out=[]
for alpha in [.1,.15,.2,.25,.3,.35,.4,.44]:
 p=run/'raptor/processed'/f'alpha-{alpha:.2f}'/'stokes.npy'
 if not p.is_file():continue
 cube=np.load(p,mmap_mode='r',allow_pickle=False);rows=[]
 for t in [0,len(cube)//2,len(cube)-1]:
  y=_evpa(cube[t,-16:]);vec=np.mean(np.exp(2j*y),axis=0);center=np.angle(vec)/2;lift=center+(y-center+np.pi/2)%np.pi-np.pi/2;fit=fit_phase_dft(y);model=fit[...,2][None]+fit[...,0][None]*np.cos(psi[:,None,None]+fit[...,1][None]);res=lift-model;rms=np.sqrt(np.mean(res**2,axis=0));valid=np.isfinite(fit).all(axis=-1)&(fit[...,0]>1e-12);relative=rms[valid]/fit[...,0][valid];rows.append(dict(t=int(t),finite_pixels=int(valid.sum()),amplitude_max_rad=float(np.nanmax(fit[...,0])),circular_resultant_min=float(np.nanmin(abs(vec))),residual_rms_max_rad=float(np.nanmax(rms)),residual_over_A_quantiles=np.quantile(relative,[.5,.9,.99,1]).tolist()))
 # Diagnostic across full saved time response, not actual selected EHT pixels.
 with np.load(run/'raptor/response'/f'alpha-{alpha:.2f}'/'response-2d.npz') as z:A=z['amplitude'];P=z['phase']
 center=np.angle(np.exp(1j*P).mean(axis=0));dp=(P-center+np.pi)%(2*np.pi)-np.pi;da=A-A.mean(axis=0);dp-=dp.mean(axis=0);cov=(da*dp).mean(axis=0);sigA=da.std(axis=0);sigP=dp.std(axis=0);corr=np.divide(cov,sigA*sigP,out=np.full(cov.shape,np.nan),where=(sigA*sigP)>0);scale=32.6868/64;xy=(np.arange(64)-31.5)*scale;rad=np.hypot(xy[:,None],xy[None,:]);geom=(rad>=2)&(rad<=9)
 record=dict(alpha=alpha,phasefit=rows,full_grid_abs_Aphi_corr_quantiles=np.nanquantile(abs(corr),[.5,.9,.99,1]).tolist(),geometric_annulus_only={'definition':'2<=screen r<=9Rg, not EHT selected mask','pixels':int(geom.sum()),'phase_std_gt_0p5':int((geom&(sigP>.5)).sum()),'phase_std_gt_1':int((geom&(sigP>1)).sum()),'max_phase_std':float(sigP[geom].max()),'abs_Aphi_corr_quantiles':np.nanquantile(abs(corr[geom]),[.5,.9,.99,1]).tolist()});out.append(record);print(json.dumps(record),flush=True)
(root/'quality.json').write_text(json.dumps(out,indent=2)+'\n');print('DONE',root,flush=True)
