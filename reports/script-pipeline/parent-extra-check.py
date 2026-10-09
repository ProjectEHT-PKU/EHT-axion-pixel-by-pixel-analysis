import os,sys,runpy,tempfile
from pathlib import Path
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',MPLBACKEND='Agg',PYTHONDONTWRITEBYTECODE='1')
repo=Path('/home/canoziia/MyProject/PHYS/eht-script-pipeline');sys.path.insert(0,str(repo));root=Path(tempfile.mkdtemp(prefix='phys-parent-extra-'));print('ROOT',root,flush=True)
# Re-run source-reviewed migration reproduction in a fresh location, no prior artifacts consumed.
probe=Path(__file__).with_name('migration-reproduction.py')
if not probe.is_file():
    probe=Path('/tmp/phys-independent-review/migration_probe.py')
src=probe.read_text().replace("dir='/tmp/phys-independent-review'",'dir='+repr(str(root)))
exec(compile(src,'reviewed_migration_reproduction','exec'))
import numpy as np
from scipy.special import logsumexp
from scipy.linalg import cho_factor,cho_solve
from eht_axion.limits import LimitOptions,calculate_from_products
from eht_axion.observations import ObservationData
from eht_axion.model import evpa_templates,cloud_amplitude_scaling
# New direct dense Cholesky oracle with N=24 and condition spread.
rng=np.random.default_rng(83);n=24;q,_=np.linalg.qr(rng.normal(size=(n,n)));C=q@np.diag(np.geomspace(1e-4,2,n))@q.T
obs=ObservationData(tuple(np.zeros(n) for _ in range(4)),tuple(np.sqrt(np.diag(C/2)) for _ in range(4)),tuple(C/2 for _ in range(4)),(np.zeros(n),np.zeros(n)),(C,C*1.3),np.arange(n))
a={.3:np.vstack([rng.uniform(.01,.1,n),rng.uniform(.001,.02,n)])};ph={.3:np.vstack([rng.uniform(-np.pi,np.pi,n),rng.uniform(.01,.1,n)])};opt=LimitOptions(phase_count=8,coupling_points=9,coupling_max=50)
_,cuts,d=calculate_from_products(obs,a,ph,[.3],mode='cor_grmhd',options=opt)
s,u=evpa_templates(.3,a,ph,phase_count=8);s*=cloud_amplitude_scaling(.3);u*=cloud_amplitude_scaling(.3)
direct=[]
for g in d['couplings']:
 ll=[]
 for p in range(8):
  v=0.
  for e in range(2):
   cov=obs.pair_covariances[e]+np.diag(np.deg2rad(opt.grmhd_sigma_degrees)**2+(g*u[p,e])**2)
   L=np.linalg.cholesky(cov);w=np.linalg.solve(L,g*s[p,e]);v-=.5*(w@w+2*np.log(np.diag(L)).sum())
  ll.append(v)
 direct.append(logsumexp(ll)-np.log(8))
print('EIGEN_DIRECT_24D_MAX_LOG_ERROR',float(np.max(abs(np.array(direct)-d['log_likelihood'][0]))),flush=True)
# True complex-coefficient expectation differs from mean-amplitude/phase with broad phase.
phi=rng.normal(.4,.9,200000);print('PHASE_EXPECTATION',float(np.mean(np.sin(phi))),float(np.sin(.4)),float(np.exp(-.9**2/2)*np.sin(.4)))
# Exposure average distinction, coherent constant P amplitude.
t=np.linspace(-1.5,1.5,10001);chi=.8*np.cos(t+.4)
print('EXPOSURE_STOKES_VS_EVPA',float(np.angle(np.trapezoid(np.exp(2j*chi),t))/2),float(np.trapezoid(chi,t)/3))
# Two spatial samples with different modulation phases: amplitude map interpolation vs Stokes interpolation.
psi=np.arange(16)*2*np.pi/16;chi1=.08*np.cos(psi);chi2=-.08*np.cos(psi);complexp=(np.exp(2j*chi1)+np.exp(2j*chi2))/2;out=np.angle(complexp)/2
from eht_axion.pipeline.raptor import fit_phase_dft
print('TRANSPORT_NONCOMMUTING',float(fit_phase_dft(out)[0]),float((.08+.08)/2))
# Central-selection illustrative reduction.
from scipy.stats import norm
fraction=1572/4800;z=norm.ppf((1+fraction)/2);std=np.sqrt(1-2*z*norm.pdf(z)/fraction)
print('CENTRAL_TRUNCATION_STD_RATIO',float(std))
print('DONE',root,flush=True)
