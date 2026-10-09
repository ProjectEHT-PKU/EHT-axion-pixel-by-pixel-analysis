import os,sys,json,struct,tempfile,math
from pathlib import Path
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
sys.path.insert(0,'/home/canoziia/MyProject/PHYS/eht-script-pipeline')
import numpy as np,h5py
from eht_axion.pipeline.constants import PhysicalScale,HBAR_EV_S
base=Path('/home/canoziia/MyProject/PHYS/vmdk/home/yuxin/Work/Work/EHT');new=base/'EHT-axion-pixel-by-pixel';root=Path(tempfile.mkdtemp(prefix='eht-mounted-constants-'));report={}
oldG=6.674e-8;oldM=1.989e33;oldc=2.99792458e10;oldpc=3.086e18
s=PhysicalScale(mass_solar=6.2e9,distance_mpc=16.8);current=PhysicalScale()
rg=oldG*oldM*6.2e9/oldc**2;tg=rg/oldc;theta=206265*rg/(16800*3.086e21)*1e6
k=np.sqrt(8*np.log(2));
report['scales']={'producer_old_constants':{'M':6.2e9,'D_Mpc':16.8,'rg_cm':rg,'tg_h':tg/3600,'theta_uas_per_rg':theta,'omega_per_hour':3600/tg,'mass_per_alpha_eV':HBAR_EV_S/tg},'producer_M_modern_units':s.metadata(),'current':current.metadata()}
report['geometry']={'nativeFOV_rg':46.2261,'nativeFOV_uas':46.2261*theta,'native_half_FOV_uas':46.2261*theta/2,'current_finalFOV_uas_producer':32.6868*theta,'current_sigmaFWHM_uas_producer':2.16871*k*theta,'matching128_rg_oldunits':128/theta,'matching120_rg_oldunits':120/theta,'matching20sigma_rg_oldunits':20/(k*theta),'matching19sigma_rg_oldunits':19/(k*theta),'matching128_rg_modernunits':128/s.angular_rg_uas,'matching20sigma_rg_modernunits':20/(k*s.angular_rg_uas),'min_nativeFOV_for_rot108_and128_rg':128/theta*(abs(np.cos(np.deg2rad(108)))+abs(np.sin(np.deg2rad(108)))),'native_div_currentfinal':46.2261/32.6868}
report['relative_current_vs_producer_percent']={'mass':(current.mass_solar/6.2e9-1)*100,'tg':(current.tg_hours/(tg/3600)-1)*100,'theta':(current.angular_rg_uas/theta-1)*100,'omega_conversion':((1/current.tg_hours)/(3600/tg)-1)*100,'mass_per_alpha':((HBAR_EV_S/current.tg_seconds)/(HBAR_EV_S/tg)-1)*100,'finalFOV_vs128':(32.6868/(128/theta)-1)*100,'sigma_vs20':(2.16871/(20/k/theta)-1)*100}
# Fluid tail according to actual producer reader, only final ~100 bytes.
fluid=base/'raptor_run/MAD_a+15o16/a094-high';tails=[]
for label in [1500,1501,1600,1610,3000]:
 path=fluid/f'data{label:06d}.dat'
 if not path.is_file():continue
 with path.open('rb') as f:
  f.seek(-40,2);v=struct.unpack('<8Id',f.read(40));ndim,neq=v[2],v[6];f.seek(-40-(ndim*4+neq*8),2);dims=struct.unpack('<'+'I'*ndim,f.read(ndim*4));params=struct.unpack('<'+'d'*neq,f.read(neq*8))
 tails.append({'label':label,'time_tg':v[-1],'grid':dims,'neqpar':params,'spin_NSPIN3':params[3]})
report['fluid_tails']=tails
report['timing']={'producer_3tg_h':3*tg/3600,'modern6p2_3tg_h':3*s.tg_hours,'current3tg_h':3*current.tg_hours,'producer_span_140tg_days':140*tg/86400}
# Check modest samples, not full 6GB hash scan.
comparisons=[]
for alpha in [.1,.15,.2,.25,.3,.35,.4,.44]:
 name=f'GRMHD_Data_1600-3000_alpha={alpha}_res=100.npy';p=new/'raptor_output'/name;local=Path('/home/canoziia/MyProject/PHYS/EHT_ZC/GRMHD/data_in')/name
 a=np.load(p,mmap_mode='r',allow_pickle=False);b=np.load(local,mmap_mode='r',allow_pickle=False);items=[]
 for t in [0,70,140]:
  for state in [0,1,16]:items.append(bool(np.array_equal(a[t,state],b[t,state])))
 comparisons.append({'alpha':alpha,'shape':a.shape,'sampled_9_frames_all_equal':all(items)})
report['cube_comparisons']=comparisons
p=new/'raptor_output/0.3/data001600.dat/img_data_0.h5'
with h5py.File(p,'r') as f:
 x=f['alpha'][()];y=f['beta'][()];ux=np.unique(x);uy=np.unique(y);output=np.empty((100,100,4));names=['I2.300000e+11','Q2.300000e+11','U2.300000e+11','V2.300000e+11']
 # Direct producer block/pixel layout reconstruction.
 B=int(np.sqrt(x.shape[0]));P=int(np.sqrt(x.shape[1]));bx,by=np.indices((B,B));px,py=np.indices((P,P));ix=(bx.ravel()[:,None]*P+px.ravel()).ravel();iy=(by.ravel()[:,None]*P+py.ravel()).ravel()
 for j,key in enumerate(names):output[ix,iy,j]=f[key][()].ravel()
 local=np.load('/home/canoziia/MyProject/PHYS/EHT_ZC/GRMHD/data_in/GRMHD_Data_1600-3000_alpha=0.3_res=100.npy',mmap_mode='r')
 report['HDF5_match']={'path':str(p),'coordinate_minmax':[float(ux.min()),float(ux.max())],'pixelstep':float(np.median(np.diff(ux))),'edgeFOV_rg':float(ux.ptp()+np.median(np.diff(ux))) if hasattr(ux,'ptp') else float(np.ptp(ux)+np.median(np.diff(ux))),'IQUV_exact_match_current_cube_t1600_state0':bool(np.array_equal(output,local[0,0]))}
# Small FITS header blocks, no image loading.
fits_files=[base/'eht-imaging/out.fits',base/'M87-image/eht/data/M87_2017_095_lo_hops_netcal_StokesI.fits']
headers=[]
for p in fits_files:
 if not p.is_file():continue
 cards=[]
 with p.open('rb') as f:
  for _ in range(50):
   block=f.read(2880)
   if not block:break
   done=False
   for i in range(0,len(block),80):
    card=block[i:i+80].decode('ascii',errors='replace');key=card[:8].strip()
    if key in ['NAXIS','NAXIS1','NAXIS2','CDELT1','CDELT2','CRPIX1','CRPIX2','CTYPE1','CTYPE2','BMAJ','BMIN','BPA','BUNIT','OBJECT']:cards.append(card.strip())
    if key=='END':done=True;break
   if done:break
 headers.append({'path':str(p),'cards':cards})
report['sample_FITS_headers']=headers
(root/'comparison.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True);print('DONE',root,flush=True)
