import pickle,json,tempfile
from pathlib import Path
import numpy as np
base=Path('/home/canoziia/MyProject/PHYS/vmdk/home/yuxin/Work/Work/EHT');root=Path(tempfile.mkdtemp(prefix='eht-fits-pickle-'))
# Trusted local historical pickle, loaded once in isolated process (~630 MB).
with (base/'eht-imaging/result/230_95.pkl').open('rb') as f:d=pickle.load(f)
def read_hdus(path):
 out=[]
 with path.open('rb') as f:
  while True:
   header={};done=False
   for _ in range(100):
    block=f.read(2880)
    if not block:return out
    for i in range(0,2880,80):
     card=block[i:i+80].decode('ascii');key=card[:8].strip()
     if key=='END':done=True;break
     if card[8:10]=='= ':header[key]=card[10:].split('/')[0].strip()
    if done:break
   naxis=int(header.get('NAXIS','0'));shape=[int(header[f'NAXIS{i}']) for i in range(1,naxis+1)];count=int(np.prod(shape)) if shape else 0;bits=int(header['BITPIX']);nbytes=count*abs(bits)//8;data=f.read(nbytes);f.seek((-nbytes)%2880,1)
   dtype={-64:'>f8',-32:'>f4',8:'u1',16:'>i2',32:'>i4',64:'>i8'}[bits];array=np.frombuffer(data,dtype=dtype).reshape(tuple(reversed(shape))) if shape else np.array([]);out.append((header,array))
results=[]
for k in [0,4799]:
 m=d['meta_data'][k];name='_'.join(str(float(m[x])) for x in ['zbl','prior_fwhm','sys_noise','MEM','TV','T2V','l1'])+'.fits';p=base/'M87-image/eht/poloutput/230_95'/name;hdus=read_hdus(p)
 row={'candidate':k,'path':str(p),'hdus':[{'stokes':h.get('STOKES'),'extname':h.get('EXTNAME'),'shape':a.shape} for h,a in hdus],'IQU_exact_match':{key:bool(np.array_equal(d[key][...,k],hdus[i][1])) for i,key in enumerate(['I_data','Q_data','U_data'])},'theta_matches_V':bool(np.array_equal(d['theta_data'][...,k],hdus[3][1])),'theta_maxabs':float(np.max(abs(d['theta_data'][...,k])))};results.append(row);print(json.dumps(row),flush=True)
(root/'fits-pickle.json').write_text(json.dumps(results,indent=2)+'\n');print('DONE',root,flush=True)
