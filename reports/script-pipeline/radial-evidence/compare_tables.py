import json,hashlib
from pathlib import Path
import numpy as np
root=Path('/tmp/eht-radial-audit');one=root/'archive1/Radial_function/Radial_function/Radial_function';two=root/'archive2/Radial_function';mounted=Path('/home/canoziia/MyProject/PHYS/EHT_ZC')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
records=[]
files=[*one.glob('result/*.txt*'),one/'NormR.csv',*two.glob('result/*.txt*'),two/'NormR.csv',mounted/'NormR.csv']
for p in files:
 x=np.loadtxt(p,delimiter=',');a=np.unique(x[:,0]);r=np.unique(x[:,1]);groups=[]
 for val in a:
  y=x[x[:,0]==val];groups.append({'alpha':float(val),'N':len(y),'Rmax':float(np.nanmax(y[:,2])),'negative':int((y[:,2]<0).sum()),'nonfinite':int((~np.isfinite(y[:,2])).sum())})
 records.append({'path':str(p),'sha256':sha(p),'shape':x.shape,'alpha_range':[float(a.min()),float(a.max())],'alpha_count':len(a),'radial_range':[float(r.min()),float(r.max())],'radial_count':len(r),'groups':groups});print(json.dumps({k:v for k,v in records[-1].items() if k!='groups'}),flush=True)
for p in [mounted/'NormR.csv']:
 matches=[r['path'] for r in records if r['path']!=str(p) and r['sha256']==sha(p)];print('MOUNTED_EXACT_MATCH',str(p),matches,flush=True)
# Compare same-a old result table numeric equivalence and C file source identities.
p=two/'result/result_a=0.9375_l1m1.txt';q=two/'result_old/result_a=0.9375_l1m1.csv';print('OLD_RESULT_SAME_2021',sha(p)==sha(q),np.array_equal(np.loadtxt(p,delimiter=','),np.loadtxt(q,delimiter=',')),flush=True)
for filename in ['Axion_Radial.h','README.txt','Interpolation.nb','extraporlation.nb','Plot_Radial_funtion.nb']:
 print('IDENTICAL_SOURCE',filename,sha(one/filename)==sha(two/filename),flush=True)
# get_omega implements interpolation of omega/mu by mu grid, then multiplies by mu.
freq=[]
for spin in [.9375,.99]:
 re=np.loadtxt(one/f'source/a={spin}_l1m1_real.txt',skiprows=1,delimiter=',');im=np.loadtxt(one/f'source/a={spin}_l1m1_ima.txt',skiprows=1,delimiter=',')
 for alpha in [.1,.15,.2,.25,.3,.35,.4,.44]:
  if not(re[0,0]<alpha<=re[-1,0] and im[0,0]<alpha<=im[-1,0]):continue
  w=float(np.interp(alpha,re[:,0],re[:,1])*alpha);wi=float(np.interp(alpha,im[:,0],im[:,1])*alpha);w4=alpha*(1-alpha**2/32);w2=alpha*(1-alpha**2/8)
  freq.append({'spin':spin,'alpha':alpha,'omega_re_interpolated':w,'omega_im_interpolated':wi,'omega4':w4,'omega2':w2,'omega4_rel_error_percent':(w4/w-1)*100,'omega2_rel_error_percent':(w2/w-1)*100,'direct_alpha_rel_error_percent':(alpha/w-1)*100});print('FREQUENCY',json.dumps(freq[-1]),flush=True)
(root/'tables-comparison.json').write_text(json.dumps(records,indent=2)+'\n');(root/'frequency-comparison.json').write_text(json.dumps(freq,indent=2)+'\n')
# Inventory all science files differing, not bundled environments.
left={str(p.relative_to(one)):p for p in one.rglob('*') if p.is_file()};right={str(p.relative_to(two)):p for p in two.rglob('*') if p.is_file()};diff=[]
for name in sorted(left.keys()|right.keys()):
 if name not in left:status='only_archive2'
 elif name not in right:status='only_archive1'
 elif sha(left[name])!=sha(right[name]):status='different'
 else:continue
 diff.append({'path':name,'status':status})
print('DIFFERENCES',json.dumps(diff),flush=True);(root/'file-differences.json').write_text(json.dumps(diff,indent=2)+'\n')
