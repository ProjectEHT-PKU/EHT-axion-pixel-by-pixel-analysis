import hashlib,json
from pathlib import Path
import numpy as np
root=Path(__file__).parent
original=Path('/home/canoziia/MyProject/PHYS/EHT_ZC/NormR.csv')
fresh=root/'result/result_a=0.9375_l1m1.txt'
old=np.loadtxt(original,delimiter=',');old=old[np.isclose(old[:,0],.1,rtol=0,atol=1e-12)]
new=np.loadtxt(fresh,delimiter=',')
np.testing.assert_array_equal(old[:,:2],new[:,:2])
delta=new[:,2]-old[:,2]
nonzero=old[:,2]!=0
report={'parameters':{'spin':.9375,'alpha':.1,'points':len(new),'CIRCLE_TIMES':1000},'change':'Only protect lower interpolation bracket against index -1; original real-frequency radial evaluation and normalization retained.', 'identical_numeric_values':bool(np.array_equal(old,new)), 'finite_old':bool(np.isfinite(old).all()),'finite_new':bool(np.isfinite(new).all()),'max_abs_difference':float(abs(delta).max()),'max_relative_difference_nonzero':float((abs(delta[nonzero])/abs(old[nonzero,2])).max()),'different_entries':int(np.count_nonzero(delta)),'old_peak':float(old[:,2].max()),'new_peak':float(new[:,2].max()),'peak_radius_old':float(old[np.argmax(old[:,2]),1]),'peak_radius_new':float(new[np.argmax(new[:,2]),1]),'source_csv_sha256_after':hashlib.sha256(original.read_bytes()).hexdigest()}
# Compare exact output bytes against the original slice, not entire multi-alpha file.
lines=original.read_bytes().splitlines(keepends=True)
old_bytes=b''.join(line for line in lines if line.split(b',',1)[0]==b'0.100000')
report['identical_slice_bytes']=fresh.read_bytes()==old_bytes
report['old_slice_sha256']=hashlib.sha256(old_bytes).hexdigest()
report['new_slice_sha256']=hashlib.sha256(fresh.read_bytes()).hexdigest()
(root/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
