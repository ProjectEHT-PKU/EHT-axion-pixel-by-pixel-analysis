import json,hashlib,tempfile
from pathlib import Path
root=Path(tempfile.mkdtemp(prefix='eht-mounted-cube-hashes-'));base=Path('/home/canoziia/MyProject/PHYS/vmdk/home/yuxin/Work/Work/EHT/EHT-axion-pixel-by-pixel/raptor_output');manifest=json.loads(Path('/home/canoziia/MyProject/PHYS/EHT_ZC/shared_data/analysis-inputs/draft/manifest.json').read_text());rows=[]
for alpha in [.1,.15,.2,.25,.3,.35,.4,.44]:
 name=f'GRMHD_Data_1600-3000_alpha={alpha}_res=100.npy';p=base/name;h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 refs=[e['sha256'] for e in manifest['files'] if e['path'].endswith(f'/alpha-{alpha:g}.npy')];row={'alpha':alpha,'filename':name,'bytes':p.stat().st_size,'sha256':h.hexdigest(),'reference_candidates':refs,'matches_curated_reference':h.hexdigest() in refs};rows.append(row);print(json.dumps(row),flush=True)
(root/'cube-hashes.json').write_text(json.dumps(rows,indent=2)+'\n');print('DONE',root,flush=True)
