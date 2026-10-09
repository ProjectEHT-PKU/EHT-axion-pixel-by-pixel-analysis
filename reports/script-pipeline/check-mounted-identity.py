import hashlib,json,tempfile
from pathlib import Path
base=Path('/home/canoziia/MyProject/PHYS/vmdk/home/yuxin/Work/Work/EHT');root=Path(tempfile.mkdtemp(prefix='eht-mounted-identity-'));out={}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
# Mounted current-ensemble identity (stream bytes only; no pickle execution).
rows=[]
for day in [95,96,100,101]:
 p=base/'eht-imaging/result'/f'230_{day}.pkl';local=Path('/home/canoziia/MyProject/PHYS/EHT_ZC/EHT-Imaging/data_in')/p.name
 # Use curated local manifest hash as reference rather than rehashing same local bytes.
 manifest=json.loads(Path('/home/canoziia/MyProject/PHYS/EHT_ZC/shared_data/analysis-inputs/draft/manifest.json').read_text());refs=[e for e in manifest['files'] if e['path'].endswith('/'+p.name)];digest=sha(p);row={'day':day,'bytes':p.stat().st_size,'sha256':digest,'reference_candidates':[e['sha256'] for e in refs],'matches_curated_reference':any(e['sha256']==digest for e in refs)};rows.append(row);print('ENSEMBLE',json.dumps(row),flush=True)
out['ensembles']=rows
headers=[]
for day in [95,96,100,101]:
 d=base/'M87-image/eht/poloutput'/f'230_{day}';ps=sorted(d.glob('*.fits'))
 for p in ([ps[0],ps[-1]] if ps else []):
  fields={};cards=[]
  with p.open('rb') as f:
   for _ in range(50):
    block=f.read(2880)
    if not block:break
    end=False
    for i in range(0,len(block),80):
     card=block[i:i+80].decode('ascii',errors='replace');key=card[:8].strip()
     if key in ['NAXIS1','NAXIS2','CDELT1','CDELT2','CRPIX1','CRPIX2','BMAJ','BMIN','BPA','BUNIT','OBJECT','DATE-OBS','FREQ']:
      value=card[10:].split('/')[0].strip();fields[key]=value;cards.append(card.strip())
     if key=='END':end=True;break
    if end:break
  row={'path':str(p),'fits_count_same_day':len(ps),'header':fields}
  row['pixel_uas']=abs(float(fields['CDELT1']))*3600*1e6;row['FOV_uas']=row['pixel_uas']*int(fields['NAXIS1']);headers.append(row);print('FITS',json.dumps(row),flush=True)
out['headers']=headers
(root/'identity.json').write_text(json.dumps(out,indent=2)+'\n');print('DONE',root,flush=True)
