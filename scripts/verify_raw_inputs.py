"""Verify known primary raw inputs and record their hashes without silently replacing mismatches."""
from pathlib import Path
import hashlib,json,gzip,sys
ROOT=Path(__file__).resolve().parents[1];raw=ROOT/'data/raw';manifest=raw/'SHA256SUMS.txt'
expected={}
if manifest.exists():
 expected={p.removeprefix('./'):h for h,p in (x.split(maxsplit=1) for x in manifest.read_text().splitlines())}
def exclude(name):
 return (name=='SHA256SUMS.txt' or '/cache_' in '/'+name or
   any(name.startswith(r) for r in ['opentargets/','ensembl/','collectri/','ldsc/GRCh38/','gse296117/converted/']) or Path(name).name=='.DS_Store')
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()
rows={};bad=[]
for p in sorted(raw.rglob('*')):
 if not p.is_file():continue
 name=str(p.relative_to(raw))
 if exclude(name):continue
 h=digest(p);rows[name]=h
 if name in expected and h!=expected[name]:bad.append(name+' differs from submitted input')
for name in expected:
 if not exclude(name) and name not in rows:bad.append(name+' missing')
if bad:raise SystemExit('Raw input verification failed:\n'+'\n'.join(bad))
manifest.write_text(''.join(h+'  ./'+name+'\n' for name,h in rows.items()))
print('PASS:',len(rows),'primary raw/reference files; matching submitted hashes where present')
