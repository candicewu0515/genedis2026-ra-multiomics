"""Verify a submission without modifying it."""
from pathlib import Path
import hashlib,sys
root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
manifest=root/'MANIFEST.tsv'
if not manifest.exists():raise SystemExit('MANIFEST.tsv missing; run this inside the unpacked submission archive.')
rows=manifest.read_text().splitlines();bad=[]
for line in rows:
    name,size,want=line.split('\t');p=root/name
    if not p.is_file():bad.append(name+' missing');continue
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    if p.stat().st_size!=int(size) or h.hexdigest()!=want:bad.append(name+' checksum mismatch')
print(f'{len(rows)} manifest entries; {len(bad)} failures')
if bad:print('\n'.join(bad));raise SystemExit(1)
