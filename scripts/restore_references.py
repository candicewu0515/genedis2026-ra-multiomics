"""Restore the saved reference snapshot and check every file in it."""
from pathlib import Path
import hashlib,json,tarfile,urllib.request,os
ROOT=Path(__file__).resolve().parents[1]
ref=ROOT/'references'
archive=ref/'reference_snapshot.tar.gz'
expected_archive=(ROOT/'env/reference_snapshot.sha256').read_text().split()[0]
if not archive.exists():
    url='https://github.com/candicewu0515/genedis2026-ra-multiomics/releases/download/v1.4-submission/reference_snapshot.tar.gz'
    ref.mkdir(parents=True,exist_ok=True)
    temporary=archive.with_suffix('.download')
    print('Downloading saved reference snapshot from',url,flush=True)
    try:
        urllib.request.urlretrieve(url,temporary)
        if hashlib.sha256(temporary.read_bytes()).hexdigest()!=expected_archive:
            raise RuntimeError('Downloaded reference archive checksum mismatch')
        os.replace(temporary,archive)
    finally:
        temporary.unlink(missing_ok=True)
if hashlib.sha256(archive.read_bytes()).hexdigest()!=expected_archive:
    raise RuntimeError('Reference archive checksum mismatch')
manifest=json.loads((ref/'snapshot_manifest.json').read_text())
expected={r['path']:r for r in manifest['files']}
raw=ROOT/'data/raw';raw.mkdir(parents=True,exist_ok=True)
with tarfile.open(archive,'r:gz') as tf:
    members=[x for x in tf if x.isfile()]
    if {x.name for x in members}!=set(expected):raise RuntimeError('Snapshot membership mismatch')
    for member in members:
        if Path(member.name).is_absolute() or '..' in Path(member.name).parts:raise RuntimeError('Unsafe archive path')
        blob=tf.extractfile(member).read();r=expected[member.name]
        if len(blob)!=r['bytes'] or hashlib.sha256(blob).hexdigest()!=r['sha256']:raise RuntimeError('Snapshot checksum failed: '+member.name)
        dst=raw/member.name
        if dst.exists() and hashlib.sha256(dst.read_bytes()).hexdigest()!=r['sha256']:
            raise RuntimeError('Existing reference differs from the saved snapshot: '+str(dst))
        dst.parent.mkdir(parents=True,exist_ok=True)
        if not dst.exists():dst.write_bytes(blob)
print('PASS: restored/verified',len(expected),'saved reference inputs; Open Targets',manifest['open_targets_release'])
