"""Build the submission archive; run from the project root after validation."""
from pathlib import Path
import shutil,hashlib,zipfile,re
ROOT=Path(__file__).resolve().parents[1]
NAME='genedis2026_RA_XiaWu_revised_20261020'
DIST=ROOT/'dist';DIST.mkdir(exist_ok=True);stage=DIST/NAME
if stage.exists():shutil.rmtree(stage)
stage.mkdir()
for name in ['README.md','environment.yml','environment_ldsc.yml','.gitignore','SUBMISSION.json']:
 if (ROOT/name).exists():shutil.copy2(ROOT/name,stage/name)
for name in ['scripts','docs','env','tests','presentation','references','data/processed','results/tables','results/figures/main','results/knowledge_graph']:
 src=ROOT/name
 if not src.exists():raise FileNotFoundError(src)
 shutil.copytree(src,stage/name,ignore=shutil.ignore_patterns('__pycache__','.DS_Store','*.pyc'))
(stage/'data/raw').mkdir(parents=True)
shutil.copy2(ROOT/'data/raw/SHA256SUMS.txt',stage/'data/raw/SHA256SUMS.txt')
for p in (ROOT/'results').glob('*'):
 if p.is_file() and p.suffix in ['.tsv','.json','.md']:shutil.copy2(p,stage/'results'/p.name)
(stage/'logs').mkdir()
# Local paths in the logs are replaced, so the archived logs do not depend on this computer.
for p in (ROOT/'logs').iterdir():
 if p.is_file() and p.suffix in ['.log','.tsv','.json']:
  text=p.read_text(errors='replace').replace(str(ROOT),'<project>')
  # Any other absolute project path (for example a separate validation copy) and the home folder are masked too.
  text=re.sub(r'/(?:Users|home|private)/[^\n\'"]*?/(?:project|genedis_ra)(?=[/\s]|$)','<project>',text)
  text=re.sub(r'/Users/[^/\s]+','<home>',text)
  (stage/'logs'/p.name).write_text(text)
files=sorted(p for p in stage.rglob('*') if p.is_file())
with (stage/'MANIFEST.tsv').open('w') as f:
 for p in files:
  h=hashlib.sha256()
  with p.open('rb') as stream:
   for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
  f.write(f'{p.relative_to(stage)}\t{p.stat().st_size}\t{h.hexdigest()}\n')
archive=DIST/(NAME+'.zip')
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(stage.rglob('*')):
  if p.is_file():z.write(p,Path(NAME)/p.relative_to(stage))
h=hashlib.sha256()
with archive.open('rb') as f:
 for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
(DIST/(NAME+'.sha256')).write_text(h.hexdigest()+'  '+archive.name+'\n')
print(archive,archive.stat().st_size,'bytes;',len(files),'manifest files')
