"""Build value-verified native LDSC Parquet annotation caches from public gzip inputs.

Uses the same pandas parser as the pinned LDSC tool. The original files remain
unchanged. Chromosome 1 is installed last because LDSC uses its suffix for all 22.
"""
from pathlib import Path
import hashlib,json,sys,tempfile,os,shutil,time
import pandas as pd
import pyarrow as pa
pa.set_cpu_count(1)
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools/polyfun'))
from ldsc_polyfun.parse import read_csv
ref=ROOT/'data/raw/ldsc/GRCh38/baselineLD_v2.2'
record=ROOT/'data/interim/ldsc/annotation_parquet_cache.json'
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def main():
 started=time.time();sources=[ref/f'baselineLD.{i}.annot.gz' for i in range(1,23)]
 raw={p.name:digest(p) for p in sources}
 if record.exists():
  old=json.loads(record.read_text())
  if old['raw_sha256']==raw and all((ref/n).exists() and digest(ref/n)==h for n,h in old['parquet_sha256'].items()):
   (ROOT/'logs').mkdir(exist_ok=True)
   shutil.copy2(record,ROOT/'logs/annotation_parquet_cache.json')
   print('PASS: reused verified native annotation storage cache',flush=True);return
 record.parent.mkdir(parents=True,exist_ok=True)
 staging=Path(tempfile.mkdtemp(prefix='annotation_parquet_',dir=ROOT/'data/interim/ldsc'))
 rows=[]
 try:
  for src in sources:
   df=read_csv(str(src),header=0,compression='gzip')
   out=staging/(src.name.removesuffix('.gz')+'.parquet')
   df.to_parquet(out,index=True,compression='snappy')
   roundtrip=pd.read_parquet(out)
   pd.testing.assert_frame_equal(df,roundtrip,check_exact=True)
   rows.append({'source':src.name,'cache':out.name,'rows':len(df),'columns':len(df.columns),'exact_roundtrip':True,'sha256':digest(out)})
   print('PASS exact values:',src.name,len(df),'rows',flush=True)
   del df,roundtrip
  # Existing jobs have already chosen their gzip suffix. Future jobs switch only
  # after all 22 cache files have been installed and chromosome 1 appears.
  for row in rows[1:]+rows[:1]:os.replace(staging/row['cache'],ref/row['cache'])
  record.parent.mkdir(parents=True,exist_ok=True)
  record.write_text(json.dumps({'status':'PASS','seconds':time.time()-started,'parser':'pinned LDSC pandas whitespace parser; original gzip unchanged','raw_sha256':raw,'parquet_sha256':{r['cache']:r['sha256'] for r in rows},'files':rows},indent=2))
  (ROOT/'logs').mkdir(exist_ok=True)
  shutil.copy2(record,ROOT/'logs/annotation_parquet_cache.json')
 finally:shutil.rmtree(staging,ignore_errors=True)
 print('PASS: 22 annotation files round-trip exactly; native cache installed',flush=True)
if __name__=='__main__':main()
