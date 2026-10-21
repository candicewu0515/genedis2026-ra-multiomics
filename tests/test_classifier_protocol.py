import importlib.util,sys
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts';sys.path.insert(0,str(p))
spec=importlib.util.spec_from_file_location('classifier',p/'38_classifier_baselines.py');m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)
rng=np.random.default_rng(4)
X=rng.poisson(30,size=(30,12)).astype(float);y=np.tile([0,1],15)
norm=m.TrainingCounts().fit(X[:20]);fields={k:v.copy() if hasattr(v,'copy') else v for k,v in vars(norm).items()}
one=norm.transform(X[20:21],np.arange(5));other=X[20:24].copy();other[1:]*=100000
assert np.allclose(one,norm.transform(other,np.arange(5))[:1]),'Held-out neighbours changed query normalization'
for k,v in fields.items():assert np.array_equal(v,vars(norm)[k]),'Transform mutated training state'
plan=m.make_plan(X,y,1)
for f in plan:
 assert not set(f.train)&set(f.test)
 for tr,va,n in f.inner:
  assert not (set(tr)|set(va))&set(f.test)
  assert not set(tr)&set(va)
  ref=m.TrainingCounts().fit(X[tr]);assert np.array_equal(ref.log_geo_,n.log_geo_)
A,_=m.prepared(X,np.arange(20),np.arange(20,30),norm,'random',groups=np.arange(3),selected=np.arange(12))
assert A.shape[1]==int(norm.keep_[:3].sum()),'Random comparator does not match retained F2 feature count'
seed=42;yp=np.random.default_rng(seed).permutation(y);groups=[np.arange(4),np.arange(4,8)]
from threadpoolctl import threadpool_limits
with threadpool_limits(limits=1):
 a=m.permutation_stat(X,y,groups,1,seed)
 b=m.evaluate(X,yp,m.make_plan(X,yp,1),'programmes',groups=groups).mean()
assert a==b,'Permutation and observed estimators differ'
print('PASS: train-only references; query independence; nested holdout separation; matched permutation algorithm')
