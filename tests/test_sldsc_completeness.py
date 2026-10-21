"""Reject incomplete/invalid model families before any multiple-testing adjustment."""
from pathlib import Path
import importlib.util,sys
import numpy as np,pandas as pd
scripts=Path(__file__).resolve().parents[1]/'scripts';sys.path.insert(0,str(scripts))
spec=importlib.util.spec_from_file_location('sldsc',scripts/'33_sldsc.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
tests=[(gw,None,f'ST_{i}') for gw in ['A','B','C'] for i in range(65)]
rows=[dict(gwas=g,annotation=n,prop_snps=.1,prop_h2=.2,enrichment=2.,tau=1e-8,tau_se=1e-9,tau_z=10.,p_tau=.01) for g,_,n in tests]
valid=pd.DataFrame(rows);m.validate_complete_results(valid,tests)
def rejected(d):
 try:m.validate_complete_results(d,tests)
 except RuntimeError:return
 raise AssertionError('Invalid family passed the check')
rejected(valid.iloc[:-1]);rejected(pd.concat([valid.iloc[:-1],valid.iloc[[0]]],ignore_index=True))
for col,value in [('tau_se',np.nan),('tau_z',np.inf),('p_tau',1.1),('tau_se',0)]:
 d=valid.copy();d.loc[0,col]=value;rejected(d)
d=valid.copy();d.loc[0,'annotation']='unknown';rejected(d)
print('PASS: complete finite family passes; missing, duplicate, unknown, nonfinite and invalid-statistic families rejected before BH')
