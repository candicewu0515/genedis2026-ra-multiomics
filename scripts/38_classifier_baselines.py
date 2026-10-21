"""Donor/biopsy-level nested-CV baselines with training-only count preprocessing.

Counts are supplied directly to each fold. Expression filtering, median-of-ratios
reference, gene scaling, variable-gene selection, imputation and classifier tuning
are fitted using training observations only, including inside the inner CV.
Observed F1 and every label permutation use exactly the same repeated nested-CV
algorithm and parameter grid. The external 15-gene signature and prioritised gene
list are fixed inputs, not learned from GSE89408 outcome labels.

Tasks and feature sets follow the earlier exploratory analysis. This revision
corrects preprocessing and permutation inference; it is not external validation.
"""
import argparse
import importlib
import json
import hashlib
import warnings
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import t as tdist
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits
from config import FIG, TAB, SEED, MAC_HOMEOSTATIC

step30 = importlib.import_module('30_gse89408_disease_stage')
PROGRAMMES = {'IFN':step30.IFN, 'homeostatic':MAC_HOMEOSTATIC, 'macrophage':step30.MAC,
 'Tph':['PDCD1','CXCL13','MAF','TOX','ICOS'],
 'Treg':['FOXP3','IL2RA','CTLA4','IKZF2','TNFRSF18'],
 'Tcell':['CD3E','CD3D','CD2','CD247'], 'Bplasma':['MS4A1','CD79A','MZB1','JCHAIN']}
C_VALUES = (0.01, 0.1, 1.0)
F1='F1 eight programme scores'

class TrainingCounts:
    """A fixed training reference; transform never estimates cohort parameters."""
    def fit(self, X):
        X=np.asarray(X,dtype=float)
        self.n_features_in_=X.shape[1]
        lib=X.sum(axis=1)
        if np.any(lib<=0):raise ValueError('Empty count library')
        self.keep_=(X>=lib[:,None]*1e-6).mean(axis=0)>=.20
        self.ref_idx_=np.flatnonzero((X>0).all(axis=0))
        if len(self.ref_idx_)==0:raise ValueError('No training-positive genes for count reference')
        self.log_geo_=np.log(X[:,self.ref_idx_]).mean(axis=0)
        train_logsf=np.median(np.log(X[:,self.ref_idx_])-self.log_geo_,axis=1)
        self.log_sf_center_=train_logsf.mean()
        return self

    def size_factors(self, X):
        X=np.asarray(X,dtype=float)
        ref=X[:,self.ref_idx_]
        with np.errstate(divide='ignore'):
            ratios=np.where(ref>0,np.log(ref)-self.log_geo_,np.nan)
        logsf=np.nanmedian(ratios,axis=1)-self.log_sf_center_
        if not np.isfinite(logsf).all():raise ValueError('No positive reference genes in held-out library')
        return np.exp(logsf)

    def transform(self, X, idx):
        idx=np.asarray(idx,dtype=int)
        L=np.log2(np.asarray(X[:,idx],dtype=float)/self.size_factors(X)[:,None]+1)
        return L

@dataclass
class Fold:
    train: np.ndarray
    test: np.ndarray
    normalizer: object
    inner: list


def make_plan(X,y,repeats):
    outer=list(RepeatedStratifiedKFold(n_splits=5,n_repeats=repeats,random_state=SEED).split(X,y))
    plan=[]
    for tr,te in outer:
        inner=[]
        for a,b in StratifiedKFold(3,shuffle=True,random_state=SEED).split(X[tr],y[tr]):
            itr,iva=tr[a],tr[b]
            inner.append((itr,iva,TrainingCounts().fit(X[itr])))
        plan.append(Fold(tr,te,TrainingCounts().fit(X[tr]),inner))
    return plan


def prepared(X, tr, te, norm, kind, groups=None, selected=None):
    if kind=='covariates':
        A,B=X[tr].copy(),X[te].copy()
        med=np.nanmedian(A[:,0]);A[np.isnan(A[:,0]),0]=med;B[np.isnan(B[:,0]),0]=med
    elif kind=='programmes':
        idx=np.unique(np.concatenate(groups));idx=idx[norm.keep_[idx]]
        Ltr,Lte=norm.transform(X[tr],idx),norm.transform(X[te],idx)
        mu,sd=Ltr.mean(0),Ltr.std(0)+1e-9
        Ztr,Zte=(Ltr-mu)/sd,(Lte-mu)/sd
        pos={int(g):j for j,g in enumerate(idx)}
        A=np.column_stack([Ztr[:,[pos[int(g)] for g in group if int(g) in pos]].mean(1) if any(int(g) in pos for g in group) else np.zeros(len(tr)) for group in groups])
        B=np.column_stack([Zte[:,[pos[int(g)] for g in group if int(g) in pos]].mean(1) if any(int(g) in pos for g in group) else np.zeros(len(te)) for group in groups])
    else:
        if kind=='variable':
            candidates=np.flatnonzero(norm.keep_)
            train=norm.transform(X[tr],candidates)
            idx=candidates[np.argsort(-train.var(axis=0),kind='stable')[:2000]]
        else:
            idx=np.asarray(selected)[norm.keep_[selected]]
            if kind=='random':idx=idx[:int(norm.keep_[groups].sum())]
        if not len(idx):
            A,B=np.zeros((len(tr),1)),np.zeros((len(te),1))
        else:A,B=norm.transform(X[tr],idx),norm.transform(X[te],idx)
    sc=StandardScaler().fit(A)
    return sc.transform(A),sc.transform(B)


def parameters(kind):
    return [(c,l1) for c in C_VALUES for l1 in ((.1,.5,.9) if kind=='variable' else (None,))]


def estimator(param):
    c,l1=param
    if l1 is None:return LogisticRegression(C=c,max_iter=5000,class_weight='balanced',random_state=SEED)
    return LogisticRegression(C=c,solver='saga',l1_ratio=l1,
        max_iter=10000,tol=1e-3,class_weight='balanced',random_state=SEED)


def fit_checked(model,A,y):
    for _ in range(3):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always',ConvergenceWarning)
            model.fit(A,y)
        if not any(issubclass(w.category,ConvergenceWarning) for w in caught):return model
        model.max_iter *= 4
    raise RuntimeError('Classifier did not converge after explicit retries')


def evaluate(X,y,plan,kind,groups=None,selected=None,return_details=False):
    params=parameters(kind);aucs=[];chosen=[]
    for fold in plan:
        inner_scores=np.zeros((len(params),len(fold.inner)))
        for k,(tr,va,norm) in enumerate(fold.inner):
            A,B=prepared(X,tr,va,norm,kind,groups,selected)
            for j,param in enumerate(params):
                model=fit_checked(estimator(param),A,y[tr])
                inner_scores[j,k]=roc_auc_score(y[va],model.predict_proba(B)[:,1])
        param=params[int(np.argmax(inner_scores.mean(1)))]
        A,B=prepared(X,fold.train,fold.test,fold.normalizer,kind,groups,selected)
        model=fit_checked(estimator(param),A,y[fold.train])
        aucs.append(roc_auc_score(y[fold.test],model.predict_proba(B)[:,1]));chosen.append(param)
    return (np.asarray(aucs),chosen) if return_details else np.asarray(aucs)


def permutation_stat(X,y,groups,repeats,seed):
    yp=np.random.default_rng(seed).permutation(y)
    # Same algorithm as observed: label-stratified 5 x repeats outer, 3-fold inner tuning.
    plan=make_plan(X,yp,repeats)
    return evaluate(X,yp,plan,'programmes',groups=groups).mean()


def corrected_t(diff,n_train,n_test):
    k=len(diff);t=diff.mean()/np.sqrt((1/k+n_test/n_train)*diff.var(ddof=1)+1e-12)
    return t,2*tdist.sf(abs(t),k-1)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--jobs',type=int,default=6)
    parser.add_argument('--permutations',type=int,default=199)
    parser.add_argument('--random-sets',type=int,default=100)
    parser.add_argument('--outer-repeats',type=int,default=5)
    args=parser.parse_args()
    # A clean project has no cache. Concurrent validation branches may reuse only
    # byte-verified outputs of this exact code, input files and resampling protocol.
    sources=[Path(__file__),Path(__file__).with_name('config.py'),Path(step30.__file__),
        step30.D/'GSE89408_GEO_count_matrix_rename.txt.gz',step30.D/'GSE89408_series_matrix.txt.gz',
        TAB/'signature_genes.tsv',TAB/'gene_prioritisation.tsv']
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    fingerprint={str(p.relative_to(Path(__file__).resolve().parents[1])):sha(p) for p in sources}
    fingerprint['parameters']={'permutations':args.permutations,'random_sets':args.random_sets,'outer_repeats':args.outer_repeats}
    cache=TAB/'classifier_cache.json'
    if cache.exists():
        old=json.loads(cache.read_text())
        if old['inputs']==fingerprint and all((TAB/n).is_file() and sha(TAB/n)==v for n,v in old['outputs'].items()):
            print('Reusing byte-verified classifier outputs from this exact code, inputs and protocol',flush=True)
            return
    rng=np.random.default_rng(SEED)
    meta=step30.metadata()
    raw=pd.read_csv(step30.D/'GSE89408_GEO_count_matrix_rename.txt.gz',sep='\t',index_col=0)[meta.index].round()
    raw=raw[~raw.index.duplicated()]
    genes=raw.index;Xall=raw.T.to_numpy(dtype=float)
    meta['log_lib']=np.log(Xall.sum(axis=1))
    sig=pd.read_csv(TAB/'signature_genes.tsv',sep='\t')
    progs=dict(PROGRAMMES,signature=sig.loc[(sig.signature=='primary_FDR')&(sig.direction=='up_in_active'),'gene'].tolist())
    groups=[np.flatnonzero(genes.isin(g)) for g in progs.values()]
    if any(not len(g) for g in groups):raise ValueError('Empty programme')
    priority=pd.read_csv(TAB/'gene_prioritisation.tsv',sep='\t',index_col=0)
    top30=[g for g in priority.index if g in genes][:30];idx30=np.flatnonzero(genes.isin(top30))
    tasks={'T1 RA vs OA':meta.group.isin(['early RA','established RA','OA']),
      'T2 ACPA+ vs ACPA- (RA)':meta.group.isin(['early RA','established RA'])&meta['acpa positivity'].isin(['Positive','Negative'])}
    rows,folds,coefs,perms,randoms=[],[],[],[],[]
    for task,sel in tasks.items():
        m=meta.loc[sel];X=Xall[sel];y=(m.group.ne('OA') if task.startswith('T1') else m['acpa positivity'].eq('Positive')).to_numpy(dtype=int)
        Xc=m[['age','female','log_lib']].to_numpy(dtype=float)
        print(f'{task}: {len(y)} observations; training-only preprocessing',flush=True)
        plan=make_plan(X,y,args.outer_repeats)
        res={}
        for name,kind,A,kw in [('F0 covariates','covariates',Xc,{}),(F1,'programmes',X,{'groups':groups}),
           ('F2 top-30 prioritised genes','fixed',X,{'selected':idx30}),('F3 genome-wide elastic net','variable',X,{})]:
            res[name]=evaluate(A,y,plan,kind,**kw)
            print(f'  {name}: mean AUC {res[name].mean():.4f}',flush=True)
        # Parallel workers share immutable count matrices and fold references, with BLAS limited to one thread.
        sets=[rng.permutation(len(genes)) for _ in range(args.random_sets)]
        rand=np.asarray(Parallel(n_jobs=args.jobs,prefer='threads',verbose=5)(delayed(lambda idx:evaluate(X,y,plan,'random',groups=idx30,selected=idx).mean())(idx) for idx in sets))
        perm_seeds=rng.integers(0,2**31-1,size=args.permutations)
        null=np.asarray(Parallel(n_jobs=args.jobs,prefer='threads',verbose=5)(delayed(permutation_stat)(X,y,groups,args.outer_repeats,int(seed)) for seed in perm_seeds))
        for j,v in enumerate(rand):randoms.append({'task':task,'draw':j,'mean_auc':v})
        for j,(seed,v) in enumerate(zip(perm_seeds,null)):perms.append({'task':task,'draw':j,'seed':int(seed),'mean_auc':v})
        for name,a in res.items():
            rows.append(dict(task=task,model=name,n=len(y),n_pos=int(y.sum()),auc_mean=a.mean(),auc_sd=a.std(ddof=1)))
            folds.extend(dict(task=task,model=name,fold=i,auc=v) for i,v in enumerate(a))
        p_rand=(1+int((rand>=res['F2 top-30 prioritised genes'].mean()).sum()))/(1+len(rand)) if len(rand) else np.nan
        p_perm=(1+int((null>=res[F1].mean()).sum()))/(1+len(null)) if len(null) else np.nan
        rows.append(dict(task=task,model=f'F4 random 30-gene sets (median of {len(rand)})',n=len(y),n_pos=int(y.sum()),auc_mean=float(np.median(rand)) if len(rand) else np.nan,auc_sd=float(rand.std(ddof=1)) if len(rand)>1 else np.nan,note=f'F2 vs random empirical p {p_rand:.4g}; SD across random-set mean AUCs'))
        rows.append(dict(task=task,model=f'F1 label-permutation null (median of {len(null)})',n=len(y),n_pos=int(y.sum()),auc_mean=float(np.median(null)) if len(null) else np.nan,auc_sd=float(null.std(ddof=1)) if len(null)>1 else np.nan,note=f'F1 permutation p {p_perm:.4g}; identical repeated nested-CV procedure'))
        diff=res[F1]-res['F3 genome-wide elastic net'];t,p=corrected_t(diff,len(y)*.8,len(y)*.2)
        rows.append(dict(task=task,model='F1 minus F3 (corrected resampled t-test)',n=len(y),n_pos=int(y.sum()),auc_mean=diff.mean(),auc_sd=np.nan,note=f't {t:.4g}, p {p:.4g}; approximate corrected resampled test'))
        norm=TrainingCounts().fit(X);idx=np.arange(len(y));A,_=prepared(X,idx,idx,norm,'programmes',groups=groups)
        fit=fit_checked(estimator((1,None)),A,y)
        coefs.extend(dict(task=task,programme=k,coef=float(v)) for k,v in zip(progs,fit.coef_[0]))
        print(f'{task} complete; matched-procedure permutation p {p_perm:.4g}',flush=True)
    out=pd.DataFrame(rows)
    out.to_csv(TAB/'classifier_baselines.tsv',sep='\t',index=False,float_format='%.6g')
    pd.DataFrame(folds).to_csv(TAB/'classifier_folds.tsv',sep='\t',index=False,float_format='%.6g')
    pd.DataFrame(coefs).to_csv(TAB/'classifier_F1_coefficients.tsv',sep='\t',index=False,float_format='%.6g')
    pd.DataFrame(perms,columns=['task','draw','seed','mean_auc']).to_csv(TAB/'classifier_permutations.tsv',sep='\t',index=False,float_format='%.9g')
    pd.DataFrame(randoms,columns=['task','draw','mean_auc']).to_csv(TAB/'classifier_random_sets.tsv',sep='\t',index=False,float_format='%.9g')
    (TAB/'classifier_protocol.json').write_text(json.dumps({'revision':'2026-10-20 training-only preprocessing and matched permutation','outer_folds':5,'outer_repeats':args.outer_repeats,'inner_folds':3,'C_values':list(C_VALUES),'elastic_net_l1_ratios':[.1,.5,.9],'permutations':args.permutations,'random_sets':args.random_sets,'random_control':'fixed random ranking per draw; training-filtered count matched to the retained genes of F2 in each inner/outer fold','seed':SEED,'filter':'CPM >=1 in >=20% of training observations','normalization':'training-positive geometric reference and training size-factor centering','all_fitted_in_training':['expression filter','normalization reference','gene means/SDs','feature scaling','variable genes','age imputation','hyperparameter tuning'],'signature_genes':progs['signature'],'priority_genes':top30,'clinical_validation':'none; internal nested CV only'},indent=2))
    print(out.to_string(index=False),flush=True)
    fig,axes=plt.subplots(1,2,figsize=(12,4))
    for ax,task in zip(axes,tasks):
        d=out[(out.task==task)&out.model.str.startswith(('F0','F1 eight','F2','F3'))]
        ax.errorbar(range(len(d)),d.auc_mean,yerr=d.auc_sd,fmt='o');ax.set_xticks(range(len(d)),d.model.str.split().str[0]);ax.axhline(.5,color='.5',ls='--');ax.set_title(task);ax.set_ylabel('Nested CV AUC ± fold SD')
    fig.tight_layout();fig.savefig(FIG/'classifier_baselines.png',dpi=150);plt.close(fig)
    names=['classifier_baselines.tsv','classifier_folds.tsv','classifier_F1_coefficients.tsv',
        'classifier_permutations.tsv','classifier_random_sets.tsv','classifier_protocol.json']
    cache.write_text(json.dumps({'inputs':fingerprint,'outputs':{n:sha(TAB/n) for n in names}},indent=2))


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
