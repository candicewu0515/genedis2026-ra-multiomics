"""Write the score table (docs/evidence_map.tsv) and check that the final tables, figures and settings agree."""
from pathlib import Path
import importlib.util,sys,json
import numpy as np,pandas as pd
from scipy.stats import mannwhitneyu
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'scripts'))
from config import TAB,PROCESSED
def load(n):return pd.read_csv(TAB/n,sep='\t')
checks=[]
def check(ok,name):
    if not ok:raise AssertionError(name)
    checks.append(name)
s=load('signature_genes.tsv')
S1=s.loc[(s.signature=='primary_FDR')&(s.direction=='up_in_active'),'gene'].tolist()
S2=['CXCL9','CXCL10','CXCL11','STAT1','GBP1','GBP5','IDO1','IRF1','CD274','APOL6']
check(len(S1)==15,'S1 has 15 active-up genes')
check(set(S1)&set(S2)=={'APOL6','GBP5'},'S1 and S2 remain distinct, overlapping at 2 genes')
states=load('state_markers_gse243917.tsv')
rows=[
 dict(score='S1',genes=';'.join(S1),definition='15 FDR genes up in active disease, E-MTAB-8322; 13 measured in GSE243917 RNA/ATAC and GSE198520, 3 in proteomics',role='discovery/transfer/treatment',boundary='multiome inflammation P~0.21; bulk RA/OA and specificity tests negative; good-vs-nonresponder change P~0.131',source='signature_genes.tsv;cross_dataset_signature_stats.tsv;gse89408_tests.tsv;gse198520_tests.tsv'),
 dict(score='S2',genes=';'.join(S2),definition='predefined 10-gene IFN programme',role='Xenium/Visium/bulk stages',boundary='associations do not validate S1; spatial/protein cohorts share blocks',source='xenium_ifng_neighbourhood.tsv;visium_tests.tsv;gse89408_tests.tsv'),
 dict(score='Hallmark_IFNG',genes='HALLMARK_INTERFERON_GAMMA_RESPONSE (MSigDB v2024.1.Hs)',definition='Hallmark genes measurable in each assay',role='GSEA/proteomics/extension sensitivity',boundary='protein random-set specificity supported in lining, not in sublining; not the exact S1 or S2 score',source='proteomics_tests.tsv'),
]
for name,g in states.groupby('programme'):
 check(len(g)==30,name+' has 30 genes')
 rows.append(dict(score=name,genes=';'.join(g.gene),definition='top 30 atlas state markers; author-defined states; ranked on pooled cells',role='GSE296117 and LuCA/TCGA exploratory transfer',boundary='different from S1/S2; no treatment response labels in the fluid cohort; no causal or clinical prediction validation',source='state_markers_gse243917.tsv;ra_fluid_tests.tsv;nsclc_extension_tests.tsv'))
(root/'docs/evidence_map.tsv').parent.mkdir(parents=True,exist_ok=True)
pd.DataFrame(rows).to_csv(root/'docs/evidence_map.tsv',sep='\t',index=False)
pr=pd.read_csv(PROCESSED/'proteomics_region_scores.tsv',sep='\t')
pr=pr[pr.Disease.isin(['NH','RA'])].copy()
pr['region_group']=np.where(pr.region.eq('lining'),'lining',np.where(pr.region.str.endswith('sublining'),'sublining','other'))
dm=pr.groupby(['donor','Disease','region_group'],observed=True).ifng_response.mean().reset_index()
for reg in ['lining','sublining']:
 ra=dm.loc[dm.region_group.eq(reg)&dm.Disease.eq('RA'),'ifng_response']
 nh=dm.loc[dm.region_group.eq(reg)&dm.Disease.eq('NH'),'ifng_response']
 check((len(ra),len(nh))==(6,3),'Protein '+reg+' uses 6 RA and 3 healthy donor means')
 p=load('proteomics_tests.tsv');r=p.loc[p['set'].eq('ifng_response')&p.comparison.eq(f'RA vs NH {reg} (donor means)')].iloc[0]
 check(np.isclose(mannwhitneyu(ra,nh).pvalue,r.p_mannwhitney,atol=5e-6),'Protein '+reg+' plotted donor groups reproduce source P')
display=pd.read_csv(root/'results/figures/main/Fig3g_protein_donor_means.tsv',sep='\t')
pd.testing.assert_frame_equal(dm[dm.region_group.ne('other')].reset_index(drop=True),display,check_dtype=False)
checks.append('Fig3g exported plotted values equal independently aggregated donor means')
tf=load('tf_M6_drivers.tsv')
from statsmodels.stats.multitest import multipletests
q=multipletests(tf.rna_p_perm,method='fdr_bh')[1]
check(np.allclose(q,tf.rna_q_BH,atol=1e-7),'TF-wide BH q agrees with source P values within stored rounding')
check((~tf.both_layers | tf.rna_q_BH.lt(.05)).all(),'Corrected TF candidate flags require BH q <0.05')
protocol=json.loads((TAB/'classifier_protocol.json').read_text())
check(protocol['outer_folds']==5 and protocol['outer_repeats']==5 and protocol['inner_folds']==3,'Classifier has 5×5 outer/3 inner CV')
check(protocol['permutations']==199 and protocol['random_sets']==100,'Classifier has 199 matched-procedure permutations/100 random controls')
null=load('classifier_permutations.tsv');random=load('classifier_random_sets.tsv');folds=load('classifier_folds.tsv')
check(null.groupby('task').size().eq(199).all() and len(null)==398,'Both label-null distributions contain 199 results')
check(random.groupby('task').size().eq(100).all() and len(random)==200,'Both random-control distributions contain 100 results')
check(folds.groupby(['task','model']).size().eq(25).all() and len(folds)==200,'Every observed model has 25 held-out AUCs')
check(np.isfinite(null.mean_auc).all() and null.mean_auc.between(0,1).all(),'Permutation AUCs are finite and in range')
for name in ['multiome_plot_cells.tsv.gz','xenium_plot_CT_21045-A-15.tsv.gz']:
 check((PROCESSED/name).is_file(),'Packaged plot input exists: '+name)
# Completion and semantic checks cover the corrected graph and spatial summaries.
sl=load('sldsc_results.tsv')
check(len(sl)==195 and sl.groupby('gwas').size().eq(65).all(),'All 195 expected S-LDSC regressions are represented')
check(np.isfinite(sl[['tau_z','p_tau','fdr_tau']].to_numpy()).all(),'S-LDSC statistics are finite')
for _,group in sl.groupby(['gwas','level']):
 check(np.allclose(multipletests(group.p_tau,method='fdr_bh')[1],group.fdr_tau,rtol=2e-3,atol=1e-8),'S-LDSC BH agrees within the complete GWAS × annotation-level family (stored rounding tolerance)')
xt=load('xenium_tests.tsv');section=load('xenium_ifng_neighbourhood.tsv')
ra=section[section.condition.eq('RA')]
absent=section[section.pct_lymph_IFNG_pos.eq(0)]
check(len(absent)==1 and absent.coef_IFNGpos_neighbour.isna().all() and absent.p_perm_one_sided.isna().all(),'Zero IFNG-positive exposure is marked non-estimable rather than a null effect')
check(int(xt.iloc[0].n_positive)==int(ra.coef_IFNGpos_neighbour.gt(0).sum()),'Xenium positive count uses the actual coefficient sign')
check('stouffer_z' not in xt.columns,'Spatial summary removes extrapolated Stouffer statistics')
tp=load('tph_treg_tests.tsv')
check(not tp.test.str.contains('Stouffer',case=False).any(),'Tph summary removes extrapolated Stouffer tail tests')
edges=pd.read_csv(root/'results/knowledge_graph/edges.tsv',sep='\t')
nodes=pd.read_csv(root/'results/knowledge_graph/nodes.tsv',sep='\t')
check('score:S2_IFN10' in nodes.id.values,'Knowledge graph contains a separate S2 score node')
s2edges=edges[edges.relation.str.startswith('S2_IFN')]
check(len(s2edges)==2 and s2edges.source.eq('score:S2_IFN10').all(),'Bulk and Visium S2 evidence points to the S2 node')
check(not edges.statistic.fillna('').str.contains('Stouffer',case=False).any(),'Knowledge graph has no extrapolated spatial z evidence')

# Figure source hashes use data paths relative to the project root.
import hashlib
hashes=json.loads((root/'results/figures/main/source_hashes.json').read_text())
for name,want in hashes.items():
 check(not Path(name).is_absolute() and hashlib.sha256((root/name).read_bytes()).hexdigest()==want,'Figure source hash: '+name)
res={'status':'PASS','checks':checks,'note':'Checks that tables, figures and protocol settings agree. It does not test biological truth or external validity.'}
(root/'results/submission_checks.json').write_text(json.dumps(res,indent=2))
print('PASS:',len(checks),'submission consistency checks')
