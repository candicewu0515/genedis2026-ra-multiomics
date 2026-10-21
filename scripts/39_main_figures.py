"""Draw the project figures in a compact Nature-style layout.
All data panels use outputs that the project has already produced.

Run from any directory in the genedis_ra environment. By default the project root is
the folder above scripts/, as this file sits at scripts/39_main_figures.py.
This script only draws figures. It does not change any analysis or statistical threshold.
"""
from pathlib import Path
import argparse
import json
import shutil
import hashlib
import textwrap
import numpy as np
import pandas as pd
import h5py
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.stats import spearmanr, norm

parser = argparse.ArgumentParser()
parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
parser.add_argument('--out', type=Path)
args = parser.parse_args()
ROOT = args.root
OUT = args.out or ROOT / 'results/figures/main'
OUT.mkdir(parents=True, exist_ok=True)
TAB = ROOT / 'results/tables'
PRO = ROOT / 'data/processed'
INT = ROOT / 'data/interim'
MM = 1/25.4
BLUE, RED, GREEN, GOLD = '#0072B2', '#D55E00', '#009E73', '#E69F00'
PURPLE, CYAN, PINK, GREY = '#7B3294', '#56B4E9', '#CC79A7', '#777777'
CT = {'Tcell':BLUE, 'NK':CYAN, 'Bcell':PURPLE, 'plasma':PINK,
      'myeloid':RED, 'fibroblast':'#555555', 'endothelial':GREEN, 'mural':GOLD}
CN = {'Tcell':'T', 'NK':'NK', 'Bcell':'B', 'plasma':'Plasma', 'myeloid':'Myeloid',
      'fibroblast':'Fibroblast', 'endothelial':'Endothelial', 'mural':'Mural'}
CC = {'healthy':BLUE, 'NH':BLUE, 'OA':GREY, 'RA':RED}
M6 = 'M-6: STAT1+ CXCL10+'
GW = ['RA_Ishigaki_EUR','RA_FinnGen_R12','T2D_FinnGen_R12_control']
GS = ['EUR','FinnGen','T2D']
SEQ = LinearSegmentedColormap.from_list('evidence', ['#FFFFFF', BLUE])
DIV = LinearSegmentedColormap.from_list('signed', [BLUE,'#FFFFFF',RED])
plt.rcParams.update({'font.family':['Arial','DejaVu Sans'],'font.size':5.5,'axes.labelsize':5.5,
 'axes.titlesize':6.5,'axes.titleweight':'bold','xtick.labelsize':5,'ytick.labelsize':5,
 'legend.fontsize':5,'legend.frameon':False,'axes.linewidth':.5,
 'axes.spines.top':False,'axes.spines.right':False,'xtick.major.size':2,
 'ytick.major.size':2,'xtick.major.width':.5,'ytick.major.width':.5,
 'axes.labelpad':1.5,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
 'savefig.dpi':600,'axes.unicode_minus':True,'mathtext.default':'regular'})
SOURCES = {}
PANEL_GUIDE = {}

def read(name, index=None, where=TAB):
    path = where / (name if name.endswith(('.tsv','.tsv.gz','.csv','.gz')) else name+'.tsv')
    if not path.is_file():
        raise FileNotFoundError(f'Missing figure input: {path}. Run the analysis step that makes it, or use the complete submission package.')
    SOURCES[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path, sep='\t', index_col=index)

def page(height=169):
    f=plt.figure(figsize=(183*MM,height*MM),facecolor='white')
    f._height_mm=height
    f._panels=[]
    return f

def cell(f, letter, title, x, y, w, h, left=10, right=2, bottom=9, top=9):
    """Add one panel. Position and size are in mm from the top-left corner, with space kept for the label."""
    H=f._height_mm
    f.text(x/183,(H-y-1)/H,letter,fontsize=8,fontweight='bold',va='top')
    f.text((x+4)/183,(H-y-1.5)/H,textwrap.fill(title,width=max(20,int((w-4)*.69))),fontsize=6.5,fontweight='bold',va='top',linespacing=1.12)
    ax=f.add_axes([(x+left)/183,(H-y-h+bottom)/H,(w-left-right)/183,(h-top-bottom)/H])
    f._panels.append((letter,title,ax))
    return ax

def standard(f, r, c, letter, title, span=1, left=10, **kw):
    return cell(f,letter,title,2+c*60,2+r*55,59+(span-1)*60,54,left=left,**kw)

def note(ax, text, x=.02,y=.98,ha='left',color='#222222'):
    ax.text(x,y,text,transform=ax.transAxes,va='top',ha=ha,fontsize=5,
            color=color,linespacing=1.15,bbox={'fc':'white','ec':'none','pad':.5,'alpha':.88})

def ptext(p):
    if p==0: return 'P < 0.001'
    if p<.001: return f'P = {p:.1e}'
    return f'P = {p:.3f}'.rstrip('0').rstrip('.')

def clean(ax):
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)

def zero(ax, x=True, y=True):
    if x: ax.axvline(0,color='.8',lw=.5,zorder=0)
    if y: ax.axhline(0,color='.8',lw=.5,zorder=0)

def matrix(ax, data, cmap=DIV, vmin=-2, vmax=2, annotate=False, fmt='.1f', labels=None):
    a=np.ma.masked_invalid(np.asarray(data,dtype=float))
    cm=plt.get_cmap(cmap).copy(); cm.set_bad('#EEEEEE')
    im=ax.imshow(a,aspect='auto',cmap=cm,vmin=vmin,vmax=vmax,interpolation='none')
    clean(ax)
    if isinstance(data,pd.DataFrame):
        ax.set_xticks(range(data.shape[1]),data.columns)
        ax.set_yticks(range(data.shape[0]),data.index if labels is None else labels)
    ax.set_xticks(np.arange(-.5,a.shape[1],1),minor=True)
    ax.set_yticks(np.arange(-.5,a.shape[0],1),minor=True)
    ax.grid(which='minor',color='white',lw=.5); ax.tick_params(which='minor',length=0)
    if annotate:
        for i in range(a.shape[0]):
            for j in range(a.shape[1]):
                if not np.ma.is_masked(a[i,j]):
                    frac=(a[i,j]-vmin)/(vmax-vmin)
                    dark=abs(frac-.5)>.34 if vmin<0 else frac>.68
                    ax.text(j,i,format(a[i,j],fmt),ha='center',va='center',fontsize=5,color='white' if dark else '#111111')
                else: ax.text(j,i,'—',ha='center',va='center',fontsize=5,color='.4')
    return im

def bar_key(ax, im, label):
    cb=ax.figure.colorbar(im,ax=ax,fraction=.045,pad=.035,aspect=18)
    cb.ax.set_title(label,fontsize=5,pad=3); cb.ax.tick_params(labelsize=5,length=1.5,width=.4)
    cb.outline.set_linewidth(.4)

def dots(ax, groups, labels, colors, ylabel, counts=True, sizes=7, showbox=False):
    rng=np.random.default_rng(0)
    for i,(v,col) in enumerate(zip(groups,colors)):
        v=np.asarray(v,dtype=float); v=v[np.isfinite(v)]
        if showbox and len(v)>3:
            b=ax.boxplot([v],positions=[i],widths=.55,showfliers=False,patch_artist=True,
                boxprops={'facecolor':'white','edgecolor':col,'lw':.6},
                medianprops={'color':'black','lw':.8},whiskerprops={'color':col,'lw':.5},
                capprops={'color':col,'lw':.5})
        ax.scatter(i+rng.uniform(-.17,.17,len(v)),v,s=sizes,c=col,lw=0,alpha=.85,rasterized=len(v)>100)
        if len(v): ax.plot([i-.22,i+.22],[np.median(v)]*2,color='black',lw=.9,zorder=4)
    ax.set_xticks(range(len(labels)),[f'{l}\nn={len(v)}' if counts else l for l,v in zip(labels,groups)])
    ax.set_xlim(-.5,len(labels)-.5); ax.set_ylabel(ylabel)


def labels_right(ax, df, xcol, ycol, names, ylim=None):
    pts=df.reindex(names).dropna(subset=[xcol,ycol])
    if pts.empty:return
    lo,hi=ax.get_ylim() if ylim is None else ylim
    ys=np.linspace(max(lo+.1*(hi-lo),pts[ycol].min()),min(hi-.06*(hi-lo),pts[ycol].max()),len(pts))
    pts=pts.sort_values(ycol)
    xx=ax.get_xlim()[1]-.03*np.ptp(ax.get_xlim())
    for (name,row),yy in zip(pts.iterrows(),ys):
        ax.annotate(name,xy=(row[xcol],row[ycol]),xytext=(xx,yy),ha='right',va='center',
          fontsize=5,fontstyle='italic',arrowprops={'arrowstyle':'-','lw':.4,'color':'.45'},
          bbox={'fc':'white','ec':'none','pad':.2})

def save(f,name,guide):
    PANEL_GUIDE[name]=guide
    f.canvas.draw()
    for ext in ['pdf','svg','png']:
        f.savefig(OUT/f'{name}.{ext}',dpi=600)
    plt.close(f)
    print(name+' done',flush=True)

def atlas():
    path=PRO/'multiome_plot_cells.tsv.gz'
    if not path.exists():
        # Older workspaces may not have this table. Step 07 writes it on a clean run.
        with h5py.File(INT/'gse243917_multiome.h5mu','r') as h:
            um=h['obsm/X_umap'][:];obs=ad.io.read_elem(h['obs']);ro=ad.io.read_elem(h['mod/rna/obs']).reindex(obs.index)
        plot=ro[['author_ct','author_subtype','donor','disease']].copy()
        plot['umap1']=um[:,0];plot['umap2']=um[:,1]
        plot.to_csv(path,sep='\t',index=True,compression={'method':'gzip','mtime':0})
    plot=read(path.name,index=0,where=PRO)
    return plot[['umap1','umap2']].to_numpy(),plot,plot


def spatial_plot_cells():
    path=PRO/'xenium_plot_CT_21045-A-15.tsv.gz'
    if not path.exists():
        # Rebuild from the documented public inputs, never from a hidden plotting cache.
        import importlib
        step=importlib.import_module('20_spatial_xenium')
        section=next(x for x in step.sections() if x['block']=='CT_21045-A-15' and x['panel']=='5K')
        a=step.annotate(step.load(section),'5K')
        positive=np.asarray(a[:,'CXCL13'].layers['counts'].todense()).ravel()>0
        plot=pd.DataFrame({'x':a.obs.x.to_numpy(),'y':a.obs.y.to_numpy(),'cell_type':a.obs.cell_type.to_numpy(),'CXCL13':positive})
        plot.to_csv(path,sep='\t',index=False,compression={'method':'gzip','mtime':0})
    return read(path.name,where=PRO)


def protein_groups(pr):
    d=pr[pr.Disease.isin(['NH','RA'])].copy()
    d['region_group']=np.where(d.region.eq('lining'),'lining',np.where(d.region.str.endswith('sublining'),'sublining','other'))
    d=d[d.region_group.isin(['lining','sublining'])]
    return d.groupby(['donor','Disease','region_group'],observed=True).ifng_response.mean().reset_index()


def fig1():
    f=page(158)
    # The data inventory is a small table, so most of the page is left for the data panels.
    ax=cell(f,'a','Study cohorts and evidence layers',2,1,179,30,left=4,top=7,bottom=0)
    ax.set_axis_off()
    headers=['Genetics','Joint RNA + ATAC','Macrophages','Bulk synovium','Spatial tissue','Protein / context']
    rows=[['3 RA GWAS','31,520 nuclei','26 donors','218 biopsies','17 Xenium sections','24 LCM regions'],
          ['EUR · EAS · FinnGen','12 donors','Healthy · active · remission','46 paired anti-TNF patients','12 Visium sections','Blood · lung cancer'],
          ['S-LDSC · PIP · QTL · MR','GSE243917','E-MTAB-8322','GSE89408 · GSE198520','Matched tissue blocks','Targets · perturbations']]
    for j,head in enumerate(headers):
        x=(j+.02)/6
        ax.text(x,.96,head,fontweight='bold',fontsize=6,va='top')
        ax.plot([x,x+.153],[.76,.76],color=[PURPLE,BLUE,RED,GOLD,GREEN,GREY][j],lw=1.3)
        for k,row in enumerate(rows):ax.text(x,.66-k*.24,row[j],fontsize=5,va='top')
    ax.set_xlim(0,1);ax.set_ylim(0,1)
    um,obs,ro=atlas(); ct=obs.author_ct.astype(str); rng=np.random.default_rng(0); ix=rng.permutation(len(ct))
    ax=cell(f,'b','Joint multimodal atlas',2,33,66,65,left=5,bottom=6)
    ax.scatter(um[ix,0],um[ix,1],c=[CT.get(v,GREY) for v in ct.iloc[ix]],s=.18,lw=0,rasterized=True)
    for k in CT:
        m=ct.values==k
        if m.any():
            x,y=np.median(um[m],axis=0)
            ax.text(x,y,f'{CN[k]}\n{m.sum():,}',ha='center',va='center',fontsize=5,
                bbox={'fc':'white','ec':'none','alpha':.85,'pad':.6})
    m=ro.author_subtype.astype(str).values==M6
    ax.scatter(um[m,0],um[m,1],s=2,c='black',lw=0)
    note(ax,f'M-6: {m.sum()} nuclei',x=.98,y=.02,ha='right')
    ax.set_xticks([]);ax.set_yticks([]);clean(ax)
    # Draw two small arrows as axes instead of a full box around the embedding.
    ax.annotate('',xy=(.17,.04),xytext=(.04,.04),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.5})
    ax.annotate('',xy=(.04,.18),xytext=(.04,.04),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.5})
    ax.text(.17,.005,'UMAP 1',transform=ax.transAxes,fontsize=5,ha='right')
    ax.text(.005,.18,'2',transform=ax.transAxes,fontsize=5)
    ax=cell(f,'c','Donor composition and yield',70,33,57,65,left=9,bottom=14)
    counts=pd.crosstab(obs.donor,obs.author_ct).reindex(columns=list(CT),fill_value=0)
    disease=ro.groupby('donor',observed=True).disease.first().astype(str)
    order=sorted(counts.index,key=lambda d:(disease[d],d))
    counts=counts.loc[order]; frac=counts.div(counts.sum(axis=1),axis=0)*100
    bottom=np.zeros(len(frac))
    for k in CT:
        ax.bar(range(len(frac)),frac[k],bottom=bottom,color=CT[k],width=.8,lw=0,label=CN[k]);bottom+=frac[k]
    ax.set_ylim(0,110);ax.set_ylabel('Nuclei (%)')
    ax.set_xticks(range(len(frac)),[d.replace('BRI-','') for d in frac.index],rotation=90)
    for i,d in enumerate(order):
        ax.text(i,102,f'{counts.loc[d].sum()/1000:.1f}',ha='center',fontsize=5)
    ax.text(.01,1.03,'Total nuclei (×10³)',transform=ax.transAxes,fontsize=5)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.25),ncol=4,handlelength=.9,columnspacing=.7,handletextpad=.3)
    # Disease squares under the donor IDs. The labels come from the GEO metadata.
    for i,d in enumerate(order):
        ax.plot(i,-.25,'s',color=CC[disease[d]],ms=2.5,transform=ax.get_xaxis_transform(),clip_on=False)
    ax=cell(f,'d','Modality contribution by state',130,33,51,65,left=14,bottom=10)
    mw=read('wnn_modality_weights'); order_ct=mw.groupby('cell_type').median_rna_weight.median().sort_values().index
    for i,k in enumerate(order_ct):
        v=mw[mw.cell_type==k].median_rna_weight
        ax.plot([v.min(),v.max()],[i,i],c=CT[k],lw=.7)
        ax.scatter(v,np.full(len(v),i)+rng.uniform(-.15,.15,len(v)),s=5,c=CT[k],lw=0)
        ax.scatter([v.median()],[i],s=13,c='black',marker='|',lw=1)
    ax.axvline(.5,c='.6',ls='--',lw=.5)
    ax.set_yticks(range(len(order_ct)),[CN[k] for k in order_ct]);ax.set_xlabel('RNA weight (ATAC = 1 − RNA)')
    ax.set_xlim(.455,.53);note(ax,f'{len(mw)} subtype summaries',y=1.04)
    ax=cell(f,'e','Integration benchmark: biology, mixing and rare-state recovery',2,103,88,53,left=22,bottom=10)
    ib=read('integration_benchmark',index=0)
    cols=['ARI_author_ct','ARI_author_subtype','ARI_ATAC_leiden','batch_neighbour_entropy','M6_neighbour_purity','M6_best_cluster_F1','overall_scIB_weights']
    short=['Cell type\nARI','Subtype\nARI','ATAC\nARI','Sample\nmixing','M-6\npurity','M-6\nF1','Overall']
    names=['Concat. + Harmony','WNN','ATAC + Harmony','RNA + Harmony','RNA uncorrected','MOFA+ + Harmony']
    dat=ib[cols].copy();dat.columns=short;dat.index=names
    im=matrix(ax,dat,cmap=SEQ,vmin=0,vmax=1,annotate=True,fmt='.2f')
    # The raw metrics all use the same true 0–1 scale. No hidden column-wise normalisation.
    ax=cell(f,'f','Joint latent-factor variance',94,103,42,53,left=9,bottom=10)
    mv=read('mofa_variance_explained',index=0).head(8)
    y=np.arange(len(mv));ax.barh(y-.16,mv.rna,.31,color=RED,label='RNA');ax.barh(y+.16,mv.atac,.31,color=BLUE,label='ATAC')
    ax.set_yticks(y,[i.replace('Factor','F') for i in mv.index]);ax.invert_yaxis();ax.set_xlabel('Variance explained (%)')
    ax.legend(loc='lower right',handlelength=.8)
    ax=cell(f,'g','Donor RNA / ATAC quality (size: nuclei)',139,103,42,53,left=9,bottom=10)
    q=read('gse243917_qc_per_sample',index=0).join(read('atac_qc_per_sample',index=0),rsuffix='_atac')
    for dis in ['RA','OA']:
        v=q[q.disease==dis];ax.scatter(v.median_genes/1000,v.median_peaks/1000,s=v.n_nuclei/150,c=CC[dis],label=dis,lw=.4,ec='white')
    ax.set_xlabel('RNA genes / nucleus (×10³)');ax.set_ylabel('ATAC peaks / nucleus (×10³)')
    ax.legend(loc='upper left',handletextpad=.2)
    save(f,'Fig1',{'a':'Datasets used in the study. Counts come from the existing study description.','b':'WNN coordinates from the reanalysis. Colours and labelled counts show eight broad cell types. Black points mark M-6.','c':'Cell-type make-up of the nuclei in each donor, with the yield printed above each bar. Squares below mark RA (vermilion) and OA (grey).','d':'RNA weights by subtype. Points are subtype summaries, black ticks are cell-type medians and lines run from the lowest to the highest subtype value.','e':'Raw benchmark metrics on a 0–1 scale. All columns use the same colour scale (white = 0, saturated blue = 1). The overall score is taken as given by the analysis.','f':'Percentage of variance explained by each of the first eight MOFA+ factors.','g':'Median genes versus median peaks for each donor. Point area is proportional to the number of retained nuclei. RNA and ATAC QC counts may differ slightly.'})

def fig2():
    f=page(); sl=read('sldsc_results'); ct=sl[sl.level=='cell type'].copy();ct['ct']=ct.annotation.str[3:]
    ax=standard(f,0,0,'a','Cell-type heritability',left=15)
    z=ct.pivot(index='ct',columns='gwas',values='tau_z').reindex(columns=GW)
    z=z.loc[z[GW[0]].sort_values(ascending=False).index]
    pv=ct.pivot(index='ct',columns='gwas',values='p_tau');qv=ct.pivot(index='ct',columns='gwas',values='fdr_tau')
    for j,g in enumerate(GW):
        for i,k in enumerate(z.index):
            p=max(pv.loc[k,g],1e-8);size=8+10*min(-np.log10(p),6)
            im=ax.scatter(j,i,s=size,c=[z.loc[k,g]],cmap=DIV,vmin=-5,vmax=5,ec='black' if qv.loc[k,g]<.05 else 'none',lw=.7)
    ax.set_xticks(range(3),['RA\nEUR','RA\nFinnGen','T2D\ncontrol']);ax.set_yticks(range(len(z)),[CN[k] for k in z.index]);ax.invert_yaxis();clean(ax)
    ax.set_xlim(-.6,2.6);bar_key(ax,im,'τ z')
    note(ax,'Area: −log₁₀ P; ring: FDR < 0.05',y=-.18)
    ax=standard(f,0,1,'b','Top genetically associated states',left=24)
    st=sl[(sl.level=='subtype') & sl.gwas.isin(GW[:2])]
    zz=st.pivot(index='annotation',columns='gwas',values='tau_z')
    top=zz.mean(axis=1).nlargest(10).index;zz=zz.loc[top,GW[:2]];zz.columns=['EUR','FinnGen']
    sd=read('scdrs_groups'); sub=sd[(sd.level=='subtype') & sd.trait.isin(GW[:2])].copy()
    safe=lambda g:''.join(c if c.isalnum() else '_' for c in g).strip('_')
    sub['annotation']=sub.group.map(lambda g:'ST_'+safe(g))
    name=sub.drop_duplicates('annotation').set_index('annotation').group
    im=matrix(ax,zz, vmin=-5,vmax=5,annotate=True,labels=[name.get(k,k).replace(': ',' ').replace('CD4+ ','').replace('HLA-DR+IgG+ plasmablast','plasmablast').replace('CD25-low','CD25low').replace('CD25-high','CD25high') for k in top]);bar_key(ax,im,'τ z')
    ax=standard(f,0,2,'c','Fine-mapped risk in accessible DNA',left=15)
    fm=read('finemap_enrichment');fc=fm[(fm.level=='cell type') & fm.study.isin(['GCST90132223','FINNGEN_R12_M13_RHEUMA','GCST90132224'])]
    for j,(sid,col,lab) in enumerate([('GCST90132223',BLUE,'EUR'),('FINNGEN_R12_M13_RHEUMA',RED,'FinnGen'),('GCST90132224',GREEN,'EAS')]):
        v=fc.set_index(['study','state'])
        for i,k in enumerate(z.index):
            if (sid,k) in v.index:
                r=v.loc[(sid,k)];ax.scatter(r.enrichment,i+(j-1)*.21,s=8+3*min(-np.log10(max(r.p_empirical,1e-6)),5),c=col,
                  marker=['o','s','^'][j],lw=.5,ec='black' if r.fdr<.05 else col)
    ax.axvline(1,c='.6',lw=.5,ls='--');ax.set_xscale('log');ax.set_xticks([.1,1,10],[.1,1,10]);ax.minorticks_off()
    ax.set_yticks(range(len(z)),[CN[k] for k in z.index]);ax.invert_yaxis();ax.set_xlabel('PIP-mass enrichment (× null)')
    ax.legend(handles=[Line2D([],[],marker=m,color=c,ls='',ms=3,label=l) for m,c,l in [('o',BLUE,'EUR'),('s',RED,'FinnGen'),('^',GREEN,'EAS')]],loc='lower right',handletextpad=.2)
    ax=standard(f,1,0,'d','Chromatin–expression convergence',left=10)
    zs=st.groupby('annotation').tau_z.mean(); ss=sub.groupby('group').assoc_mcz.mean()
    cv=pd.DataFrame([{'group':g,'sldsc':zs.get('ST_'+safe(g),np.nan),'scdrs':v} for g,v in ss.items()]).set_index('group')
    lineage={'T':'Tcell','N':'NK','B':'Bcell','M':'myeloid','F':'fibroblast','E':'endothelial'}
    finite=cv.dropna()
    ax.scatter(finite.sldsc,finite.scdrs,s=9,c=[CT.get(lineage.get(k[0],''),GOLD) for k in finite.index],ec='white',lw=.3)
    zero(ax);rho,p=spearmanr(finite.sldsc,finite.scdrs)
    ax.set_xlabel('Chromatin: mean RA τ z');ax.set_ylabel('Expression: mean RA scDRS z')
    for k in finite.nlargest(3,'scdrs').index:
        r=finite.loc[k];ax.annotate(k.split(':')[0],(r.sldsc,r.scdrs),xytext=(2,3),textcoords='offset points',fontsize=5)
    note(ax,f'ρ = {rho:.2f}; {ptext(p)}\nn = {len(finite)} states',y=.98)
    # M-6 has no S-LDSC estimate, so it is left out of this scatter on purpose.
    ax=standard(f,1,1,'e','Cross-cohort fine-mapped variant support',left=24)
    cand=read('finemap_candidates');fv=read('finemap_variants')
    top=cand.sort_values('pip',ascending=False).drop_duplicates('nearest_protein_coding_gene').head(10)
    studies={'GCST90132223':'EUR','FINNGEN_R12_M13_RHEUMA':'FinnGen','GCST90132224':'EAS'}
    P=pd.DataFrame(np.nan,index=top.variant,columns=list(studies.values()))
    for v in P.index:
        for sid,n in studies.items():
            vals=fv.loc[(fv.variant==v)&(fv.study==sid),'pip'];P.loc[v,n]=vals.max() if len(vals) else np.nan
    labs=[f'{r.nearest_protein_coding_gene} · {r.state.split(":")[0]}' for _,r in top.iterrows()]
    im=matrix(ax,P,cmap=SEQ,vmin=0,vmax=1,annotate=True,fmt='.2f',labels=labs);bar_key(ax,im,'PIP')
    ax=standard(f,1,2,'f','B-cell neighbours of CXCL13⁺ T cells',left=10)
    sp=read('tph_spatial')
    for j,(_,r) in enumerate(sp.iterrows()):
        ax.plot([0,1],[r.B_neighbours_CXCL13neg_T,r.B_neighbours_CXCL13pos_T],c='.6',lw=.6)
        ax.scatter([0,1],[r.B_neighbours_CXCL13neg_T,r.B_neighbours_CXCL13pos_T],s=12,c=[CYAN,BLUE],zorder=3,lw=0)

    last=-100
    for _,r in sp.sort_values('B_neighbours_CXCL13pos_T').iterrows():
        yy=max(np.log10(r.B_neighbours_CXCL13pos_T),last+.13);last=yy
        ax.annotate(f'{r.B_neighbours_CXCL13pos_T/r.B_neighbours_CXCL13neg_T:.1f}×',xy=(1,r.B_neighbours_CXCL13pos_T),xytext=(1.12,10**yy),fontsize=5,va='center',arrowprops={'arrowstyle':'-','lw':.4,'color':'.5'})
    ax.set_yscale('log');ax.set_xlim(-.3,1.65);ax.set_xticks([0,1],['CXCL13⁻','CXCL13⁺']);ax.set_xlabel('T-cell neighbourhood (30 µm)');ax.set_ylabel('Mean neighbouring B cells')
    note(ax,f'{len(sp)} RA sections; paired',y=.98)
    ax=standard(f,2,0,'g','RA synovium: full section',left=3,bottom=8)
    xe=spatial_plot_cells()
    pos=xe.CXCL13.astype(str).str.lower().isin(['true','1'])
    cats=[('Other',~xe.cell_type.isin(['T','B','macrophage']),'#C9C9C9',.15),
          ('Macrophage',xe.cell_type=='macrophage',RED,.2),('T, CXCL13⁻',(xe.cell_type=='T')&~pos,CYAN,.25),
          ('B',xe.cell_type=='B',PURPLE,.45),('T, CXCL13⁺',(xe.cell_type=='T')&pos,BLUE,1)]
    # A fixed, reproducible ROI: the densest 500-µm B-cell bin, centred in a 1.2-mm square.
    bs=xe[xe.cell_type=='B'];edgesx=np.arange(xe.x.min(),xe.x.max()+500,500);edgesy=np.arange(xe.y.min(),xe.y.max()+500,500)
    hist,ex,ey=np.histogram2d(bs.x,bs.y,bins=[edgesx,edgesy]);ii,jj=np.unravel_index(hist.argmax(),hist.shape)
    xc,yc=(ex[ii]+ex[ii+1])/2,(ey[jj]+ey[jj+1])/2;half=600
    for lab,m,col,size in cats:ax.scatter(xe.x[m],-xe.y[m],s=size,c=col,lw=0,rasterized=True,label=lab)
    ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([]);clean(ax)
    ax.add_patch(Rectangle((xc-half,-yc-half),2*half,2*half,fill=False,lw=.7,ec='black'))
    ax.text(xc+half,-yc+half,'h',fontsize=7,fontweight='bold',ha='left',va='bottom')
    x0,y0=xe.x.min(),-xe.y.max();ax.plot([x0,x0+1000],[y0-100,y0-100],c='black',lw=1);ax.text(x0+500,y0-250,'1 mm',ha='center',va='top',fontsize=5)
    ax=standard(f,2,1,'h','B-cell aggregate: linked 1.2-mm ROI',left=3,bottom=10)
    roi=(abs(xe.x-xc)<=half)&(abs(xe.y-yc)<=half)
    for lab,m,col,size in cats:ax.scatter(xe.x[m&roi],-xe.y[m&roi],s=size*6,c=col,lw=0,rasterized=True,label=lab)
    ax.set_xlim(xc-half,xc+half);ax.set_ylim(-yc-half,-yc+half);ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([]);clean(ax)
    ax.plot([xc-half+50,xc-half+250],[-yc-half+60]*2,c='black',lw=1);ax.text(xc-half+150,-yc-half+100,'200 µm',ha='center',va='bottom',fontsize=5,bbox={'fc':'white','ec':'none','alpha':.85,'pad':.5})
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.02),ncol=3,markerscale=1.5,handlelength=.5,handletextpad=.2,columnspacing=.6,labelspacing=.2,fontsize=5)
    ax=standard(f,2,2,'i','Spatial tests across RA sections',left=22)
    y=np.arange(len(sp));b_log=-np.log10(sp.p_B.clip(lower=1/1001));m_log=-np.log10(sp.p_mac.clip(lower=1/1001));ax.scatter(b_log,y,s=13,c=PURPLE,label='B-cell proximity');ax.scatter(m_log,y,s=13,c=RED,marker='s',label='Macrophage IFN')
    for i,(_,r) in enumerate(sp.iterrows()):ax.plot([m_log.iloc[i],b_log.iloc[i]],[i,i],c='.75',lw=.6,zorder=0)
    ax.axvline(0,c='.7',lw=.5);ax.set_yticks(y,[r.replace('CT_','') for r in sp.block]);ax.set_xlabel('−log₁₀ empirical permutation P');ax.set_xlim(-.1,3.3)
    ax.legend(loc='upper left',bbox_to_anchor=(0,1.08),handletextpad=.3);ax.set_ylim(-.5,len(sp)+.6)
    save(f,'Fig2',{'a':'Colour shows S-LDSC tau z and point area shows nominal −log10 P, capped at 6. A black ring marks FDR < 0.05.','b':'Top 10 subtypes by mean RA S-LDSC statistic in the two European RA cohorts.','c':'Fine-mapped PIP mass relative to a permutation null. Shape and colour show the GWAS, point area shows nominal P and a black outline marks FDR < 0.05.','d':'Mean RA S-LDSC versus mean RA scDRS for states that have both estimates, with the Spearman correlation across states. M-6 has no S-LDSC estimate and is left out.','e':'Top variants, chosen by the highest supplied PIP with one variant per nearest protein-coding gene. Labels include the peak state. Grey dashes mark missing study records, not zero PIP. Exact variant IDs are in the source table saved with this figure.','f':'Paired B-cell neighbourhood means within each section, with the positive/negative fold ratio. Lines join values from the same section.','g':'Xenium 5K section CT_21045-A-15, read from a processed plotting table made by step 20. The outlined square is the ROI shown in h. A cell counts as CXCL13-positive when its raw count is >0.','h':'1.2-mm ROI centred on the densest 500-µm B-cell bin. It is chosen by a reproducible rule to show the structure of the aggregate, not to measure enrichment. Same cells and coordinates as g.','i':'Empirical permutation P for each section, capped at 1/1001, the smallest value possible with 1,000 permutations. Global label shuffling does not keep the spatial clustering of labels or the tissue compartments. These associations do not show causal organisation.'})

def fig3():
    f=page()
    ax=standard(f,0,0,'a','S1 discovery: active vs remission',left=10)
    de=read('emtab8322_DE_active_vs_remission',index=0).dropna(subset=['stat'])
    de=de[~de.index.str.startswith('ENSG')].copy();de['nlp']=-np.log10(np.maximum(2*norm.sf(np.abs(de.stat)),1e-300))
    sig=de.padj<.05
    ax.scatter(de.log2FoldChange,de.nlp,s=.6,c='.75',lw=0,rasterized=True)
    ax.scatter(de.loc[sig,'log2FoldChange'],de.loc[sig,'nlp'],s=7,c=np.where(de.loc[sig,'log2FoldChange']>0,RED,BLUE),lw=0)
    zero(ax,y=False);ax.set_xlim(-7.5,10);ax.set_ylim(0,de.nlp.max()*1.13)
    labels_right(ax,de,'log2FoldChange','nlp',['MT1H','MT1G','MT2A','PKM','SLAMF7','GBP5'])
    ax.set_xlabel('log₂ fold change (active / remission)');ax.set_ylabel('−log₁₀ P (Wald test)')
    note(ax,f'{sig.sum()} FDR < 0.05\n{len(de):,} tested genes')
    ax=standard(f,0,1,'b','S1 signature across myeloid states',left=25,right=11)
    sm=read('signature_by_myeloid_subset',index=0).sort_values('rna_signature_z',ascending=False)
    dat=sm.rename(columns={'rna_signature_z':'RNA','atac_signature_z':'ATAC'})
    im=matrix(ax,dat,vmin=-1.5,vmax=1.5,annotate=True,labels=[k.replace(': ',' ').replace('SELENOP+ ','') for k in sm.index]);bar_key(ax,im,'z')
    ax=standard(f,0,2,'c','S1 donor RNA–chromatin concordance',left=10)
    md=read('multiome_donor_signature',index=0);md=md[md.disease=='RA'].dropna(subset=['rna_signature_score','atac_gene_activity_signature_score'])
    x,y=md.rna_signature_score,md.atac_gene_activity_signature_score
    ax.scatter(x,y,s=16,c=RED,lw=.4,ec='white');coef=np.polyfit(x,y,1);xx=np.linspace(x.min(),x.max(),50);ax.plot(xx,np.polyval(coef,xx),c='black',lw=.6)
    rho,p=spearmanr(x,y);stats=read('multiome_signature_stats');rp=stats.loc[stats.score=='atac_gene_activity_signature_score','p_vs_random_genesets'].iloc[0]
    note(ax,f'n = {len(md)} RA donors\nρ = {rho:.2f}; {ptext(p)}\nRandom gene sets: {ptext(rp)}')
    ax.set_xlabel('Donor RNA signature');ax.set_ylabel('Donor gene-activity signature')
    ax=standard(f,1,0,'d','Disease and treatment pathway effects',left=25)
    g=read('gsea_results');g=g[g.library=='hallmark'];nes=g.pivot(index='pathway',columns='contrast',values='nes');q=g.pivot(index='pathway',columns='contrast',values='fdr')
    ca=next(k for k in nes if k.startswith('A'));cb=next(k for k in nes if k.startswith('B'))
    pathways=nes[ca].nlargest(8).index
    dat=nes.loc[pathways,[ca,cb]];dat.columns=['Disease','Anti-TNF']
    labs=[p.replace('HALLMARK_','').replace('_',' ').lower().replace('interferon gamma response','IFN-γ response').replace('interferon alpha response','IFN-α response').replace('tnfa signaling via nfkb','TNF-α / NF-κB').replace('oxidative phosphorylation','OXPHOS') for p in pathways]
    im=matrix(ax,dat,vmin=-3.5,vmax=3.5,labels=labs)
    for i,pw in enumerate(pathways):
        for j,contrast in enumerate([ca,cb]):
            ax.text(j,i,f'{dat.iloc[i,j]:.1f}'+('*' if q.loc[pw,contrast]<.05 else ''),ha='center',va='center',fontsize=5,color='white' if abs(dat.iloc[i,j])>2.3 else 'black')
    bar_key(ax,im,'NES');note(ax,'* FDR < 0.05',y=-.2)
    ax=standard(f,1,1,'e','IFNG⁺ neighbourhoods (4 RA sections)',left=10)
    xi=read('xenium_ifng_neighbourhood');xi=xi[xi.condition=='RA']
    for j,(_,r) in enumerate(xi.iterrows()):
        yy=[r.coef_IFNGneg_neighbour,r.coef_IFNGpos_neighbour]
        ax.plot([0,1],yy,c='.6',lw=.6);ax.scatter([0,1],yy,s=12,c=[CYAN,BLUE],lw=0,zorder=3)

    last=-100
    for _,r in xi.sort_values('coef_IFNGpos_neighbour').iterrows():
        yy=max(r.coef_IFNGpos_neighbour,last+.012);last=yy
        ax.annotate(ptext(r.p_perm_one_sided),xy=(1,r.coef_IFNGpos_neighbour),xytext=(1.12,yy),fontsize=5,va='center',arrowprops={'arrowstyle':'-','lw':.4,'color':'.5'})
    zero(ax,x=False);ax.set_xticks([0,1],['IFNG⁻','IFNG⁺']);ax.set_xlim(-.25,2.1)
    ax.set_xlabel('T / NK neighbour within 30 µm');ax.set_ylabel('Macrophage IFN score / neighbour')

    ax=standard(f,1,2,'f','S2 IFN-high tissue: Visium',left=10)
    vs=read('visium_sections');groups=[vs.loc[vs.condition==k,'pct_spots_ifn_high'] for k in ['NH','OA','RA']]
    dots(ax,groups,['Healthy','OA','RA'],[BLUE,GREY,RED],'IFN-high Visium spots (%)',sizes=10)
    vt=read('visium_tests');note(ax,f'RA vs healthy: {ptext(vt.p.iloc[0])}\nDonor-level test; OA: 1 donor')
    ax=standard(f,2,0,'g','Hallmark IFN-γ proteins in matched tissue',left=10)
    pr=read('proteomics_region_scores',where=PRO);dm=protein_groups(pr)
    groups=[dm.loc[(dm.region_group==reg)&(dm.Disease==dis),'ifng_response'] for reg in ['lining','sublining'] for dis in ['NH','RA']]
    dots(ax,groups,['Healthy','RA','Healthy','RA'],[BLUE,RED,BLUE,RED],'IFN-γ protein score (z)',sizes=10)
    pt=read('proteomics_tests');pt=pt[pt.set=='ifng_response']
    for x,reg in [(0,'lining'),(2,'sublining')]:
        row=pt.loc[pt.comparison.eq(f'RA vs NH {reg} (donor means)')].iloc[0]
        from scipy.stats import mannwhitneyu
        ra=dm.loc[(dm.region_group==reg)&dm.Disease.eq('RA'),'ifng_response']
        nh=dm.loc[(dm.region_group==reg)&dm.Disease.eq('NH'),'ifng_response']
        assert len(ra)==int(row.n_RA) and len(nh)==int(row.n_NH), 'Protein plot/test sample mismatch'
        assert np.isclose(mannwhitneyu(ra,nh).pvalue,row.p_mannwhitney,atol=5e-6), 'Protein plot/test P mismatch'
        dm.to_csv(OUT/'Fig3g_protein_donor_means.tsv',sep='\t',index=False)
        ax.text(x+.5,1.02,f'{reg.capitalize()}: {ptext(row.p_mannwhitney)}',transform=ax.get_xaxis_transform(),ha='center',fontsize=5)
    ax=standard(f,2,1,'h','S2 IFN programme by disease stage',left=10)
    st=read('gse89408_sample_scores',where=PRO,index=0)
    order=['healthy','arthralgia','UA','early RA','established RA','OA']
    groups=[st.loc[st.group==k,'S2'] for k in order]
    dots(ax,groups,['Healthy','Arthr.','UA','Early\nRA','Estab.\nRA','OA'],[BLUE,CYAN,GOLD,RED,PURPLE,GREY],'S2 IFN programme (z)',sizes=1.8,showbox=True)
    gt=read('gse89408_tests');pv=gt.loc[(gt.score=='S2')&(gt.test=='H2 RA vs OA'),'p_vs_random'].iloc[0]
    note(ax,f'RA vs OA: random-set {ptext(pv)}')
    ax=standard(f,2,2,'i','S1 signature change on anti-TNF',left=10)
    ss=read('gse198520_sample_scores',where=PRO,index=0)
    w=ss.pivot_table(index='patient',columns='timepoint',values='signature_score_tech_resid').join(ss.groupby('patient').response.first()).dropna(subset=['Pre','Post'])
    w['delta']=w.Post-w.Pre
    groups=[w.loc[w.response==k,'delta'] for k in ['non','moderate','good']]
    dots(ax,groups,['Non','Moderate','Good'],[GREY,GOLD,RED],'Δ S1 signature (week 12 − baseline)',sizes=5,showbox=True)
    zero(ax,x=False)
    tt=read('gse198520_tests')
    for i,k in enumerate(['non','moderate','good']):
        row=tt[(tt.hypothesis=='H3-adjusted') & (tt.test==f'Wilcoxon technical-residual score change ({k})')].iloc[0]
        ax.text(i,1.02,ptext(row.p),transform=ax.get_xaxis_transform(),ha='center',fontsize=5)
    save(f,'Fig3',{'a':'DESeq2 active/remission effects and Wald statistics from the original analysis. P is recomputed from the stored Wald statistic so that no P value is rounded to zero. FDR status uses the supplied adjusted P.','b':'S1 is defined by 15 discovery genes. The 13 measurable genes give the RNA and ATAC gene-activity z scores shown for all 14 myeloid states, on one shared centred scale.','c':'S1 discovery-signature scores for each donor from RNA and ATAC gene activity (13 of 15 genes measurable), with the Spearman correlation and the supplied random-gene-set test. The line is a descriptive least-squares fit.','d':'The eight Hallmark pathways with the highest disease NES, with treatment NES for the same pathways. * marks supplied GSEA FDR < 0.05. Disease = active/remission macrophages; treatment = post/pre bulk in good responders.','e':'Paired regression coefficients for each section, for IFNG-positive and IFNG-negative T/NK neighbours. P values are one-sided spatial permutation tests from the existing table.','f':'S2 is the predefined 10-gene IFN programme. Each dot is a section and black bars are medians. The comparison shown is the supplied donor-level RA/healthy test, so section counts are not independent donor counts.','g':'Hallmark IFN-gamma protein score. Sublining joins immune_sublining, fibrous_sublining and sublining before the average within each donor, as in the statistical test (RA n=6, NH n=3 in each compartment). Black bars are medians. Two-sided Mann-Whitney P values are recomputed and checked against the source table. Protein-set specificity is supported in lining but not in sublining (random-set P=0.339).','h':'S2 is the predefined 10-gene IFN programme and is different from S1. Points are single biopsy scores. Boxes show the IQR and median, and whiskers 1.5×IQR. Counts are biopsies. The supplied random-gene-set comparison is shown.','i':'S1 is the 15-gene discovery signature (13 measured in this cohort). Points are the paired change for each patient after technical residualisation. Boxes show the IQR and median, and whiskers 1.5×IQR. The P values shown are supplied paired Wilcoxon tests, not comparisons between response groups.'})

def fig4():
    f=page();drv=read('tf_M6_drivers',index=0);drv=drv[~drv.complex].dropna(subset=['motif_diff_z'])
    ax=standard(f,0,0,'a','Candidate regulatory activity',left=10)
    ax.scatter(drv.rna_activity_diff,drv.motif_diff_z,s=3,c='.75',lw=0)
    hi=drv[drv.both_layers];ax.scatter(hi.rna_activity_diff,hi.motif_diff_z,s=12,c=RED,lw=.3,ec='white')
    zero(ax);ax.set_xlabel('Δ RNA regulon activity (M-6 − other)');ax.set_ylabel('Δ ATAC motif activity (z)')
    for k in ['RELA','REL','STAT1','IRF5','IRF7','IRF9']:
        if k in drv.index:
            r=drv.loc[k];ax.annotate(k,(r.rna_activity_diff,r.motif_diff_z),xytext=(2,3),textcoords='offset points',fontsize=5)
    note(ax,f'{len(hi)} two-layer candidates\nExploratory: 2 donors')
    ax=standard(f,0,1,'b','Exploratory TF evidence (2 donors)',left=11)
    top=drv.loc[drv.rna_rank.add(drv.motif_rank).nsmallest(10).index]
    dat=top[['rna_activity_diff','motif_diff_z']].copy();dat.columns=['RNA Δ','ATAC Δz']
    # The two columns have different units, so each one has its own scale.
    for j,(col,vmax) in enumerate([('RNA Δ',max(7,dat.iloc[:,0].max())),('ATAC Δz',max(6,dat.iloc[:,1].max()))]):
        for i,(g,r) in enumerate(top.iterrows()):
            val=dat.iloc[i,j];ax.scatter(j,i,s=15+2*np.sqrt(r.n_targets),c=[val],cmap=SEQ,vmin=0,vmax=vmax,ec='black' if r.rna_q_BH<.05 else 'none',lw=.5)
            ax.text(j+.22,i,f'{val:.1f}',fontsize=5,va='center')
    ax.set_yticks(range(len(top)),top.index);ax.set_xticks([0,1],['RNA Δ\n0–7','ATAC Δz\n0–6']);ax.invert_yaxis();ax.set_xlim(-.5,1.8);clean(ax)
    note(ax,'Area: target count; ring: RNA BH q < 0.05',y=-.17)
    ax=standard(f,0,2,'c','IRF5 signal-specific colocalisation',left=25)
    cm=read('coloc_multisignal_all');ir=cm[(cm.gene=='IRF5') & cm.cis & (cm.category=='myeloid')].copy()
    cond={'naive':'naive','IFNg_18h':'IFN-γ','LPS_6h':'LPS','Listeria_5h':'Listeria','Salmonella_5h':'Salmonella','R848_6h':'R848','Pam3CSK4_6h':'Pam3CSK4','Influenza_6h':'Flu'}
    ir['ds']=ir.qtl_project.str.split('_').str[0]+' '+ir.biosample.map(lambda b:'mono' if 'monocyte' in b else 'mac' if 'macrophage' in b else b)+' '+ir.condition.map(lambda c:cond.get(c,c) if isinstance(c,str) else '')
    best=ir.groupby(['ra_lead_rsid','ds']).h4.max().unstack(0)
    cols=[s for s in ['rs3778754','rs4731532','rs13246321'] if s in best]
    best=best.loc[best.max(axis=1).nlargest(10).index,cols]
    im=matrix(ax,best,cmap=SEQ,vmin=0,vmax=1,annotate=True,fmt='.2f');ax.set_xticklabels(['rs3778754\nEUR','rs4731532\nEUR','rs13246321\nFinnGen'],rotation=35,ha='right')
    ax=standard(f,1,0,'d','Expression direction at drug-target loci',span=2,left=23,right=25)
    mr=read('mr_primary');keep=['CTLA4','IL6ST','CD40','FCGR2A','CCR6','IRF5'];mr=mr[mr.gene.isin(keep)].copy()
    mr['order']=mr.gene.map({g:i for i,g in enumerate(keep)});mr=mr.sort_values(['order','gwas'],ascending=[False,True]);y=np.arange(len(mr))
    for i,(_,r) in enumerate(mr.iterrows()):
        col=BLUE if r.OR_per_SD<1 else RED
        if r.share_same_sign<.75:col=GREY
        ax.plot([r.OR_lo,r.OR_hi],[i,i],c=col,lw=.8)
        ax.scatter(r.OR_per_SD,i,s=13,c=col,marker='o' if r.gwas.startswith('Ishigaki') else 's',lw=0,zorder=3)
        agree=int(round(r.share_same_sign*r.n_datasets))
        ax.text(1.03,i,f'{agree}/{int(r.n_datasets)}',transform=ax.get_yaxis_transform(),fontsize=5,va='center')
    ax.text(1.03,1.03,'Direction agreement',transform=ax.transAxes,fontsize=5)
    labs=[];dup=mr.groupby(['gene','gwas']).size();seen={}
    for _,r in mr.iterrows():
        key=(r.gene,r.gwas);seen[key]=seen.get(key,0)+1
        labs.append(r.gene+' · '+('EUR' if r.gwas.startswith('Ishigaki') else 'FinnGen')+(f' [{seen[key]}]' if dup[key]>1 else ''))
    ax.set_yticks(y,labs)
    mr.to_csv(OUT/'Fig4d_MR_estimates.tsv',sep='\t',index=False)
    ax.axvline(1,c='.5',lw=.5,ls='--');ax.set_xscale('log');ax.set_xlim(.6,1.8);ax.set_xticks([.7,1,1.5],[.7,1,1.5]);ax.minorticks_off()
    ax.set_xlabel('RA odds ratio per SD higher expression (95% CI)');ax.set_ylim(-.7,len(mr)+1.8)
    note(ax,'Circle: EUR; square: FinnGen\nGrey: <75% same-direction QTL datasets',x=.01,y=.98)
    ax=standard(f,1,2,'e','Risk-allele direction and replication',left=12)
    di=read('coloc_direction',index=0);sel=['CTLA4','IL6ST','CD40','IRF5','ANKRD55','CCR6','FCGR2A','STAT4']
    di=di.reindex(sel).dropna(subset=['frac_risk_raises_expression']);y=np.arange(len(di))
    for i,(gene,r) in enumerate(di.iterrows()):
        col=RED if r.frac_risk_raises_expression>=.75 else BLUE if r.frac_risk_raises_expression<=.25 else GREY
        ax.plot([.5,r.frac_risk_raises_expression],[i,i],c='.75',lw=.6)
        ax.scatter(r.frac_risk_raises_expression,i,s=8+.5*r.immune_datasets,c=col,lw=.5,ec='black' if r.consistent else 'white')
        ax.text(1.1,i,str(int(r.immune_datasets)),fontsize=5,va='center')
    ax.axvline(.5,c='.5',ls='--',lw=.5);ax.set_yticks(y,di.index);ax.invert_yaxis();ax.set_xlim(-.05,1.4);ax.set_xticks([0,.5,1],[0,.5,1]);ax.set_xlabel('Fraction: risk allele raises expression')
    note(ax,'Area / right: QTL datasets',y=-.17)
    ax=standard(f,2,0,'f','Top-gene evidence across eleven layers',span=2,left=21,right=3,bottom=14)
    pz=read('gene_prioritisation',index=0);layers=[c for c in pz if c.startswith('L') and '_' in c][:11];top=pz.head(24)
    dat=top[layers].T;dat.index=['scRNA disease','M-6 RNA marker','Anti-TNF change','Chromatin link','GWAS locus','Variant motif','IFN-γ pathway','Communication','Protein','Colocalisation','Fine-mapping']
    im=matrix(ax,dat,cmap=SEQ,vmin=0,vmax=1)
    ax.set_xticks(range(len(top)),top.index,rotation=60,ha='right',fontstyle='italic')
    for j,(gene,r) in enumerate(top.iterrows()):
        ax.text(j,-1.1,f'{r.evidence_score:g}',ha='center',fontsize=5)
        if str(r.novelty).startswith('new'):ax.plot(j,11.25,'o',ms=1.7,c=GREEN,clip_on=False)
    ax.set_ylim(10.5,-1.8);ax.text(-.8,-1.1,'Score',ha='right',fontsize=5)
    # Green markers sit on a separate track, so the gene labels stay black.
    ax=standard(f,2,2,'g','Evidence and prior RA knowledge',left=10)
    ax.scatter(pz.ot_ra_score,pz.evidence_score,s=5+pz.n_layers_with_evidence*2,c=[GREEN if str(v).startswith('new') else GREY for v in pz.novelty],alpha=.8,lw=.3,ec='white')
    for gene in ['CD40','STAT4','SLAMF7','CD38','IL4I1']:
        if gene in pz.index:
            r=pz.loc[gene];ax.annotate(gene,(r.ot_ra_score,r.evidence_score),xytext=(2,3),textcoords='offset points',fontsize=5,fontstyle='italic')
    ax.set_xlabel('Open Targets RA association score');ax.set_ylabel('Integrated evidence score')
    ax.set_ylim(-.3,7.6);note(ax,'Green: little prior RA evidence\nArea: number of evidence layers',y=.98)
    save(f,'Fig4',{'a':'Differences in transcription-factor RNA regulon activity and chromVAR motif activity. Candidates pass BH q <0.05 across all tested RNA gene-label-permutation P values and rank in the top motif decile. Data come from two donors, and gene-label permutations are not donor replication. These are candidate regulators.','b':'Top 10 factors by summed RNA and motif rank. Point area = 15 + 2 sqrt(regulon target count). Each column has its own labelled colour scale in its own units (RNA 0–7, ATAC 0–6). A black outline marks RNA gene-label-permutation BH q <0.05. Printed numbers are unscaled effects.','c':'Highest IRF5 H4 for each project/biosample/condition and RA lead signal. The top 10 conditions by highest H4 are shown. Grey dashes mark missing records.','d':'Existing MR odds ratios and supplied 95% CI for six genes. All supplied estimates are shown, including more than one lead signal for some gene/GWAS pairs; numbers in brackets tell these apart, in source-table order. The right track gives same-direction/total datasets. Colour shows direction and grey marks agreement <75%. These are exploratory estimates, not measures of drug efficacy.','e':'Fraction of immune QTL datasets in which the risk allele raises expression. Point area and the number on the right give the number of datasets. A black ring marks the supplied consistency flag.','f':'Top 24 genes across the 11 original evidence layers, on a shared 0–1 scale. The original evidence score is printed above. White = 0, mid blue = 0.5, saturated blue = 1. A green marker below a gene shows low prior RA evidence, based on the supplied novelty label. The evidence layers are not assumed to be statistically independent.','g':'All prioritised genes: original integrated score versus Open Targets association score. Point area = 5 + 2 × number of evidence layers. Green marks genes with the supplied low-prior-evidence label.'})

def ext1():
    f=page(169)
    ax=cell(f,'a','Classification against matched baselines',2,2,87,54,left=23,bottom=9)
    cl=read('classifier_baselines');models=['F0 covariates','F1 eight programme scores','F2 top-30 prioritised genes','F3 genome-wide elastic net','F4 random 30-gene sets (median of 100)']
    short=['Covariates','8 programmes','Top 30 genes','Genome-wide EN','Random 30 genes']
    for j,(task,col) in enumerate([('T1 RA vs OA',RED),('T2 ACPA+ vs ACPA- (RA)',BLUE)]):
        x=cl[cl.task==task].set_index('model').reindex(models);y=np.arange(5)+(j-.5)*.22
        ax.errorbar(x.auc_mean,y,xerr=x.auc_sd,fmt='o',ms=3,c=col,elinewidth=.7,label=task.split(' ',1)[1])
    ax.set_yticks(range(5),short);ax.invert_yaxis();ax.set_xlim(.2,1.05);ax.axvline(.5,c='.6',lw=.5,ls='--');ax.set_xlabel('AUC ± SD, 5 × 5 nested CV');ax.legend(loc='lower left',handletextpad=.3)
    ax=cell(f,'b','Anti-TNF response prediction and confounding',93,2,88,54,left=36,bottom=9)
    au=read('gse198520_auc');y=np.arange(len(au))
    ax.errorbar(au.auc_good_vs_non,y,xerr=[au.auc_good_vs_non-au.ci_low,au.ci_high-au.auc_good_vs_non],fmt='o',ms=3,c=BLUE,elinewidth=.7)
    ax.axvline(.5,c='.6',lw=.5,ls='--');ax.set_yticks(y,au.predictor);ax.invert_yaxis();ax.set_xlabel('AUC (95% CI)');ax.set_xlim(0,1)

    ax=cell(f,'c','Robustness across perturbations',2,58,64,54,left=30,bottom=9)
    rb=read('robustness_summary');labs=['Signature FDR retention','Signature direction','Genome-wide agreement','Inflammatory fraction','Homeostatic score','RNA–ATAC concordance','Transfer / histology','T-cell enrichment','GWAS empirical P','Pathway concordance']
    if len(labs)!=len(rb):labs=rb.analysis.str.slice(0,22).tolist()
    y=np.arange(len(rb));ax.barh(y,rb.fraction_meeting_rule,color=GREEN,height=.5)
    for i,r in rb.iterrows():ax.text(1.03,i,str(int(r.n_perturbations)),fontsize=5,va='center')
    ax.text(1.03,1.03,'Runs',transform=ax.transAxes,fontsize=5);ax.set_yticks(y,labs);ax.invert_yaxis();ax.set_xlim(0,1);ax.set_xlabel('Fraction meeting original rule')
    ax=cell(f,'d','Drug-class reversal versus null',69,58,53,54,left=9,bottom=10)
    l=read('l1000_reversal_class_enrichment');x=np.arange(len(l));ax.bar(x-.16,l.observed_share,.3,color=RED,label='Observed');ax.bar(x+.16,l.random_median,.3,color=GREY,label='Random median')
    for i,r in l.iterrows():ax.text(i,max(r.observed_share,r.random_median)+.025,ptext(r.empirical_p),fontsize=5,ha='center')
    ax.set_xticks(x,['Active /\nremission','IFN-γ /\nhomeostatic']);ax.set_ylabel('Anti-inflammatory class share');ax.legend(loc='upper right');ax.set_ylim(0,max(l.observed_share.max(),l.random_median.max())*1.4)
    ax=cell(f,'e','M-6 graph neighbourhood',125,58,56,54,left=11,bottom=9)
    pg=read('myeloid_paga',index=0).drop(index='M-6').sort_values('rna_paga_to_M6',ascending=False);y=np.arange(len(pg))
    ax.scatter(pg.rna_paga_to_M6,y,s=11,c=RED,label='RNA');ax.scatter(pg.atac_paga_to_M6,y,s=11,c=BLUE,marker='s',label='ATAC')
    for i,(_,r) in enumerate(pg.iterrows()):ax.plot([r.rna_paga_to_M6,r.atac_paga_to_M6],[i,i],c='.75',lw=.6,zorder=0)
    ax.set_yticks(y,pg.index);ax.invert_yaxis();ax.set_xlabel('PAGA connectivity to M-6');ax.legend(loc='lower right',handletextpad=.3)
    ax=cell(f,'f','Boundary findings and negative controls',2,116,179,51,left=4,top=8,bottom=2)
    ax.set_axis_off()
    blood=read('blood_tests');bp=blood.loc[blood.test.str.contains('paired Wilcoxon'),'p'].iloc[0]
    items=[('Blood monocytes',f'Paired signature test: P = {bp:.2f}; 13 pairs'),
           ('Bulk deconvolution','Failed the original sanity checks'),
           ('Active / remission classification','AUC 0.50–0.60 in the original comparison'),
           ('RA / OA: prioritised versus random genes',cl.loc[(cl.task=='T1 RA vs OA') & cl.model.str.startswith('F4'),'note'].iloc[0].split(';')[0]),
           ('Bulk Tph score, adjusted for T-cell content','Random-set P = 0.81'),
           ('Myeloid S-LDSC','Negative conditional τ z in both RA GWAS'),
           ('Trajectory timing','Chromatin-before-RNA is not interpretable'),
           ('East Asian T-cell peak enrichment','1.26×; P = 0.22')]
    for i,(label,result) in enumerate(items):
        col=i//4;row=i%4;xx=col*.51;yy=.98-row*.24
        ax.text(xx,yy,label,fontsize=6,fontweight='bold',va='top')
        ax.text(xx,yy-.10,result,fontsize=5.5,va='top')
        ax.plot([xx,xx+.47],[yy-.19,yy-.19],c='.85',lw=.5)
    ax.set_xlim(0,1);ax.set_ylim(0,1)
    save(f,'ExtData_Fig1',{'a':'Mean AUC ± fold SD from the revised repeated nested CV, with all fitting done on training data only. Error bars are SD, not CI. For random controls, the error is the SD across the means of the draws. Observed and label-permuted F1 use the same 5 × 5 outer CV, 3-fold tuning and preprocessing. These are internal estimates from one cohort.','b':'AUC for good response versus nonresponse at baseline, with the supplied bootstrap 95% CI. The predictive signal seen at baseline is confounded by technical factors.','c':'Fraction of perturbations that meet each original rule. The right column gives the number of perturbations. Metric definitions and ranges are in robustness_summary.tsv.','d':'Observed share of anti-inflammatory compound classes versus the random median. P values are supplied empirical null comparisons.','e':'Existing RNA and ATAC PAGA connectivity to M-6. Closeness in the graph does not show developmental or temporal direction.','f':'Boundary findings. The blood P and the revised classifier random-set P are read from their result tables. The other entries come from the named source analyses.'})

def paired(ax, df, col, colour_by=None, ylabel=''):
    """Draw one before/after line per patient. df is indexed by unit and has columns patient, treatment/origin and col."""
    for p, g in df.groupby('patient'):
        g = g.set_index('state')
        if {'before', 'after'} <= set(g.index):
            c = colour_by.get(p, GREY) if colour_by else GREY
            ax.plot([0, 1], [g.loc['before', col], g.loc['after', col]], c=c, lw=.7, alpha=.85, marker='o', ms=2.2)
    w = df.pivot_table(index='patient', columns='state', values=col).dropna()
    for i, s in enumerate(['before', 'after']):
        ax.plot([i - .18, i + .18], [w[s].median()] * 2, c='black', lw=1.1, zorder=5)
    lo, hi = np.nanmin(df[col]), np.nanmax(df[col]); r = hi - lo
    ax.set_xlim(-.4, 1.4); ax.set_ylim(lo - .05 * r, hi + .42 * r); ax.set_ylabel(ylabel)

def ext2():
    f = page(169)
    rt = read('ra_fluid_tests').set_index('test'); nt = read('nsclc_extension_tests').set_index('test')
    def row(t, key): return t[t.index.str.startswith(key)].iloc[0]
    u = read('ra_fluid_unit_scores', index=0, where=PRO).copy()
    u['state'] = u.treatment.map({'BT': 'before', 'AT': 'after'})
    drug_col = {p: (RED if d == 'Tofacitinib' else BLUE) for p, d in u.groupby('patient').drug.first().items()}

    ax = standard(f, 0, 0, 'a', 'RA synovial fluid: macrophage M-6 programme on treatment', left=11)
    paired(ax, u.dropna(subset=['M6']), 'M6', drug_col, 'M-6 programme (z)')
    ax.set_xticks([0, 1], ['Naive', 'On treatment'])
    r = row(rt, 'V2 macrophage M6')
    note(ax, f"Fell in {r.detail.split(';')[0].split()[-1]} patients\nPaired {ptext(r.p)}\nvs 200 random sets {ptext(r.p_vs_random)}", x=.98, ha='right')
    ax.legend(handles=[Line2D([], [], c=BLUE, lw=.8, label='Adalimumab'), Line2D([], [], c=RED, lw=.8, label='Tofacitinib')],
              loc='upper left', handlelength=1.2)
    ax = standard(f, 0, 1, 'b', 'Share of M-6-like fluid macrophages', left=11)
    v = u.dropna(subset=['m6_like_share']).copy(); v['pct'] = 100 * v.m6_like_share
    paired(ax, v, 'pct', drug_col, 'M-6-like macrophages (%)'); ax.set_xticks([0, 1], ['Naive', 'On treatment'])
    r = row(rt, 'V3'); note(ax, f"Fell in {r.detail.split(';')[0].split()[-1]} patients\nPaired {ptext(r.p)}", x=.98, ha='right')
    ax = standard(f, 0, 2, 'c', 'Tph programme of fluid T cells', left=11)
    paired(ax, u.dropna(subset=['tcell_Tph']), 'tcell_Tph', drug_col, 'Tph programme (z)'); ax.set_xticks([0, 1], ['Naive', 'On treatment'])
    r = row(rt, 'V4'); note(ax, f"Fell in {r.detail.split()[-1]} patients\nPaired {ptext(r.p)}", x=.98, ha='right')

    ax = standard(f, 1, 0, 'd', 'NSCLC tumour macrophages: M-6-like cluster markers', left=13)
    mk = read('luca_m6_cluster_markers').head(14).iloc[::-1]
    ax.barh(range(len(mk)), mk.logfoldchanges, color=[RED if x else GREY for x in mk.in_programme], height=.65)
    ax.set_yticks(range(len(mk)), mk.gene); ax.set_xlabel('log2 fold change vs other macrophages')
    r = row(nt, 'N1'); note(ax, f"{r.detail.split(';')[0].split(': ')[1]}\nvs random sets {ptext(r.p)}\nvermilion: RA M-6 genes", x=.98, y=.50, ha='right')
    ax = standard(f, 1, 1, 'e', 'NSCLC: macrophage M-6 programme, tumour vs adjacent lung', left=11)
    lu = read('luca_unit_scores', index=0, where=PRO).copy()
    lu = lu[lu.origin.isin(['normal_adjacent', 'tumor_primary'])].rename(columns={'donor_id': 'patient'})
    lu['state'] = lu.origin.map({'normal_adjacent': 'before', 'tumor_primary': 'after'})
    both = lu.groupby('patient').state.nunique(); lu = lu[lu.patient.isin(both[both == 2].index)]
    paired(ax, lu, 'macrophage_M6', None, 'M-6 programme (z)'); ax.set_xticks([0, 1], ['Adjacent\nnormal', 'Primary\ntumour'])
    r = row(nt, 'N2'); note(ax, f"Higher in tumour: {r.detail.split()[-2]}\nPaired {ptext(r.p)}\nvs random sets {ptext(r.p_vs_random)}", x=.02)
    ax = standard(f, 1, 2, 'f', 'NSCLC: macrophage M-6 vs CD4 T-cell Tph', left=11)
    tt = read('luca_unit_scores', index=0, where=PRO)
    tt = tt[(tt.origin == 'tumor_primary')].dropna(subset=['cd4_Tph'])
    rx = tt.macrophage_M6 - tt.groupby('study').macrophage_M6.transform('mean')
    ry = tt.cd4_Tph - tt.groupby('study').cd4_Tph.transform('mean')
    ax.scatter(rx, ry, s=7, c=BLUE, lw=0, alpha=.85); zero(ax)
    ax.set_xlabel('Macrophage M-6 (study-centred)'); ax.set_ylabel('CD4 T Tph (study-centred)')
    r = row(nt, 'N3'); note(ax, f"ρ = {r.statistic:.2f}, {ptext(r.p)}\nvs random sets {ptext(r.p_vs_random)}\n{int(r.n)} tumours, {r.detail}")

    ax = cell(f, 'g', 'TCGA lung cancer: programme scores and overall survival', 2, 112, 92, 55, left=26, bottom=10)
    cx = read('nsclc_extension_cox')
    cx['label'] = cx.cancer + ' ' + cx.score.str.replace('M6', 'M-6') + np.where(cx.adjusted_for_PTPRC, ' + PTPRC', '')
    y = np.arange(len(cx))[::-1]
    ax.errorbar(cx.HR_per_SD, y, xerr=[cx.HR_per_SD - cx.CI_low, cx.CI_high - cx.HR_per_SD], fmt='o', ms=2.5,
                c=BLUE, elinewidth=.6, ecolor=BLUE)
    for yy, (_, r) in zip(y, cx.iterrows()):
        ax.text(2.12, yy, ptext(r.p), fontsize=5, va='center')
    ax.axvline(1, c='.6', lw=.5, ls='--'); ax.set_yticks(y, cx.label); ax.set_xlim(.6, 2.55); ax.set_xticks([.6, 1, 1.4, 1.8])
    ax.set_xlabel('Hazard ratio per SD (age and stage adjusted)')
    ax = cell(f, 'h', 'Pre-specified tests and verdicts', 97, 112, 84, 55, left=2, top=9, bottom=2)
    ax.set_axis_off()
    items = [('RA fluid, 9 patients (GSE296117)', None),
             ('V1 M-6-like cluster exists', f"borderline ({ptext(row(rt, 'V1').p)})"),
             ('V2 M-6 programme falls on treatment', f"vs random {ptext(row(rt, 'V2 macrophage M6').p_vs_random)}; paired {ptext(row(rt, 'V2 macrophage M6').p)}"),
             ('V3 M-6-like share falls', f"supported ({ptext(row(rt, 'V3').p)})"),
             ('V4 Tph programme falls', f"borderline ({ptext(row(rt, 'V4').p)})"),
             ('NSCLC, LuCA atlas + TCGA', None),
             ('N1 M-6-like tumour macrophages', f"supported ({ptext(row(nt, 'N1').p)})"),
             ('N2 higher in tumour than adjacent lung', f"supported ({ptext(row(nt, 'N2').p)})"),
             ('N3 couples with CD4 Tph', f"suggestive (ρ {row(nt, 'N3').statistic:.2f}, {ptext(row(nt, 'N3').p)})"),
             ('N4 bulk survival', 'no association; LUAD HR 1.6 after PTPRC')]
    yy = .97
    for label, res in items:
        if res is None:
            ax.text(0, yy, label, fontsize=6, fontweight='bold', va='top'); yy -= .09; continue
        ax.text(.02, yy, label, fontsize=5.5, va='top'); ax.text(.98, yy, res, fontsize=5.5, va='top', ha='right', color='#333333')
        yy -= .085
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    save(f, 'ExtData_Fig2', {'a': 'M-6 programme score (30 genes defined in RA synovial tissue) in synovial-fluid macrophage pseudobulk for each patient, treatment-naive and on treatment. Lines join the same patient, colour = drug and black bars are medians. The paired Wilcoxon test and the comparison with 200 expression-matched random gene sets come from ra_fluid_tests.tsv.',
                              'b': 'Share of macrophages above the 90th percentile of treatment-naive M-6 cell scores, for each patient and time point.',
                              'c': 'Tph programme score in T-cell pseudobulk. Eight patients have enough T cells at both time points.',
                              'd': 'Top markers of the Leiden cluster with the highest M-6 score among 39,260 primary-tumour macrophages (LuCA core atlas, Harmony by study). Red marks genes in the RA M-6 programme. P compares this cluster with the best cluster of 200 random gene sets.',
                              'e': 'Macrophage pseudobulk M-6 score for each donor in adjacent normal lung and primary tumour, for the 29 donors with >= 50 macrophages in both.',
                              'f': 'Macrophage M-6 versus CD4 T-cell Tph programme score for each tumour, each centred within study. Spearman correlation on the study-adjusted residuals.',
                              'g': 'Cox proportional hazards per SD of the bulk programme score in TCGA primary tumours, adjusted for age and stage, with and without adjustment for leukocyte abundance (PTPRC). Bars are 95% CI.',
                              'h': 'Hypotheses as set in scripts/40_ra_fluid_treatment.py and 41_nsclc_extension.py before scoring. Verdicts use the stated two-sided tests.'})

def main():
    fig1();fig2();fig3();fig4();ext1();ext2()
    (OUT/'panel_guide.json').write_text(json.dumps(PANEL_GUIDE,ensure_ascii=False,indent=2))
    notes=['# Figure panel guide', '', 'Panel notes for the Nature-style data figures, written by scripts/39_main_figures.py.', '']
    for name,panels in PANEL_GUIDE.items():
        notes.extend(['## '+name, ''])
        notes.extend('- **'+letter+'**: '+description for letter,description in panels.items())
        notes.append('')
    (OUT/'Figure_notes.md').write_text('\n'.join(notes))

    (OUT/'source_hashes.json').write_text(json.dumps(SOURCES,indent=2))
    # These composite figures are the project's final figure set. Single analysis steps
    # may write diagnostic PNGs. Remove those only after all final figures have saved without error.
    if OUT.resolve() == (ROOT / 'results/figures/main').resolve():
        for old in (ROOT / 'results/figures').iterdir():
            if old.is_file() and old.suffix.lower() in {'.png', '.pdf', '.svg', '.jpg', '.jpeg', '.tif', '.tiff'}:
                old.unlink()
    # Save the exact variant IDs to go with the short gene/state labels in Fig2e.
    read('finemap_candidates').sort_values('pip',ascending=False).drop_duplicates('nearest_protein_coding_gene').head(10).to_csv(OUT/'Fig2e_variants.tsv',sep='\t',index=False)

if __name__=='__main__':main()
