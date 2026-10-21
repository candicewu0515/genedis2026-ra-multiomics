"""Single-cell disease relevance (scDRS): which nuclei, and which small states S-LDSC could not test (M-6), express
RA polygenic risk genes more than matched control gene sets?

Gene-level GWAS scores: MAGMA v1.10 (SNP-wise mean model; MHC excluded; --seed 1, without which a few gene p-values
vary between runs), SNPs assigned to NCBI GRCh38 genes with a 10-kb window,
LD from the 1000 Genomes EUR GRCh38 panel (per chromosome), full summary statistics with effective N
(Ishigaki 2022 RA EUR, FinnGen R12 RA, FinnGen R12 T2D as control).
scDRS (Zhang et al. 2022 Nat Genet): top 1,000 genes by MAGMA z outside the extended MHC (chr6:25-34 Mb, whose
long-range LD lifts flanking genes such as ITPR3 and ZSCAN31), weights = z (capped at 10); GSE243917 nuclei
(log-normalised), covariates log counts and mito fraction; 1,000 control gene sets matched on mean and variance;
group tests (Monte Carlo) for author cell types and subtypes; heterogeneity within groups.
Hypotheses fixed before analysis:
  D1 T cells (and Tph/Treg subtypes) are associated with RA (both GWAS), not with T2D.
  D2 M-6 STAT1+ CXCL10+ macrophages: associated or not (the state S-LDSC could not test)?
  D3 within T cells, per-nucleus RA scores are heterogeneous (scDRS heterogeneity test).

Output: data/interim/magma/*, results/tables/magma_genes_<trait>.tsv, scdrs_groups.tsv; results/figures/scdrs.png
"""
import os
import hashlib,json
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from config import ROOT, RAW, INTERIM, FIG, TAB, SEED

MAGMA = ROOT / "tools" / "magma" / "magma"
REF = RAW / "ldsc" / "GRCh38" / "plink_files"
GENELOC = RAW / "magma" / "NCBI38.gene.loc"
W = INTERIM / "magma"
TRAITS = {"RA_Ishigaki_EUR": dict(file=RAW / "gwas_sumstats" / "GCST90132223.h.tsv.gz", fmt="gwascat", ncas=22350, ncon=74823),
          "RA_FinnGen_R12": dict(file=RAW / "gwas_sumstats" / "finngen_R12_M13_RHEUMA.gz", fmt="finngen", ncas=16314, ncon=315115),
          "T2D_FinnGen_R12_control": dict(file=RAW / "gwas_sumstats" / "finngen_R12_T2D.gz", fmt="finngen", ncas=82878, ncon=403489)}
N_TOP, N_CTRL, N_JOBS = 1000, 1000, 4
M6 = "M-6: STAT1+ CXCL10+"


def run(cmd, log):
    env = dict(os.environ, OMP_NUM_THREADS="2")
    with open(log, "w") as fh:
        return subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env).returncode


def annotate():
    out = W / "annot.genes.annot"
    if out.exists():
        return out
    W.mkdir(parents=True, exist_ok=True)
    loc = W / "snp.loc"
    with open(loc, "w") as fo:
        for c in range(1, 23):
            b = pd.read_csv(REF / f"1000G.EUR.hg38.{c}.bim", sep="\t", header=None, usecols=[0, 1, 3])
            b[[1, 0, 3]].to_csv(fo, sep="\t", header=False, index=False)
    if run([str(MAGMA), "--annotate", "window=10,10", "--snp-loc", str(loc), "--gene-loc", str(GENELOC),
            "--out", str(W / "annot")], str(W / "annot.runlog")):
        raise RuntimeError("MAGMA annotation failed")
    return out


def pval_file(name, spec):
    out = W / f"{name}.pval"
    if out.exists():
        return out
    keep = []
    for d in pd.read_csv(spec["file"], sep="\t", chunksize=2_000_000, dtype=str):
        if spec["fmt"] == "gwascat":
            x = pd.DataFrame({"SNP": d["rsid"], "P": pd.to_numeric(d.p_value, errors="coerce")})
        else:
            x = pd.DataFrame({"SNP": d.rsids.str.split(",").str[0], "P": pd.to_numeric(d.pval, errors="coerce")})
        keep.append(x[x.SNP.str.startswith("rs", na=False) & x.P.between(0, 1)])
    s = pd.concat(keep).drop_duplicates("SNP")
    b6 = pd.read_csv(REF / "1000G.EUR.hg38.6.bim", sep="\t", header=None, usecols=[1, 3], names=["SNP", "BP"])
    mhc = set(b6.SNP[b6.BP.between(28_510_120, 33_480_577)])
    n_tiny = int(((s.P < 1e-300) & ~s.SNP.isin(mhc)).sum())
    s = s[~s.SNP.isin(mhc)].assign(P=lambda d: d.P.clip(lower=1e-300))   # MHC excluded as elsewhere; MAGMA cannot parse p < 1e-308
    print(f"  {name}: {len(s):,} SNPs after MHC exclusion; non-MHC p < 1e-300 clipped: {n_tiny}")
    s.to_csv(out, sep="\t", index=False)
    return out


def gene_analysis(name, spec, annot):
    final = TAB / f"magma_genes_{name}.tsv"
    if final.exists():
        return pd.read_csv(final, sep="\t")
    pv = pval_file(name, spec)
    neff = int(round(4 / (1 / spec["ncas"] + 1 / spec["ncon"])))
    jobs = [([str(MAGMA), "--seed", "1", "--bfile", str(REF / f"1000G.EUR.hg38.{c}"), "--pval", str(pv), f"N={neff}",
              "--gene-annot", str(annot), "--out", str(W / f"{name}.chr{c}")], str(W / f"{name}.chr{c}.runlog"))
            for c in range(1, 23) if not (W / f"{name}.chr{c}.genes.out").exists()]
    with ThreadPoolExecutor(N_JOBS) as ex:
        codes = list(ex.map(lambda j: run(*j), jobs))
    if any(codes):
        raise RuntimeError(f"MAGMA gene analysis failed for {name}; see {W}")
    g = pd.concat([pd.read_csv(W / f"{name}.chr{c}.genes.out", sep=r"\s+") for c in range(1, 23)])
    sym = pd.read_csv(GENELOC, sep="\t", header=None, names=["GENE", "CHR", "START", "STOP", "STRAND", "SYMBOL"]).set_index("GENE").SYMBOL
    g["SYMBOL"] = g.GENE.map(sym)
    g = g.sort_values("ZSTAT", ascending=False)
    g.to_csv(final, sep="\t", index=False)
    print(f"  MAGMA {name}: {len(g)} genes; top: {', '.join(g.SYMBOL.head(12))}")
    # Use the same canonical serialized input whether MAGMA was just run or cached.
    return pd.read_csv(final, sep="\t")


def main():
    if not hasattr(np, "float_"):          # scDRS 1.0.2 still uses np.float_, removed in NumPy 2
        np.float_ = np.float64
    import scdrs
    annot = annotate()
    genes = {t: gene_analysis(t, s, annot) for t, s in TRAITS.items()}

    def sha(p):
        h=hashlib.sha256()
        with p.open('rb') as stream:
            for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
        return h.hexdigest()
    inputs=[Path(__file__),Path(__file__).with_name('config.py'),INTERIM/'gse243917_annotated.h5ad']
    inputs += [TAB/f'magma_genes_{trait}.tsv' for trait in TRAITS]
    fingerprint={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    cache=TAB/'scdrs_cache.json'
    if cache.exists():
        old=json.loads(cache.read_text())
        if old['inputs']==fingerprint and all((ROOT/n).exists() and sha(ROOT/n)==v for n,v in old['outputs'].items()):
            print('Reusing byte-verified scDRS outputs from this exact code and inputs',flush=True)
            return
    adata = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    adata.obs["log_counts"] = np.log(adata.obs["total_counts"])
    cov = pd.DataFrame({"const": 1.0, "log_counts": adata.obs.log_counts, "pct_mt": adata.obs.pct_counts_mt}, index=adata.obs_names)
    scdrs.preprocess(adata, cov=cov, n_mean_bin=20, n_var_bin=20, copy=False)
    rows = []
    for t, g in genes.items():
        xmhc = (g.CHR.astype(str) == "6") & (g.STOP >= 25_000_000) & (g.START <= 34_000_000)   # extended MHC: long-range LD
        g = g[~xmhc & g.SYMBOL.isin(adata.var_names)].head(N_TOP)
        w = np.minimum(g.ZSTAT.values, 10)
        df = scdrs.score_cell(adata, gene_list=g.SYMBOL.tolist(), gene_weight=list(w), ctrl_match_key="mean_var",
                              n_ctrl=N_CTRL, weight_opt="vs", return_ctrl_raw_score=False, return_ctrl_norm_score=True,
                              random_seed=SEED, verbose=False)
        df.loc[:, ["norm_score", "pval"]].to_csv(INTERIM / "magma" / f"scdrs_cells_{t}.tsv.gz", sep="\t")
        for col in ("author_ct", "author_subtype"):
            adata.obs[col] = adata.obs[col].astype(str)
            res = scdrs.method.downstream_group_analysis(adata=adata, df_full_score=df, group_cols=[col])[col]
            res = res.reset_index().rename(columns={col: "group", "index": "group"})
            res.insert(0, "level", "cell type" if col == "author_ct" else "subtype")
            res.insert(0, "trait", t)
            rows.append(res)
        print(f"  scDRS {t}: done", flush=True)
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(TAB / "scdrs_groups.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 220)
    keep = ["level", "group", "n_cell", "n_ctrl", "assoc_mcp", "assoc_mcz", "hetero_mcp", "n_fdr_0.1"]
    keep = [k for k in keep if k in out.columns]
    for t in TRAITS:
        x = out[out.trait == t]
        print(f"\n== {t}")
        print(x[x.level == "cell type"].sort_values("assoc_mcz", ascending=False)[keep].round(4).to_string(index=False))
        print(x[x.level == "subtype"].sort_values("assoc_mcz", ascending=False)[keep].head(10).round(4).to_string(index=False))
        m6 = x[(x.level == "subtype") & (x.group == M6)]
        if len(m6):
            r = m6.iloc[0]
            print(f"  D2 M-6: assoc MC p {r.assoc_mcp:.3g} (z {r.assoc_mcz:.2f}); cells FDR < 0.1: {r['n_fdr_0.1']}")

    fig, axes = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={"width_ratios": [1, 2.2]})
    cols = {"RA_Ishigaki_EUR": "C3", "RA_FinnGen_R12": "C1", "T2D_FinnGen_R12_control": "C7"}
    for ax, level in zip(axes, ("cell type", "subtype")):
        p = out[out.level == level].pivot_table(index="group", columns="trait", values="assoc_mcz")
        if level == "subtype":
            p = p[(out[out.level == level].groupby("group").n_cell.first() >= 50).reindex(p.index).fillna(False)]
        p = p.sort_values("RA_Ishigaki_EUR")
        y = np.arange(len(p)); h = 0.8 / len(TRAITS)
        for k, t in enumerate(TRAITS):
            ax.barh(y + h * (1 - k), p[t], height=h, color=cols[t], label=t)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y, p.index, fontsize=6 if level == "subtype" else 9)
        ax.set_xlabel("scDRS group association (MC z)"); ax.set_title(f"scDRS, {level}", fontsize=10); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "scdrs.png", dpi=150)
    outputs=[TAB/'scdrs_groups.tsv']+[W/f'scdrs_cells_{trait}.tsv.gz' for trait in TRAITS]
    cache.write_text(json.dumps({'inputs':fingerprint,'outputs':{str(p.relative_to(ROOT)):sha(p) for p in outputs}},indent=2))


if __name__ == "__main__":
    from pipeline_step_lock import step_lock
    with step_lock("36"):
        main()
