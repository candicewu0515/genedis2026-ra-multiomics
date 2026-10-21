"""Robustness and ablation of the key results (reuses outputs of earlier steps; no new data).

R1  active-RA macrophage signature (E-MTAB-8322 DESeq2, 12 active vs 7 remission): leave-one-donor-out refits.
R2  active RA vs healthy macrophage features: leave-one-donor-out Mann-Whitney p.
R3  donor-level RNA vs ATAC signature agreement (GSE243917, n = 10): leave-one-out rho.
R4  signature transfer to GSE243917 vs inflammation grade (n = 10): leave-one-out rho.
R5  RA GWAS enrichment in T-cell-specific peaks: number of specific peaks 2,000 / 4,000 / 8,000.
R6  GSEA cross-contrast agreement: ranking by Wald statistic vs log2 fold change.

Output: results/tables/robustness_summary.tsv, robustness_R1_loo_genes.tsv; results/figures/robustness.png
"""
import importlib.util

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

from config import PROCESSED, FIG, TAB, SEED, ROOT

rows = []


def add(analysis, metric, full, loo_vals, sig_rule=None, note=""):
    v = np.asarray(loo_vals, dtype=float)
    rows.append(dict(analysis=analysis, metric=metric, full_data=full, min=np.nanmin(v), max=np.nanmax(v),
                     n_perturbations=len(v),
                     fraction_meeting_rule=(np.mean(sig_rule(v)) if sig_rule else np.nan), note=note))


def load_step(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def r1_signature():
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    pb = pd.read_csv(PROCESSED / "emtab8322_macrophage_pseudobulk_counts.tsv.gz", sep="\t", index_col=0).T
    f = pd.read_csv(PROCESSED / "emtab8322_donor_features.tsv", sep="\t", index_col=0)
    md = f.loc[f["group"].isin(["RA_naive_active", "RA_resistant_active", "RA_remission"]), ["group", "sort_panel"]].copy()
    md["condition"] = np.where(md["group"] == "RA_remission", "remission", "active")
    full = pd.read_csv(TAB / "emtab8322_DE_active_vs_remission.tsv", sep="\t", index_col=0)
    sig = full[(full.padj < 0.05) & (full.log2FoldChange > 0.5) & ~full.index.str.match(r"^(ENSG|MT-)")].index.tolist()
    inf = DefaultInference(n_cpus=8)
    keep_frac, sign_ok, rho_stat, detail = [], [], [], []
    for d in md.index:
        m = md.drop(d)
        c = pb.loc[m.index]
        c = c.loc[:, (c >= 10).sum(axis=0) >= 5]
        dds = DeseqDataSet(counts=c, metadata=m, design="~sort_panel + condition", refit_cooks=True, inference=inf, quiet=True)
        dds.deseq2()
        st = DeseqStats(dds, contrast=["condition", "active", "remission"], inference=inf, quiet=True)
        st.summary()
        r = st.results_df
        g = [x for x in sig if x in r.index]
        keep_frac.append(np.mean(r.loc[g, "padj"] < 0.05))
        sign_ok.append(np.mean(r.loc[g, "log2FoldChange"] > 0))
        common = r.index.intersection(full.index)
        rho_stat.append(spearmanr(r.loc[common, "stat"], full.loc[common, "stat"], nan_policy="omit")[0])
        detail.append(pd.Series(r.loc[g, "padj"].values, index=g, name=f"drop_{d}"))
    pd.DataFrame(detail).T.round(5).to_csv(TAB / "robustness_R1_loo_genes.tsv", sep="\t")
    add("R1 signature (E-MTAB-8322)", "fraction of the 15 FDR genes still FDR < 0.05", 1.0, keep_frac,
        lambda v: v >= 0.5, "rule: at least half retained")
    add("R1 signature (E-MTAB-8322)", "fraction of the 15 genes still up in active", 1.0, sign_ok, lambda v: v == 1.0,
        "rule: all directions retained")
    add("R1 signature (E-MTAB-8322)", "Spearman of genome-wide Wald statistics vs full data", 1.0, rho_stat,
        lambda v: v >= 0.8)
    return keep_frac


def r2_active_vs_healthy():
    f = pd.read_csv(PROCESSED / "emtab8322_donor_features.tsv", sep="\t", index_col=0)
    act = f["group"].isin(["RA_naive_active", "RA_resistant_active"])
    sub = f[act | (f["group"] == "healthy")]
    for col in ("mac_frac_inflammatory", "mean_score_homeostatic"):
        full_p = mannwhitneyu(sub.loc[act[sub.index], col], sub.loc[~act[sub.index], col]).pvalue
        ps = []
        for d in sub.index:
            s = sub.drop(d)
            a = act[s.index]
            ps.append(mannwhitneyu(s.loc[a, col], s.loc[~a, col]).pvalue)
        add("R2 active RA vs healthy (E-MTAB-8322)", f"Mann-Whitney p, {col}", full_p, ps, lambda v: v < 0.05)


def loo_rho(x, y):
    return [spearmanr(np.delete(x, i), np.delete(y, i))[0] for i in range(len(x))]


def r3_r4():
    d = pd.read_csv(TAB / "multiome_donor_signature.tsv", sep="\t", index_col=0)
    ra = d[d["disease"] == "RA"]
    x, y = ra["rna_signature_score"].values, ra["atac_gene_activity_signature_score"].values
    add("R3 RNA vs ATAC signature (GSE243917)", "donor-level Spearman rho", spearmanr(x, y)[0], loo_rho(x, y),
        lambda v: v >= 0.69, "rule: above the 95th percentile of random gene sets (0.69)")
    s = pd.read_csv(TAB / "cross_dataset_signature_scores.tsv", sep="\t", index_col=0)
    s = s[(s["disease"] == "RA")].dropna(subset=["score_primary_FDR"])
    x, y = s["score_primary_FDR"].values, s["lymphocytic_inflammation_grade"].values
    add("R4 signature transfer (GSE243917)", "Spearman rho vs lymphocytic inflammation grade", spearmanr(x, y)[0],
        loo_rho(x, y), lambda v: v > 0, "rule: direction positive")


def r5_gwas():
    g = load_step("10_gwas_peak_enrichment")
    import scanpy as sc
    rng = np.random.default_rng(SEED)
    atac = sc.read_h5ad(ROOT / "data" / "interim" / "gse243917_atac.h5ad")
    peaks = atac.var[["chrom", "start", "end", "frac_cells"]].copy()
    mhc = g.MHC
    ok = ~((peaks["chrom"] == mhc[0]) & (peaks["end"] >= mhc[1]) & (peaks["start"] <= mhc[2])).values
    peaks = peaks[ok]
    tss = g.tss_table()
    mid = ((peaks["start"] + peaks["end"]) // 2).values
    near = np.zeros(len(peaks), bool)
    for chrom, t in tss.groupby("chrom"):
        idx = np.where(peaks["chrom"].values == chrom)[0]
        if len(idx):
            tt = np.sort(t["tss"].values)
            j = np.clip(np.searchsorted(tt, mid[idx]), 1, len(tt) - 1)
            near[idx] = np.minimum(np.abs(mid[idx] - tt[j - 1]), np.abs(mid[idx] - tt[j])) <= 2000
    strata = pd.qcut(peaks["frac_cells"].rank(method="first"), 10, labels=False).values * 2 + near.astype(int)
    pools = {s: np.where(strata == s)[0] for s in np.unique(strata)}
    loci = g.merge_loci(g.lead_snps("RA"), 25_000)
    lid = g.peak_locus_ids(peaks, loci)
    out = {}
    for n_top in (2000, 4000, 8000):
        g.N_TOP = n_top
        sets = {k: v for k, v in g.specific_sets(atac, ok).items() if k == ("cell type", "Tcell")}
        idx = sets[("cell type", "Tcell")]
        obs = len(np.unique(lid[idx][lid[idx] >= 0]))
        need = pd.Series(strata[idx]).value_counts()
        null = []
        for _ in range(1000):
            pick = np.concatenate([rng.choice(pools[s], n, replace=False) for s, n in need.items()])
            h = lid[pick]
            null.append(len(np.unique(h[h >= 0])))
        null = np.array(null)
        out[n_top] = (obs / null.mean(), (1 + np.sum(null >= obs)) / 1001)
        print(f"  R5 N_TOP={n_top}: T-cell enrichment {out[n_top][0]:.2f}, p = {out[n_top][1]:.4f}")
    add("R5 RA GWAS enrichment, T-cell peaks (25 kb)", "enrichment (loci hit / matched null)", out[4000][0],
        [out[k][0] for k in (2000, 8000)], lambda v: v > 1, "perturbations: 2,000 and 8,000 specific peaks")
    add("R5 RA GWAS enrichment, T-cell peaks (25 kb)", "empirical p", out[4000][1],
        [out[k][1] for k in (2000, 8000)], lambda v: v < 0.05)


def r6_gsea():
    import gseapy as gp
    g12 = load_step("12_pathway_gsea")
    sets = g12.load_gmt("hallmark", "h.all")
    a = pd.read_csv(TAB / "emtab8322_DE_active_vs_remission.tsv", sep="\t", index_col=0).dropna(subset=["stat"])
    a = a[~a.index.str.match(r"^(ENSG|MT-)")]
    b = pd.read_csv(TAB / "gse198520_DE_post_vs_pre_good.tsv", sep="\t", index_col=0).dropna(subset=["stat"])
    res = {}
    for metric in ("stat", "log2FoldChange"):
        nes = {}
        for tag, de in (("A", a), ("B", b)):
            rnk = de[metric].sort_values(ascending=False).reset_index()
            rnk.columns = ["gene", "score"]
            r = gp.prerank(rnk=rnk, gene_sets=sets, min_size=15, max_size=500, permutation_num=1000, seed=SEED,
                           threads=8, outdir=None, verbose=False).res2d
            nes[tag] = r.set_index("Term")["NES"].astype(float)
        common = nes["A"].index.intersection(nes["B"].index)
        res[metric] = spearmanr(nes["A"][common], nes["B"][common])[0]
        print(f"  R6 ranking by {metric}: Hallmark NES rho(A, B) = {res[metric]:.2f}")
    add("R6 GSEA cross-contrast (Hallmark)", "Spearman rho of NES, A vs B", res["stat"], [res["log2FoldChange"]],
        lambda v: v < -0.3, "perturbation: rank by log2 fold change instead of Wald statistic")


def main():
    print("R1 ..."); keep = r1_signature()
    print("R2 ..."); r2_active_vs_healthy()
    print("R3/R4 ..."); r3_r4()
    print("R5 ..."); r5_gwas()
    print("R6 ..."); r6_gsea()
    summ = pd.DataFrame(rows)
    summ.round(4).to_csv(TAB / "robustness_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
    print(summ.round(3).to_string(index=False))

    fig, ax = plt.subplots(figsize=(10, 0.45 * len(summ) + 1.5))
    for i, r in summ.iterrows():
        ax.plot([r["min"], r["max"]], [i, i], color="C0", lw=3)
        ax.scatter(r["full_data"], i, color="C3", zorder=3)
    ax.set_yticks(range(len(summ)), [f"{a.split(' (')[0]}: {m}"[:80] for a, m in zip(summ.analysis, summ.metric)], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("value (red = full data; blue = range under leave-one-out / parameter perturbation)")
    ax.set_title("Robustness of key results", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "robustness.png", dpi=150)


if __name__ == "__main__":
    main()
