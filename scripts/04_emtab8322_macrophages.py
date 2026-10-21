"""E-MTAB-8322: synovial macrophage states across healthy / active RA / remission (donor-level).

The Atlas count matrix labels every barcode "macrophage" at sample level, but it also contains
fibroblasts, T/NK, B/plasma, endothelial, mural and mast cells. So this runs in two stages:
  1. all cells: QC -> Harmony -> Leiden -> broad cell type by marker score; keep myeloid only
  2. myeloid only: re-embed -> Leiden -> name subclusters by markers and a priori programs
Then donor-level features, group tests and pseudobulk DESeq2 (active vs remission,
~ sort_panel + condition) on macrophages only (DC-like clusters excluded).

Input : data/raw/E-MTAB-8322/*
Output: data/interim/emtab8322_qc.h5ad (cache), emtab8322_myeloid.h5ad
        data/processed/emtab8322_donor_features.tsv, emtab8322_macrophage_pseudobulk_counts.tsv.gz
        results/tables/emtab8322_qc_per_sample.tsv, emtab8322_broad_composition.tsv,
                       emtab8322_clusters.tsv, emtab8322_group_tests.tsv,
                       emtab8322_DE_active_vs_remission.tsv
        results/figures/emtab8322_umap.png, emtab8322_donor_features.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from config import RAW, INTERIM, PROCESSED, FIG, TAB, SEED, BROAD_MARKERS, MAC_HOMEOSTATIC, MAC_INFLAMMATORY

D = RAW / "E-MTAB-8322"
GROUP_ORDER = ["healthy", "UPA", "RA_naive_active", "RA_resistant_active", "RA_remission"]
GROUP_TICKS = ["Healthy", "UPA", "RA naive\nactive", "RA resistant\nactive", "RA\nremission"]
DC_PROGRAM = ["CLEC10A", "CD1C", "FCER1A", "HLA-DQA1", "CLEC9A", "LAMP3"]
MARKERS = {**BROAD_MARKERS, "mast": ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HPGDS"]}
NOISE_GENES = r"^(MT-|RP[SL]\d|ENSG|MALAT1|NEAT1|FTL$|FTH1$)"


def run_harmony(pcs, obs, key):
    import harmonypy
    ho = harmonypy.run_harmony(pcs, obs, [key], max_iter_harmony=20, random_state=SEED, verbose=False)
    z = np.asarray(ho.Z_corr)
    return z if z.shape[0] == pcs.shape[0] else z.T


def load_raw():
    X = scipy.io.mmread(D / "E-MTAB-8322.aggregated_counts.mtx.gz").T.tocsr().astype(np.float32)
    genes = pd.read_csv(D / "E-MTAB-8322.aggregated_counts.decorated.mtx_rows", sep="\t", header=None,
                        names=["ensembl", "symbol"])
    cells = pd.read_csv(D / "E-MTAB-8322.aggregated_counts.mtx_cols.gz", header=None)[0].values
    var = pd.DataFrame(index=genes["symbol"].fillna(genes["ensembl"]).values)
    var["ensembl"] = genes["ensembl"].values
    adata = sc.AnnData(X=X, obs=pd.DataFrame(index=cells), var=var)
    adata.var_names_make_unique()

    meta = pd.read_csv(D / "E-MTAB-8322.cell_metadata.tsv", sep="\t", low_memory=False,
                       usecols=["id", "individual", "disease", "response_to_treatment", "immunophenotype",
                                "sex", "age"]).set_index("id")
    adata = adata[adata.obs_names.isin(meta.index)].copy()
    adata.obs = adata.obs.join(meta)
    # SA149P is a second sample of SA149 (same age/sex, same UPA label): treat as one donor
    adata.obs["donor"] = adata.obs["individual"].str.replace(r"P$", "", regex=True)
    resp = adata.obs["response_to_treatment"].str.lower()
    adata.obs["group"] = np.select(
        [adata.obs["disease"] == "normal",
         adata.obs["disease"].str.startswith("undifferentiated"),
         resp.str.contains("na") & resp.str.contains("ve"),
         resp.str.contains("resistant"),
         resp.str.contains("remission")],
        GROUP_ORDER, default="other")
    adata.obs["sort_panel"] = np.where(adata.obs["immunophenotype"].str.contains("MARCO"), "extended", "basic")
    assert (adata.obs["group"] != "other").all()
    return adata


def qc_or_load():
    cache = INTERIM / "emtab8322_qc.h5ad"
    if cache.exists():
        return sc.read_h5ad(cache)
    adata = load_raw()
    n0 = adata.n_obs
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata = adata[(adata.obs["n_genes_by_counts"] >= 300) & (adata.obs["pct_counts_mt"] < 15)].copy()
    sc.pp.filter_genes(adata, min_cells=20)
    sc.pp.scrublet(adata, batch_key="individual", random_state=SEED)
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    print(f"QC kept {adata.n_obs}/{n0} cells")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.write_h5ad(cache, compression="gzip")
    return adata


def embed(adata, n_hvg, n_pcs, res):
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat", batch_key="individual")
    hv = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hv, max_value=10)
    sc.tl.pca(hv, n_comps=n_pcs, random_state=SEED)
    adata.obsm["X_pca"] = hv.obsm["X_pca"]
    adata.obsm["X_pca_harmony"] = run_harmony(hv.obsm["X_pca"], adata.obs, "individual")
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=20, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.tl.leiden(adata, resolution=res, key_added="leiden", flavor="igraph", n_iterations=2, random_state=SEED)


def top_markers(adata, n=3):
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)
    mk = sc.get.rank_genes_groups_df(adata, None)
    mk = mk[(mk["pvals_adj"] < 0.01) & (mk["logfoldchanges"] > 0.5) & ~mk["names"].str.match(NOISE_GENES)]
    return mk, mk.groupby("group").head(n).groupby("group")["names"].apply(lambda s: "+ ".join(s) + "+")


def stage1_broad(adata):
    adata.raw = adata
    embed(adata, n_hvg=2500, n_pcs=30, res=0.5)
    for ct, genes in MARKERS.items():
        sc.tl.score_genes(adata, [g for g in genes if g in adata.raw.var_names], score_name=f"score_{ct}",
                          use_raw=True, random_state=SEED)
    cols = [f"score_{ct}" for ct in MARKERS]
    cl = adata.obs.groupby("leiden", observed=True)[cols].mean()
    z = (cl - cl.mean()) / cl.std()
    broad = z.idxmax(axis=1).str.replace("score_", "")
    _, top = top_markers(adata)
    summary = pd.DataFrame({"broad_type": broad, "n_cells": adata.obs["leiden"].value_counts().reindex(broad.index),
                            "top_markers": top.reindex(broad.index)})
    print(summary.to_string())
    adata.obs["broad_type"] = adata.obs["leiden"].map(broad).astype("category")
    comp = pd.crosstab(adata.obs["individual"], adata.obs["broad_type"])
    comp.insert(0, "group", adata.obs.groupby("individual", observed=True)["group"].first())
    comp.to_csv(TAB / "emtab8322_broad_composition.tsv", sep="\t")
    return summary


def stage2_myeloid(adata):
    my = adata[adata.obs["broad_type"] == "myeloid"].copy()
    my.obs = my.obs.drop(columns=["leiden"])
    embed(my, n_hvg=2000, n_pcs=30, res=0.6)
    for name, genes in (("homeostatic", MAC_HOMEOSTATIC), ("inflammatory", MAC_INFLAMMATORY), ("dc", DC_PROGRAM)):
        sc.tl.score_genes(my, [g for g in genes if g in my.raw.var_names], score_name=f"score_{name}",
                          use_raw=True, random_state=SEED)
    mk, top = top_markers(my)
    expr = sc.get.obs_df(my, ["MERTK", "LYVE1", "FOLR2", "TREM2", "SPP1", "S100A12", "CLEC10A"], use_raw=True)
    expr["leiden"] = my.obs["leiden"].values
    progs = ["score_homeostatic", "score_inflammatory", "score_dc"]
    cl = my.obs.groupby("leiden", observed=True)[progs].mean()
    cl["pct_MERTK_pos"] = expr.groupby("leiden", observed=True)["MERTK"].apply(lambda s: 100 * (s > 0).mean())
    for g in ["LYVE1", "FOLR2", "TREM2", "SPP1", "S100A12", "CLEC10A"]:
        cl[f"mean_{g}"] = expr.groupby("leiden", observed=True)[g].mean()
    cl["n_cells"] = my.obs["leiden"].value_counts()
    cl["top_markers"] = top.reindex(cl.index)
    z = (cl[progs] - cl[progs].mean()) / cl[progs].std()
    cl["program"] = z.idxmax(axis=1).str.replace("score_", "")
    # Resident/homeostatic = clusters whose dominant program is the MERTK/LYVE1/FOLR2/SELENOP set.
    # (A MERTK-detection threshold was tried first and dropped: MERTK is detected in <=21% of cells,
    # and the threshold mislabelled the LYVE1-high SELENOP+ resident cluster as MERTK-negative.)
    cl["label"] = [f"C{i}: {t} ({p})" for i, t, p in zip(cl.index, cl["top_markers"], cl["program"])]
    cl.round(3).to_csv(TAB / "emtab8322_clusters.tsv", sep="\t")
    mk.groupby("group").head(25).to_csv(TAB / "emtab8322_cluster_markers_top25.tsv", sep="\t", index=False)
    print(cl[["n_cells", "pct_MERTK_pos", "mean_LYVE1", "mean_SPP1", "mean_S100A12", "program", "top_markers"]]
          .round(2).to_string())
    my.obs["cluster_label"] = my.obs["leiden"].map(cl["label"]).astype("category")
    my.obs["program"] = my.obs["leiden"].map(cl["program"]).astype("category")
    return my


def donor_features(my):
    mac = my[my.obs["program"] != "dc"]
    obs = mac.obs
    dlab = my.obs.groupby("donor", observed=True)[["group", "sort_panel", "sex", "age"]].first()
    feats = pd.DataFrame(index=dlab.index)
    feats["n_myeloid"] = my.obs.groupby("donor", observed=True).size()
    feats["n_macrophages"] = obs.groupby("donor", observed=True).size()
    feats["frac_DC_of_myeloid"] = 1 - feats["n_macrophages"] / feats["n_myeloid"]
    is_hom = (obs["program"] == "homeostatic").groupby(obs["donor"], observed=True)
    is_inf = (obs["program"] == "inflammatory").groupby(obs["donor"], observed=True)
    feats["mac_frac_homeostatic"] = is_hom.mean()
    feats["mac_frac_inflammatory"] = is_inf.mean()
    feats["mac_homeostasis_index"] = np.log2((is_hom.sum() + 1) / (is_inf.sum() + 1))
    for s in ("score_homeostatic", "score_inflammatory"):
        feats[f"mean_{s}"] = obs.groupby("donor", observed=True)[s].mean()
    clfrac = pd.crosstab(obs["donor"], obs["leiden"], normalize="index").add_prefix("mac_frac_C")
    out = dlab.join(feats).join(clfrac)
    out.round(4).to_csv(PROCESSED / "emtab8322_donor_features.tsv", sep="\t")
    return out, mac


def group_tests(out):
    rows = []
    test_cols = [c for c in out.columns if c.startswith(("mac_", "mean_score", "frac_DC"))]
    act = out["group"].isin(["RA_naive_active", "RA_resistant_active"])
    for c in test_cols:
        groups = [out.loc[out["group"] == g, c].values for g in GROUP_ORDER]
        a, r, hlt = out.loc[act, c], out.loc[out["group"] == "RA_remission", c], out.loc[out["group"] == "healthy", c]
        rows.append(dict(feature=c,
                         **{f"median_{g}": out.loc[out["group"] == g, c].median() for g in GROUP_ORDER},
                         p_kruskal=kruskal(*groups).pvalue,
                         p_active_vs_remission=mannwhitneyu(a, r).pvalue,
                         p_active_vs_healthy=mannwhitneyu(a, hlt).pvalue))
    res = pd.DataFrame(rows)
    for col in ("p_kruskal", "p_active_vs_remission", "p_active_vs_healthy"):
        res[col.replace("p_", "fdr_")] = multipletests(res[col], method="fdr_bh")[1]
    res = res.sort_values("p_active_vs_remission")
    res.round(4).to_csv(TAB / "emtab8322_group_tests.tsv", sep="\t", index=False)
    print(res[["feature", "median_healthy", "median_RA_naive_active", "median_RA_resistant_active",
               "median_RA_remission", "p_active_vs_remission", "fdr_active_vs_remission",
               "p_active_vs_healthy", "fdr_active_vs_healthy"]].round(3).to_string(index=False))
    return res


def pseudobulk_de(mac, out):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    C = mac.layers["counts"]
    donors = out.index.tolist()
    pb = np.vstack([np.asarray(C[(mac.obs["donor"] == d).values].sum(axis=0)).ravel() for d in donors])
    pb = pd.DataFrame(pb.astype(int), index=donors, columns=mac.var_names)
    pb.T.to_csv(PROCESSED / "emtab8322_macrophage_pseudobulk_counts.tsv.gz", sep="\t")
    md = out.loc[out["group"].isin(["RA_naive_active", "RA_resistant_active", "RA_remission"]), ["group", "sort_panel"]].copy()
    md["condition"] = np.where(md["group"] == "RA_remission", "remission", "active")
    counts = pb.loc[md.index]
    counts = counts.loc[:, (counts >= 10).sum(axis=0) >= 5]
    inf = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(counts=counts, metadata=md, design="~sort_panel + condition", refit_cooks=True,
                       inference=inf, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "active", "remission"], inference=inf, quiet=True)
    st.summary()
    de = st.results_df.sort_values("pvalue")
    de.round(5).to_csv(TAB / "emtab8322_DE_active_vs_remission.tsv", sep="\t")
    print(f"DE active vs remission, macrophages only (n={int((md.condition == 'active').sum())} vs "
          f"{int((md.condition == 'remission').sum())}): {(de.padj < 0.05).sum()} genes padj<0.05 "
          f"({((de.padj < 0.05) & (de.log2FoldChange > 0)).sum()} up in active)")
    print(de.head(30)[["baseMean", "log2FoldChange", "pvalue", "padj"]].round(4).to_string())


def figures(adata, my, out, res):
    fig, axes = plt.subplots(1, 4, figsize=(28, 6))
    sc.pl.umap(adata, color="broad_type", ax=axes[0], show=False, frameon=False, legend_loc="on data",
               title="Stage 1: all barcodes, broad type")
    sc.pl.umap(my, color="leiden", ax=axes[1], show=False, frameon=False, legend_loc="on data",
               title="Stage 2: myeloid subclusters")
    sc.pl.umap(my, color="score_homeostatic", ax=axes[2], show=False, frameon=False, cmap="viridis",
               title="Homeostatic program (MERTK/LYVE1/FOLR2...)")
    sc.pl.umap(my, color="score_inflammatory", ax=axes[3], show=False, frameon=False, cmap="magma",
               title="Inflammatory program (S100A8/9/12, SPP1, IL1B...)")
    fig.tight_layout(); fig.savefig(FIG / "emtab8322_umap.png", dpi=150)

    show = ["mac_homeostasis_index", "mac_frac_inflammatory", "mean_score_homeostatic", "mean_score_inflammatory"]
    fig, axes = plt.subplots(1, len(show), figsize=(5 * len(show), 4.5))
    rng = np.random.default_rng(SEED)
    for ax, c in zip(axes, show):
        for i, g in enumerate(GROUP_ORDER):
            v = out.loc[out["group"] == g, c]
            ax.boxplot(v, positions=[i], widths=0.5, showfliers=False)
            ax.scatter(i + rng.uniform(-.12, .12, len(v)), v, s=25, zorder=3,
                       c=["C1" if p == "extended" else "C0" for p in out.loc[v.index, "sort_panel"]])
        r = res.set_index("feature").loc[c]
        ax.set_xticks(range(len(GROUP_ORDER)), GROUP_TICKS, fontsize=8)
        ax.set_title(f"{c}\nactive vs healthy p={r.p_active_vs_healthy:.3f}; vs remission p={r.p_active_vs_remission:.3f}",
                     fontsize=9)
    fig.suptitle("E-MTAB-8322 donor-level macrophage features (blue = basic sort panel, orange = extended)")
    fig.tight_layout(); fig.savefig(FIG / "emtab8322_donor_features.png", dpi=150)


def main():
    adata = qc_or_load()
    adata.obs.groupby("individual", observed=True).agg(
        group=("group", "first"), n_cells=("group", "size"),
        median_genes=("n_genes_by_counts", "median"), median_pct_mt=("pct_counts_mt", "median")
    ).to_csv(TAB / "emtab8322_qc_per_sample.tsv", sep="\t")
    stage1_broad(adata)
    my = stage2_myeloid(adata)
    out, mac = donor_features(my)
    res = group_tests(out)
    pseudobulk_de(mac, out)
    figures(adata, my, out, res)
    my.write_h5ad(INTERIM / "emtab8322_myeloid.h5ad", compression="gzip")


if __name__ == "__main__":
    main()
