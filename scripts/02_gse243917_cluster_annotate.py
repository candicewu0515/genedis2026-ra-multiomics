"""GSE243917: normalise, integrate across donors (Harmony), cluster, annotate broad cell types.

Our annotation is marker-driven and blind to the authors' labels; the authors' labels are used
only afterwards to measure agreement (ARI + confusion matrix).

Input : data/interim/gse243917_qc.h5ad
Output: data/interim/gse243917_annotated.h5ad
        results/tables/gse243917_cluster_annotation.tsv, gse243917_annotation_vs_author.tsv
        results/figures/gse243917_umap.png, gse243917_marker_dotplot.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

from config import INTERIM, FIG, TAB, SEED, BROAD_MARKERS

sc.settings.figdir = FIG


def run_harmony(pcs, obs, key):
    """Call harmonypy directly; its Z_corr orientation changed across versions, so normalise it here."""
    import harmonypy
    ho = harmonypy.run_harmony(pcs, obs, [key], max_iter_harmony=20, random_state=SEED, verbose=False)
    z = np.asarray(ho.Z_corr)
    return z if z.shape[0] == pcs.shape[0] else z.T


def main():
    adata = sc.read_h5ad(INTERIM / "gse243917_qc.h5ad")
    adata = adata[adata.obs["qc_pass"] & ~adata.obs["predicted_doublet"]].copy()
    sc.pp.filter_genes(adata, min_cells=10)

    # Low-mito vs high-mito samples look like two nuclei-prep batches; keep it as a covariate
    mt_by_donor = adata.obs.groupby("donor", observed=True)["pct_counts_mt"].median()
    adata.obs["prep_batch"] = adata.obs["donor"].map(lambda d: "lowMT" if mt_by_donor[d] < 3 else "highMT")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # seurat_v3 loess fails on the 137-nucleus sample; the log-based seurat flavour is stable per batch
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", batch_key="donor")
    hv = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hv, max_value=10)
    sc.tl.pca(hv, n_comps=40, random_state=SEED)
    adata.obsm["X_pca"] = hv.obsm["X_pca"]
    adata.obsm["X_pca_harmony"] = run_harmony(hv.obsm["X_pca"], adata.obs, "donor")

    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=20, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    for res in (0.3, 0.6, 1.0):
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}", flavor="igraph",
                     n_iterations=2, random_state=SEED)
    adata.obs["leiden"] = adata.obs["leiden_0.6"]

    # Score every broad marker set; each cluster takes the set with the highest mean score
    for ct, genes in BROAD_MARKERS.items():
        genes = [g for g in genes if g in adata.raw.var_names]
        sc.tl.score_genes(adata, genes, score_name=f"score_{ct}", use_raw=True, random_state=SEED)
    score_cols = [f"score_{ct}" for ct in BROAD_MARKERS]
    cl_scores = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    z = (cl_scores - cl_scores.mean()) / cl_scores.std()
    best = z.idxmax(axis=1).str.replace("score_", "")
    margin = z.apply(lambda r: r.nlargest(2).iloc[0] - r.nlargest(2).iloc[1], axis=1)
    ann = pd.DataFrame({"cell_type": best, "z_margin": margin.round(2),
                        "n_cells": adata.obs["leiden"].value_counts().reindex(best.index)})
    ann.to_csv(TAB / "gse243917_cluster_annotation.tsv", sep="\t")
    adata.obs["cell_type"] = adata.obs["leiden"].map(best).astype("category")

    ari = adjusted_rand_score(adata.obs["author_ct"], adata.obs["cell_type"])
    conf = pd.crosstab(adata.obs["cell_type"], adata.obs["author_ct"])
    conf.to_csv(TAB / "gse243917_annotation_vs_author.tsv", sep="\t")
    agree = (adata.obs["cell_type"].astype(str) == adata.obs["author_ct"].astype(str)).mean()
    print(ann)
    print(conf)
    print(f"ARI vs author labels = {ari:.3f}; per-cell agreement = {agree:.1%}")

    fig, axes = plt.subplots(1, 4, figsize=(24, 5.5))
    for ax, col, title in zip(axes, ["cell_type", "author_ct", "donor", "prep_batch"],
                              ["Our annotation", "Author annotation", "Sample", "Nuclei-prep batch (inferred)"]):
        sc.pl.umap(adata, color=col, ax=ax, show=False, frameon=False, title=title,
                   legend_loc="on data" if col in ("cell_type", "author_ct") else "right margin",
                   legend_fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "gse243917_umap.png", dpi=150)

    order = list(BROAD_MARKERS)
    adata.obs["cell_type"] = adata.obs["cell_type"].cat.reorder_categories(
        [c for c in order if c in adata.obs["cell_type"].cat.categories])
    dp = sc.pl.dotplot(adata, {k: [g for g in v if g in adata.raw.var_names] for k, v in BROAD_MARKERS.items()},
                       groupby="cell_type", use_raw=True, standard_scale="var", show=False, return_fig=True)
    dp.savefig(FIG / "gse243917_marker_dotplot.png", dpi=150, bbox_inches="tight")

    adata.uns["ari_vs_author"] = float(ari)
    adata.write_h5ad(INTERIM / "gse243917_annotated.h5ad", compression="gzip")


if __name__ == "__main__":
    main()
