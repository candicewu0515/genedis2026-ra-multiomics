"""GSE243917 RNA + ATAC joint embedding with weighted nearest neighbours (WNN, muon).

Each modality keeps its own Harmony-corrected embedding (RNA PCA, ATAC LSI); WNN learns a per-nucleus
weight for each modality. We report (i) how well RNA-only, ATAC-only and joint graphs recover author
subtypes, and (ii) the per-subtype RNA weight, which shows the cell states defined more by chromatin.

Input : data/interim/gse243917_annotated.h5ad, data/interim/gse243917_atac.h5ad
Output: data/interim/gse243917_multiome.h5mu
        results/tables/wnn_modality_weights.tsv, wnn_recovery.tsv
        results/figures/wnn_umap.png, wnn_modality_weights.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

from config import INTERIM, PROCESSED, FIG, TAB, SEED


def knn_accuracy(conn, labels):
    """Fraction of each nucleus's graph neighbours (weighted) sharing its label, averaged per nucleus."""
    conn = conn.tocsr()
    lab = np.asarray(labels)
    acc = np.zeros(conn.shape[0])
    for i in range(conn.shape[0]):
        idx = conn.indices[conn.indptr[i]:conn.indptr[i + 1]]
        w = conn.data[conn.indptr[i]:conn.indptr[i + 1]]
        acc[i] = np.average(lab[idx] == lab[i], weights=w) if len(idx) else np.nan
    return acc


def main():
    rna = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    atac = sc.read_h5ad(INTERIM / "gse243917_atac.h5ad")
    common = rna.obs_names.intersection(atac.obs_names)
    rna, atac = rna[common].copy(), atac[common].copy()

    sc.pp.neighbors(rna, use_rep="X_pca_harmony", n_neighbors=20, random_state=SEED)
    sc.pp.neighbors(atac, use_rep="X_lsi_harmony", n_neighbors=20, metric="cosine", random_state=SEED)
    mdata = mu.MuData({"rna": rna, "atac": atac})
    mu.pp.neighbors(mdata, key_added="wnn", n_neighbors=20, random_state=SEED)
    mu.tl.umap(mdata, neighbors_key="wnn", random_state=SEED)
    # muon's leiden clusters per modality; cluster the joint WNN graph directly instead
    tmp = sc.AnnData(obs=pd.DataFrame(index=mdata.obs_names))
    sc.tl.leiden(tmp, adjacency=mdata.obsp["wnn_connectivities"], resolution=1.0, flavor="igraph",
                 n_iterations=2, random_state=SEED)
    mdata.obs["wnn_leiden"] = tmp.obs["leiden"].values

    obs = mdata.obs
    sub = obs["rna:author_subtype"].astype(str)
    ct = obs["rna:author_ct"].astype(str)
    w_rna = obs["rna:mod_weight"].values

    # Recovery of author labels by each graph
    graphs = {"RNA": rna.obsp["connectivities"], "ATAC": atac.obsp["connectivities"],
              "WNN": mdata.obsp["wnn_connectivities"]}
    rows = []
    for name, g in graphs.items():
        for level, lab in (("cell type", ct), ("subtype", sub)):
            acc = knn_accuracy(g, lab.values)
            rows.append(dict(graph=name, level=level, mean_neighbour_label_agreement=np.nanmean(acc)))
    rows.append(dict(graph="WNN leiden", level="subtype",
                     mean_neighbour_label_agreement=adjusted_rand_score(sub, obs["wnn_leiden"]), note="ARI"))
    rows.append(dict(graph="RNA leiden (step 02)", level="subtype",
                     mean_neighbour_label_agreement=adjusted_rand_score(sub, obs["rna:leiden"]), note="ARI"))
    rec = pd.DataFrame(rows)
    rec.round(3).to_csv(TAB / "wnn_recovery.tsv", sep="\t", index=False)
    print(rec.round(3).to_string(index=False))

    # Per-subtype modality weight
    wt = pd.DataFrame({"cell_type": ct, "subtype": sub, "rna_weight": w_rna})
    per = wt.groupby(["cell_type", "subtype"]).agg(n=("rna_weight", "size"), median_rna_weight=("rna_weight", "median")) \
        .reset_index().sort_values("median_rna_weight")
    per = per[per["n"] >= 30]
    per.round(3).to_csv(TAB / "wnn_modality_weights.tsv", sep="\t", index=False)
    print("Lowest RNA weight (chromatin-defined) subtypes:\n", per.head(8).round(3).to_string(index=False))
    print("Highest RNA weight subtypes:\n", per.tail(5).round(3).to_string(index=False))
    print("Median RNA weight by broad type:\n", wt.groupby("cell_type")["rna_weight"].median().round(3).to_string())

    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    mdata.obs["author_ct"] = ct.values
    mdata.obs["donor"] = obs["rna:donor"].astype(str).values
    mu.pl.umap(mdata, color="author_ct", ax=axes[0], show=False, frameon=False, legend_loc="on data",
               title="WNN (RNA+ATAC) embedding")
    mu.pl.umap(mdata, color="rna:mod_weight", ax=axes[1], show=False, frameon=False, cmap="coolwarm",
               title="Per-nucleus RNA weight (1 - ATAC weight)")
    mu.pl.umap(mdata, color="donor", ax=axes[2], show=False, frameon=False, title="Sample")
    fig.tight_layout(); fig.savefig(FIG / "wnn_umap.png", dpi=150)

    fig, ax = plt.subplots(figsize=(7, max(4, 0.22 * len(per))))
    colors = {c: f"C{i}" for i, c in enumerate(sorted(per["cell_type"].unique()))}
    ax.barh(range(len(per)), per["median_rna_weight"], color=[colors[c] for c in per["cell_type"]])
    ax.set_yticks(range(len(per)), per["subtype"], fontsize=7)
    ax.axvline(0.5, color="k", lw=0.8, ls="--")
    ax.set_xlabel("Median WNN RNA weight (<0.5 = chromatin contributes more)")
    for c, col in colors.items():
        ax.barh([], [], color=col, label=c)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "wnn_modality_weights.png", dpi=150)

    mdata.write(INTERIM / "gse243917_multiome.h5mu")
    plot = rna.obs[['author_ct', 'author_subtype', 'donor', 'disease']].reindex(mdata.obs_names).copy()
    plot['umap1'] = mdata.obsm['X_umap'][:, 0]; plot['umap2'] = mdata.obsm['X_umap'][:, 1]
    plot.to_csv(PROCESSED / 'multiome_plot_cells.tsv.gz', sep='\t', index=True,
                compression={'method': 'gzip', 'mtime': 0})


if __name__ == "__main__":
    main()
