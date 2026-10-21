"""GSE243917 ATAC layer: QC, TF-IDF + LSI, Harmony, and how well chromatin alone recovers RNA-defined states.

Only the peak x cell count matrix is deposited (no fragments file), so TSS enrichment and FRiP cannot be
computed; QC uses peaks and counts per nucleus. The same nuclei as the RNA analysis are used.

Input : data/raw/GSE243917/GSE243917_peaksXcells_ALL_raw.h5ad.gz, data/interim/gse243917_annotated.h5ad
Output: data/interim/gse243917_atac.h5ad
        results/tables/atac_qc_per_sample.tsv, atac_label_transfer.tsv
        results/figures/atac_qc.png, atac_umap.png
"""
import gzip
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score

from config import RAW, INTERIM, FIG, TAB, SEED

N_LSI = 30


def run_harmony(pcs, obs, key):
    import harmonypy
    ho = harmonypy.run_harmony(pcs, obs, [key], max_iter_harmony=20, random_state=SEED, verbose=False)
    z = np.asarray(ho.Z_corr)
    return z if z.shape[0] == pcs.shape[0] else z.T


def lsi_seeded(adata, n_comps):
    """muon's LSI (truncated SVD of the TF-IDF matrix, components in descending order, embeddings scaled) with a fixed
    ARPACK start and a fixed sign per component; muon calls svds without a seed, which made reruns differ slightly."""
    from scipy.sparse.linalg import svds
    u, sv, vt = svds(adata.X, k=n_comps, random_state=SEED)
    u, sv, vt = u[:, ::-1], sv[::-1], vt[::-1, :]
    sign = np.sign(vt[np.arange(vt.shape[0]), np.abs(vt).argmax(axis=1)])   # largest peak loading positive
    u, vt = u * sign, vt * sign[:, None]
    adata.obsm["X_lsi"] = (u - u.mean(axis=0)) / u.std(axis=0)
    adata.uns["lsi"] = {"stdev": sv / np.sqrt(adata.X.shape[0] - 1)}
    adata.varm["LSI"] = vt.T


def knn_transfer(adata, label, k=15):
    """Majority label among each nucleus's k nearest ATAC neighbours (self excluded)."""
    D = adata.obsp["distances"].tocsr()
    lab = adata.obs[label].astype(str).values
    pred = np.empty(len(lab), dtype=object)
    for i in range(D.shape[0]):
        nb = D.indices[D.indptr[i]:D.indptr[i + 1]]
        vals, cnt = np.unique(lab[nb], return_counts=True)
        pred[i] = vals[cnt.argmax()]
    return pred


def main():
    src = RAW / "GSE243917" / "GSE243917_peaksXcells_ALL_raw.h5ad.gz"
    dst = INTERIM / "GSE243917_peaksXcells_ALL_raw.h5ad"
    if not dst.exists():
        with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    atac = sc.read_h5ad(dst)
    rna = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad", backed="r")
    keep = ["donor", "disease", "lymphocytic_inflammation_grade", "lining_hyperplasia_grade", "prep_batch",
            "cell_type", "author_ct", "author_subtype"]
    atac = atac[atac.obs_names.isin(rna.obs_names)].copy()
    atac.obs = rna.obs.loc[atac.obs_names, keep].copy()
    rna.file.close()

    chrom = atac.var_names.str.split(":").str[0]
    coords = atac.var_names.str.split(":").str[1].str.split("-")
    atac.var["chrom"] = chrom
    atac.var["start"] = coords.str[0].astype(int)
    atac.var["end"] = coords.str[1].astype(int)
    atac.layers["counts"] = atac.X.copy()

    sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
    atac.obs = atac.obs.rename(columns={"n_genes_by_counts": "n_peaks", "total_counts": "atac_counts"})
    qc = atac.obs.groupby("donor", observed=True).agg(
        n_nuclei=("n_peaks", "size"), median_peaks=("n_peaks", "median"), median_counts=("atac_counts", "median"),
        prep_batch=("prep_batch", "first"))
    qc.to_csv(TAB / "atac_qc_per_sample.tsv", sep="\t")
    print(qc)
    atac = atac[atac.obs["n_peaks"] >= 1000].copy()     # floor: >= 1,000 accessible peaks per nucleus
    sc.pp.filter_genes(atac, min_cells=int(0.005 * atac.n_obs))
    print(f"ATAC after QC: {atac.n_obs} nuclei x {atac.n_vars} peaks")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    order = qc.index.tolist()
    for ax, col, title in zip(axes, ["n_peaks", "atac_counts"], ["Accessible peaks per nucleus", "ATAC counts in peaks"]):
        ax.boxplot([atac.obs.loc[atac.obs["donor"] == d, col] for d in order], showfliers=False)
        ax.set_xticks(range(1, len(order) + 1), [d.replace("BRI-", "") + ("*" if qc.loc[d, "prep_batch"] == "lowMT" else "")
                                                for d in order], rotation=90)
        ax.set_yscale("log"); ax.set_title(title)
    fig.suptitle("GSE243917 ATAC QC per sample (* = lowMT nuclei-prep batch)")
    fig.tight_layout(); fig.savefig(FIG / "atac_qc.png", dpi=150)

    # TF-IDF on all retained peaks, LSI on the most frequently accessible ones
    mu.atac.pp.tfidf(atac, scale_factor=1e4, log_tf=True, log_idf=True)
    frac = np.asarray((atac.layers["counts"] > 0).mean(axis=0)).ravel()
    atac.var["frac_cells"] = frac
    atac.var["highly_variable"] = frac >= np.sort(frac)[::-1][min(50000, len(frac)) - 1]
    lsi_seeded(atac, N_LSI + 1)
    depth_r = [abs(np.corrcoef(atac.obsm["X_lsi"][:, i], np.log1p(atac.obs["atac_counts"]))[0, 1]) for i in range(5)]
    print("LSI |r| with depth, comps 1-5:", np.round(depth_r, 2))
    drop = [i for i, r in enumerate(depth_r) if r > 0.7]
    lsi = np.delete(atac.obsm["X_lsi"], drop, axis=1)[:, :N_LSI]
    atac.uns["lsi_dropped_components"] = drop
    atac.obsm["X_lsi_used"] = lsi
    atac.obsm["X_lsi_harmony"] = run_harmony(lsi, atac.obs, "donor")

    sc.pp.neighbors(atac, use_rep="X_lsi_harmony", n_neighbors=15, metric="cosine", random_state=SEED)
    sc.tl.umap(atac, random_state=SEED)
    sc.tl.leiden(atac, resolution=0.6, key_added="atac_leiden", flavor="igraph", n_iterations=2, random_state=SEED)

    # How much of the RNA-defined structure does chromatin alone recover?
    rows = []
    for label in ("author_ct", "author_subtype"):
        pred = knn_transfer(atac, label)
        atac.obs[f"atac_knn_{label}"] = pred
        ok = pred == atac.obs[label].astype(str).values
        per = pd.Series(ok).groupby(atac.obs[label].astype(str).values).agg(["mean", "size"])
        for lv, r in per.iterrows():
            rows.append(dict(level=label, label=lv, n=int(r["size"]), knn_accuracy=r["mean"]))
        print(f"ATAC kNN transfer of {label}: overall accuracy {ok.mean():.3f}; "
              f"ARI(ATAC leiden, {label}) = {adjusted_rand_score(atac.obs[label], atac.obs['atac_leiden']):.3f}")
    lt = pd.DataFrame(rows).sort_values(["level", "knn_accuracy"])
    lt.round(3).to_csv(TAB / "atac_label_transfer.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    sc.pl.umap(atac, color="author_ct", ax=axes[0], show=False, frameon=False, legend_loc="on data",
               title="ATAC-only embedding, coloured by RNA cell type")
    sc.pl.umap(atac, color="donor", ax=axes[1], show=False, frameon=False, title="Sample")
    sc.pl.umap(atac, color="prep_batch", ax=axes[2], show=False, frameon=False, title="Nuclei-prep batch")
    fig.tight_layout(); fig.savefig(FIG / "atac_umap.png", dpi=150)

    atac.write_h5ad(INTERIM / "gse243917_atac.h5ad", compression="gzip")


if __name__ == "__main__":
    main()
