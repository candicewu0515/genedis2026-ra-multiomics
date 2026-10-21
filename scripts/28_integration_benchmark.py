"""Integration baselines: is the WNN joint graph (step 07) better than single-modality and other integration methods?

Same 31,520 nuclei and the same input features for every method; every embedding is Harmony-corrected by sample
except the uncorrected RNA reference.
  RNA PCA (uncorrected)   3,000 HVG, 40 PCs (step 02, before Harmony)
  RNA PCA + Harmony       step 02
  ATAC LSI + Harmony      TF-IDF on 50,000 peaks, 30 LSI components (step 06)
  Concatenation           z-scored RNA Harmony PCs and ATAC Harmony LSI side by side (simple early-integration baseline)
  WNN                     step 07 (muon weighted nearest neighbours)
  MOFA+                   20 factors on the 3,000 HVG (log-normalised) + 20,000 most accessible of the 50,000 peaks
                          (TF-IDF), views scaled to equal variance; factors then Harmony-corrected by sample
Each method gives a 20-nearest-neighbour graph. Metrics (scIB-style, computed on the graphs so WNN is comparable):
  biology    Leiden (resolutions 0.2-2.0) best NMI and ARI against three label sets: author cell types and author
             subtypes (both defined on RNA, so they favour RNA) and our ATAC-only Leiden clusters (step 06, favour ATAC);
             neighbour label agreement for subtypes; recovery of the rare M-6 STAT1+ CXCL10+ state (65 nuclei):
             neighbour purity and best cluster F1
  batch      neighbour sample entropy relative to the entropy expected under perfect mixing within the nucleus's
             cell type; graph connectivity of each cell type (largest connected component share)
  modality   label-free: Jaccard overlap of each nucleus's neighbours with its RNA-only and ATAC-only neighbours
Overall = 0.6 x mean scaled biology + 0.4 x mean scaled batch (scIB weights), min-max scaled across methods.

Output: data/interim/gse243917_mofa.hdf5, gse243917_mofa_factors.tsv.gz
        results/tables/integration_benchmark.tsv, mofa_factors_top_features.tsv, mofa_variance_explained.tsv
        results/figures/integration_benchmark.png
"""
import gzip
import re

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

from config import RAW, INTERIM, FIG, TAB, SEED

K = 20
RESOLUTIONS = (0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0)
N_FACTORS = 20
N_PEAKS_MOFA = 20_000
M6 = "M-6: STAT1+ CXCL10+"
IFN = ["CXCL9", "CXCL10", "CXCL11", "STAT1", "GBP1", "GBP5", "IDO1", "IRF1", "CD274", "APOL6"]


def run_harmony(x, obs, key):
    import harmonypy
    ho = harmonypy.run_harmony(x, obs, [key], max_iter_harmony=20, random_state=SEED, verbose=False)
    z = np.asarray(ho.Z_corr)
    return z if z.shape[0] == x.shape[0] else z.T


def tss_table():
    rows = []
    with gzip.open(RAW / "annotation" / "Homo_sapiens.GRCh38.98.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene" or 'gene_biotype "protein_coding"' not in f[8]:
                continue
            rows.append(("chr" + f[0], int(f[3]) if f[6] == "+" else int(f[4]), re.search(r'gene_name "([^"]+)"', f[8]).group(1)))
    return pd.DataFrame(rows, columns=["chrom", "tss", "gene"])


def nearest_gene(peak_names, tss):
    out = []
    for p in peak_names:
        c, se = p.split(":"); s, e = map(int, se.split("-")); m = (s + e) // 2
        t = tss[tss.chrom == c]
        out.append(t.gene.values[np.argmin(np.abs(t.tss.values - m))] if len(t) else "")
    return out


def train_mofa(rna, atac):
    """MOFA+ on log-normalised HVGs and TF-IDF peaks; cached because training takes ~1 h on CPU."""
    out = INTERIM / "gse243917_mofa.hdf5"
    if not out.exists():
        import muon as mu
        r = rna[:, rna.var["highly_variable"].values].copy()
        a = atac[:, atac.var["highly_variable"].values].copy()
        a.X = a.layers["counts"].copy()
        mu.atac.pp.tfidf(a, scale_factor=1e4, log_tf=True, log_idf=True)
        top = np.argsort(-a.var["frac_cells"].values)[:N_PEAKS_MOFA]
        a = a[:, np.sort(top)].copy()
        md = mu.MuData({"rna": r, "atac": a})
        mu.tl.mofa(md, n_factors=N_FACTORS, scale_views=True, convergence_mode="fast", use_float32=True,
                   seed=SEED, outfile=str(out), save_data=False, verbose=False)
    with h5py.File(out, "r") as f:
        Z = np.concatenate([f["expectations/Z"][g][:].T for g in f["expectations/Z"].keys()])
        samples = np.concatenate([f["samples"][g][:].astype(str) for g in f["samples"].keys()])
        W = {v: f["expectations/W"][v][:].T for v in f["expectations/W"].keys()}
        feats = {v: f["features"][v][:].astype(str) for v in f["features"].keys()}
        views = f["views/views"][:].astype(str)
        g0 = list(f["variance_explained/r2_per_factor"].keys())[0]
        r2 = pd.DataFrame(np.asarray(f["variance_explained/r2_per_factor"][g0][:]).T, columns=views)   # factors x views
    Z = pd.DataFrame(Z, index=samples, columns=[f"Factor{i + 1}" for i in range(Z.shape[1])])
    return Z, W, feats, r2


def graph_from(rep, metric="euclidean"):
    tmp = sc.AnnData(obs=pd.DataFrame(index=np.arange(rep.shape[0]).astype(str)))
    tmp.obsm["X"] = np.asarray(rep, dtype=np.float32)
    sc.pp.neighbors(tmp, use_rep="X", n_neighbors=K, metric=metric, random_state=SEED)
    return tmp.obsp["distances"].tocsr(), tmp.obsp["connectivities"].tocsr()


def nbr_lists(D):
    D = D.tocsr()
    return [D.indices[D.indptr[i]:D.indptr[i + 1]] for i in range(D.shape[0])]


def entropy(counts):
    p = counts[counts > 0] / counts.sum()
    return float(-(p * np.log(p)).sum())


def evaluate(name, D, C, truths, ct, donor, nb_rna, nb_atac, rng):
    nb = nbr_lists(D)
    row = dict(method=name)
    tmp = sc.AnnData(obs=pd.DataFrame(index=np.arange(D.shape[0]).astype(str)))
    clus = {}
    for r in RESOLUTIONS:
        sc.tl.leiden(tmp, adjacency=C, resolution=r, flavor="igraph", n_iterations=2, random_state=SEED, key_added="l")
        clus[r] = tmp.obs["l"].astype(str).values
    for tname, lab in truths.items():
        nmi = {r: normalized_mutual_info_score(lab, c) for r, c in clus.items()}
        best = max(nmi, key=nmi.get)
        row[f"NMI_{tname}"] = nmi[best]
        row[f"ARI_{tname}"] = adjusted_rand_score(lab, clus[best])
    sub = truths["author_subtype"]
    row["knn_agreement_subtype"] = float(np.mean([np.mean(sub[n] == sub[i]) for i, n in enumerate(nb)]))
    m6 = np.where(sub == M6)[0]
    row["M6_neighbour_purity"] = float(np.mean([np.mean(sub[nb[i]] == M6) for i in m6]))
    f1 = 0.0
    for c in clus.values():
        for k in np.unique(c[m6]):
            inc = c == k
            tp = np.sum(inc[m6]); prec = tp / inc.sum(); rec = tp / len(m6)
            f1 = max(f1, 2 * prec * rec / (prec + rec))
    row["M6_best_cluster_F1"] = f1
    # batch: neighbour sample entropy vs expectation under perfect mixing inside the nucleus's cell type
    dcodes = pd.Categorical(donor).codes
    nd = dcodes.max() + 1
    expected = {}
    for c in np.unique(ct):
        pool = dcodes[ct == c]
        expected[c] = np.mean([entropy(np.bincount(rng.choice(pool, K - 1), minlength=nd)) for _ in range(200)])
    ratio = [min(1.0, entropy(np.bincount(dcodes[n], minlength=nd)) / max(expected[ct[i]], 1e-9)) for i, n in enumerate(nb)]
    row["batch_neighbour_entropy"] = float(np.mean(ratio))
    gc = []
    for c in np.unique(ct):
        idx = np.where(ct == c)[0]
        _, comp = connected_components(C[idx][:, idx], directed=False)
        gc.append(np.bincount(comp).max() / len(idx))
    row["batch_graph_connectivity"] = float(np.mean(gc))
    row["jaccard_with_RNA_graph"] = float(np.mean([len(np.intersect1d(a, b)) / len(np.union1d(a, b)) for a, b in zip(nb, nb_rna)]))
    row["jaccard_with_ATAC_graph"] = float(np.mean([len(np.intersect1d(a, b)) / len(np.union1d(a, b)) for a, b in zip(nb, nb_atac)]))
    print(f"  {name}: " + ", ".join(f"{k} {v:.3f}" for k, v in row.items() if k != "method"), flush=True)
    return row


def main():
    rng = np.random.default_rng(SEED)
    rna = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    atac = sc.read_h5ad(INTERIM / "gse243917_atac.h5ad")
    common = rna.obs_names.intersection(atac.obs_names)
    rna, atac = rna[common].copy(), atac[common].copy()
    print(f"Nuclei: {len(common)}")

    Z, W, feats, r2 = train_mofa(rna, atac)
    Z = Z.loc[common]
    Zh = run_harmony(Z.values, rna.obs, "donor")
    pd.DataFrame(Zh, index=common, columns=Z.columns).round(5).to_csv(INTERIM / "gse243917_mofa_factors.tsv.gz", sep="\t")

    pca_h, lsi_h = rna.obsm["X_pca_harmony"], atac.obsm["X_lsi_harmony"]
    zs = lambda m: (m - m.mean(0)) / (m.std(0) + 1e-9)
    concat = np.hstack([zs(pca_h) / np.sqrt(pca_h.shape[1]), zs(lsi_h) / np.sqrt(lsi_h.shape[1])])
    with h5py.File(INTERIM / "gse243917_multiome.h5mu", "r") as f:
        idx = f["obs/_index"][:].astype(str)
        def csr(g):
            return sp.csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g.attrs["shape"]))
        Dw, Cw = csr(f["obsp/wnn_distances"]), csr(f["obsp/wnn_connectivities"])
    order = pd.Index(idx).get_indexer(common)
    Dw, Cw = Dw[order][:, order], Cw[order][:, order]

    graphs = {"RNA PCA (uncorrected)": graph_from(rna.obsm["X_pca"]),
              "RNA PCA + Harmony": graph_from(pca_h),
              "ATAC LSI + Harmony": graph_from(lsi_h, metric="cosine"),
              "Concatenation + Harmony": graph_from(concat),
              "WNN (step 07)": (Dw, Cw),
              "MOFA+ + Harmony": graph_from(Zh)}
    nb_rna, nb_atac = nbr_lists(graphs["RNA PCA + Harmony"][0]), nbr_lists(graphs["ATAC LSI + Harmony"][0])
    truths = {"author_ct": rna.obs["author_ct"].astype(str).values,
              "author_subtype": rna.obs["author_subtype"].astype(str).values,
              "ATAC_leiden": atac.obs["atac_leiden"].astype(str).values}
    ct, donor = truths["author_ct"], rna.obs["donor"].astype(str).values
    rows = [evaluate(n, D, C, truths, ct, donor, nb_rna, nb_atac, rng) for n, (D, C) in graphs.items()]
    res = pd.DataFrame(rows).set_index("method")

    bio = ["NMI_author_ct", "NMI_author_subtype", "NMI_ATAC_leiden", "ARI_author_ct", "ARI_author_subtype", "ARI_ATAC_leiden",
           "knn_agreement_subtype", "M6_neighbour_purity", "M6_best_cluster_F1"]
    batch = ["batch_neighbour_entropy", "batch_graph_connectivity"]
    mm = lambda s: (s - s.min()) / (s.max() - s.min() + 1e-12)
    res["bio_score"] = res[bio].apply(mm).mean(axis=1)
    res["batch_score"] = res[batch].apply(mm).mean(axis=1)
    res["overall_scIB_weights"] = 0.6 * res["bio_score"] + 0.4 * res["batch_score"]
    res = res.sort_values("overall_scIB_weights", ascending=False)
    res.round(4).to_csv(TAB / "integration_benchmark.tsv", sep="\t")
    pd.set_option("display.width", 250)
    print(res.round(3).T.to_string())

    # MOFA+ interpretation: variance explained per view, and the factor that tracks the IFN macrophage programme
    ve = r2.copy(); ve.index = Z.columns[:len(ve)]
    sub = truths["author_subtype"]
    gi = [rna.var_names.get_loc(g) for g in IFN if g in rna.var_names]
    X = rna.X[:, gi]
    X = X.toarray() if sp.issparse(X) else X
    ifn = ((X - X.mean(0)) / (X.std(0) + 1e-9)).mean(1)
    mac = (ct == "myeloid")
    cor = pd.Series({f: np.corrcoef(Z.loc[common, f].values[mac], ifn[mac])[0, 1] for f in Z.columns})
    fbest = cor.abs().idxmax(); sign = np.sign(cor[fbest]); k = int(fbest.replace("Factor", "")) - 1
    tss = tss_table()
    top_rows = []
    for view in W:
        w = W[view][:, k] * sign
        o = np.argsort(-w)[:25]
        names = feats[view][o]
        genes = names if view == "rna" else nearest_gene(names, tss)
        top_rows += [dict(factor=fbest, view=view, rank=i + 1, feature=n, nearest_gene=g, weight=float(w[j]))
                     for i, (n, g, j) in enumerate(zip(names, genes, o))]
    top = pd.DataFrame(top_rows)
    top.to_csv(TAB / "mofa_factors_top_features.tsv", sep="\t", index=False)
    in_m6 = Z.loc[common, fbest].values * sign
    print(f"\nMOFA+ factor tracking the IFN programme in myeloid nuclei: {fbest} (r = {cor[fbest]:.2f}); "
          f"median score M-6 {np.median(in_m6[sub == M6]):.2f} vs other myeloid {np.median(in_m6[mac & (sub != M6)]):.2f}")
    print("  top genes:", ", ".join(top[top.view == "rna"].feature.head(15)))
    print("  top peaks' nearest genes:", ", ".join(top[top.view == "atac"].nearest_gene.head(15)))
    ve.round(3).to_csv(TAB / "mofa_variance_explained.tsv", sep="\t")
    print(f"Variance explained (%) by {fbest}:", ve.loc[fbest].round(2).to_dict(), "; first 6 factors:\n", ve.head(6).round(2).to_string())

    fig, axes = plt.subplots(1, 2, figsize=(18, 5.5), gridspec_kw={"width_ratios": [3, 1]})
    cols = bio + batch + ["jaccard_with_RNA_graph", "jaccard_with_ATAC_graph", "overall_scIB_weights"]
    M = res[cols]
    ax = axes[0]
    S = M.apply(mm).values
    ax.imshow(S, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(cols)), [c.replace("_", " ") for c in cols], rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(M)), M.index, fontsize=9)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M.values[i, j]:.2f}", ha="center", va="center", fontsize=7, color="k" if S[i, j] > 0.6 else "w")
    ax.set_title("Integration benchmark (colour = scaled within column; text = raw value)", fontsize=10)
    ax = axes[1]
    order_s = pd.Series(in_m6).groupby(sub).median().loc[lambda s: s.index.str.startswith("M-")].sort_values()
    ax.barh(range(len(order_s)), order_s.values, color=["C3" if s == M6 else "C0" for s in order_s.index])
    ax.set_yticks(range(len(order_s)), order_s.index, fontsize=7)
    ax.set_xlabel(f"MOFA+ {fbest} (median per myeloid subtype)")
    ax.set_title("MOFA+ IFN factor", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "integration_benchmark.png", dpi=150)


if __name__ == "__main__":
    from pipeline_step_lock import step_lock
    with step_lock("28"):
        main()
