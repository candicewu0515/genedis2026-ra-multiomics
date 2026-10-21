"""Marker programmes of two synovial states, defined once in GSE243917 and reused unchanged in new cohorts (steps 40, 41).

  M6_programme   M-6 STAT1+ CXCL10+ vs the other macrophage/monocyte subtypes (dendritic cells excluded)
  Tph_programme  T-7 CD4+ Tph vs the other CD4+ T-cell subtypes (T-3 Tfh/Tph excluded from the reference)
Rule (fixed before any new cohort was scored): Wilcoxon on log-normalised RNA; log2 fold change > 1, detected in >= 10 % of
the state's nuclei, no mitochondrial / ribosomal genes; top 30 by Wilcoxon score.
Output: results/tables/state_markers_gse243917.tsv
"""
import numpy as np
import pandas as pd
import scanpy as sc

from config import INTERIM, TAB

N_TOP = 30


def derive(adata, group, reference, name):
    a = adata[adata.obs["author_subtype"].isin([group] + reference)].copy()
    a.obs["grp"] = np.where(a.obs["author_subtype"] == group, "state", "ref")
    sc.tl.rank_genes_groups(a, "grp", groups=["state"], reference="ref", method="wilcoxon", use_raw=True, pts=True)
    df = sc.get.rank_genes_groups_df(a, "state")
    ref = a[a.obs.grp == "ref"].raw.to_adata()
    df["pct_nz_reference"] = np.asarray((ref[:, df["names"].values].X > 0).mean(axis=0)).ravel()
    df = df[~df["names"].str.match(r"^(MT-|RPL|RPS)")]
    df = df[(df["logfoldchanges"] > 1) & (df["pct_nz_group"] >= 0.10)].sort_values("scores", ascending=False).head(N_TOP)
    return df.assign(programme=name, state=group, n_state=int((a.obs.grp == "state").sum()),
                     n_reference=int((a.obs.grp == "ref").sum()))


def markers(force=False):
    out = TAB / "state_markers_gse243917.tsv"
    if out.exists() and not force:
        return pd.read_csv(out, sep="\t")
    adata = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    sub = adata.obs["author_subtype"].astype(str)
    mac = sorted(s for s in sub.unique() if s.startswith("M-") and "DC" not in s and "LAMP3" not in s)
    cd4 = sorted(s for s in sub.unique() if s.startswith("T-") and "CD4+" in s)
    m6 = "M-6: STAT1+ CXCL10+"
    tph = "T-7: CD4+ Tph"
    res = pd.concat([
        derive(adata, m6, [s for s in mac if s != m6], "M6_programme"),
        derive(adata, tph, [s for s in cd4 if s not in (tph, "T-3: CD4+ Tfh/Tph")], "Tph_programme"),
    ])
    res = res.rename(columns={"names": "gene", "scores": "wilcoxon_score", "logfoldchanges": "log2fc",
                              "pvals_adj": "padj", "pct_nz_group": "frac_state", "pct_nz_reference": "frac_reference"})
    res = res[["programme", "state", "gene", "wilcoxon_score", "log2fc", "padj", "frac_state", "frac_reference",
               "n_state", "n_reference"]]
    res.to_csv(out, sep="\t", index=False, float_format="%.4g")
    return res


# ---------- scoring helpers shared by steps 40 and 41 ----------
def matched_sets(genes, mean_expr, n, rng, bins=20):
    """Random gene sets matched to `genes` on mean expression (quantile bins of the universe)."""
    q = pd.qcut(mean_expr.rank(method="first"), bins, labels=False)
    pool = {b: q.index[q == b].to_numpy() for b in range(bins)}
    want = q.reindex(genes).dropna().astype(int)
    return [[rng.choice(pool[b]) for b in want] for _ in range(n)]


def unit_scores(logcpm, gene_sets):
    """logcpm: units x genes. Score = mean per-gene z across units."""
    z = (logcpm - logcpm.mean()) / logcpm.std().replace(0, np.nan)
    return pd.DataFrame({k: z[[g for g in v if g in z.columns]].mean(axis=1) for k, v in gene_sets.items()})


def residualise(y, groups):
    return y - y.groupby(groups).transform("mean")


def logcpm(pb):
    return np.log2(pb.div(pb.sum(axis=1), axis=0) * 1e6 + 1)


def cluster_best_state(X, obs, genes, gene_set, batch_key, donor_key, rng, n_random, seed=0):
    """Cluster cells (Harmony on batch_key, Leiden 1.0) and ask whether the cluster with the highest programme score
    stands out more than the best cluster of n_random expression-matched random gene sets.
    Score per cell = mean log-normalised expression of the set minus that of a matched control set (2x size)."""
    import harmonypy
    import scanpy as sc
    a = sc.AnnData(X.astype(np.float32), obs=obs.copy())
    a.var_names = pd.Index(np.asarray(genes, dtype=str), name=None); a.var_names_make_unique()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    hv = a.copy()
    sc.pp.highly_variable_genes(hv, n_top_genes=2000, batch_key=batch_key)
    hv = hv[:, hv.var.highly_variable].copy()
    sc.pp.scale(hv, max_value=10); sc.tl.pca(hv, n_comps=30, random_state=seed)
    ho = harmonypy.run_harmony(hv.obsm["X_pca"], hv.obs, [batch_key], max_iter_harmony=20, random_state=seed, verbose=False)
    zc = np.asarray(ho.Z_corr); zc = zc if zc.shape[0] == hv.n_obs else zc.T
    a.obsm["X_h"] = zc
    sc.pp.neighbors(a, use_rep="X_h", n_neighbors=15, random_state=seed)
    sc.tl.leiden(a, resolution=1.0, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    Xc = a.X.tocsc()
    gi = pd.Series(np.arange(a.n_vars), index=a.var_names)
    mean_expr = pd.Series(np.asarray(Xc.mean(axis=0)).ravel(), index=a.var_names)
    cl = a.obs["leiden"].to_numpy()

    def best(gs):
        ctrl = matched_sets(gs, mean_expr, 1, rng)[0] + matched_sets(gs, mean_expr, 1, rng)[0]
        s = np.asarray(Xc[:, gi[gs].values].mean(axis=1)).ravel() - np.asarray(Xc[:, gi[ctrl].values].mean(axis=1)).ravel()
        cm = pd.Series(s).groupby(cl).mean()
        return (cm.max() - s.mean()) / s.std(), cm.idxmax(), s

    gs = [g for g in gene_set if g in gi.index]
    stat, top, score = best(gs)
    null = np.array([best(r)[0] for r in matched_sets(gs, mean_expr, n_random, rng)])
    inb = cl == top
    sub = a.copy(); sub.obs["grp"] = np.where(inb, "best", "rest")
    sc.tl.rank_genes_groups(sub, "grp", groups=["best"], reference="rest", method="wilcoxon")
    mk = sc.get.rank_genes_groups_df(sub, "best").head(30).rename(columns={"names": "gene"})
    mk["in_programme"] = mk["gene"].isin(gs)
    per_donor = pd.Series(inb).groupby(a.obs[donor_key].to_numpy()).mean()
    row = dict(n=int(a.n_obs), n_clusters=int(len(set(cl))), statistic=float(stat), null_median=float(np.median(null)),
               p=(1 + (null >= stat).sum()) / (1 + n_random),
               detail=f"cluster {top}: {inb.mean():.1%} of cells; >=1% of cells in {(per_donor >= 0.01).sum()}/"
                      f"{per_donor.size} donors; {mk.in_programme.sum()}/30 top markers are programme genes")
    return row, mk, dict(score=score, best=inb, obs=a.obs)


def cell_scores(X, genes, gene_set, rng):
    """Per-cell score as in cluster_best_state: mean log-normalised expression of the set minus a matched control set."""
    import scipy.sparse as sp
    X = sp.csr_matrix(X, dtype=np.float32)
    lib = np.asarray(X.sum(axis=1)).ravel()
    Xn = sp.diags(1e4 / np.maximum(lib, 1)) @ X
    Xn.data = np.log1p(Xn.data)
    Xc = Xn.tocsc()
    names = pd.Index(genes).astype(str)
    gi = pd.Series(np.arange(len(names)), index=names)
    gi = gi[~gi.index.duplicated()]
    mean_expr = pd.Series(np.asarray(Xc.mean(axis=0)).ravel(), index=names)
    mean_expr = mean_expr[~mean_expr.index.duplicated()]
    gs = [g for g in gene_set if g in gi.index]
    ctrl = matched_sets(gs, mean_expr, 1, rng)[0] + matched_sets(gs, mean_expr, 1, rng)[0]
    return np.asarray(Xc[:, gi[gs].values].mean(axis=1)).ravel() - np.asarray(Xc[:, gi[ctrl].values].mean(axis=1)).ravel()


if __name__ == "__main__":
    r = markers(force=True)
    for p, g in r.groupby("programme"):
        print(p, f"({g.n_state.iloc[0]} vs {g.n_reference.iloc[0]} nuclei):", " ".join(g.gene))
