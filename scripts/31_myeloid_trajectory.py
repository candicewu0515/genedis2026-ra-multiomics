"""Where do STAT1+ CXCL10+ (M-6) macrophages sit among synovial myeloid states, and does chromatin open first?

Cells: GSE243917 myeloid nuclei without DC-like states (M-0..M-8; M-6 = 65 nuclei). RNA: 2,000 HVG, 30 PCs,
Harmony by sample, 30-NN graph; PAGA on author subtypes; diffusion pseudotime (DPT).
Hypotheses fixed before analysis:
  T1 M-6 is more strongly connected (PAGA) to the monocyte-derived IL1B+ FCN1+ state (M-7) than to any resident
     MERTK+ state (M-0..M-3). Robustness: leave-one-sample-out, and a PAGA graph built on chromatin (ATAC LSI).
  T2 (chromatin priming) non-M-6 nuclei carry more M-6 programme chromatin than RNA ("primed"), and more so in M-7
     than in resident MERTK+ nuclei. Programme = top 50 RNA markers of M-6 vs other macrophages (log2FC > 1,
     detected in >= 20 % of M-6). Per nucleus: RNA score = mean z of log expression; chromatin score = mean z of
     gene activity (ATAC counts in gene body + 2 kb upstream, as step 08, log depth-normalised); both smoothed over
     the RNA graph; priming = z(chromatin) - z(RNA). Test: per-sample mean priming, M-7 vs resident (paired across
     samples with >= 10 nuclei in both), against 500 random gene sets matched on RNA level and peak count.
  T3 (secondary) along DPT from the most monocyte-like M-7 nucleus, restricted to the M-7 -> M-6 path in the PAGA
     spanning tree, the chromatin score's final rise to its maximum crosses half its range earlier than the RNA
     score's (lag > 0), beyond random gene sets.
Roots are the nuclei with the highest neighbourhood-smoothed monocyte (M-7) or homeostatic (M-0) score; sanity check:
the root state should come first in pseudotime. Sensitivity: DPT from the M-0 root.

Output: results/tables/myeloid_paga.tsv, myeloid_trajectory_tests.tsv, myeloid_priming_by_state.tsv, myeloid_lag_sensitivity.tsv,
        myeloid_lag_controls.tsv
        results/figures/myeloid_trajectory.png
"""
import importlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.sparse.csgraph import shortest_path
from scipy.stats import wilcoxon

from config import INTERIM, FIG, TAB, SEED, MAC_HOMEOSTATIC

step08 = importlib.import_module("08_gse243917_gene_activity_links")
DC_LIKE = ["M-9", "M-10", "M-11", "M-12", "M-13", "M-14"]
RESIDENT = ["M-0", "M-1", "M-2", "M-3"]
MONO = ["FCN1", "S100A8", "S100A9", "VCAN", "LYZ", "CD14", "IL1B"]
N_RAND = 500


def run_harmony(x, obs, key):
    import harmonypy
    ho = harmonypy.run_harmony(x, obs, [key], max_iter_harmony=20, random_state=SEED, verbose=False)
    z = np.asarray(ho.Z_corr)
    return z if z.shape[0] == x.shape[0] else z.T


def paga_conn(ad, rep, metric="euclidean"):
    sc.pp.neighbors(ad, use_rep=rep, n_neighbors=30, metric=metric, random_state=SEED)
    sc.tl.paga(ad, groups="code")
    cats = ad.obs["code"].cat.categories
    return pd.DataFrame(ad.uns["paga"]["connectivities"].toarray(), index=cats, columns=cats)


def t1_stat(conn):
    m6 = conn.loc["M-6"].drop("M-6")
    res = m6[[c for c in RESIDENT if c in m6.index]]
    return m6.get("M-7", np.nan), res.max(), res.idxmax(), m6.idxmax()


def half(s, thr=0.5):
    """Pseudotime quantile where the final rise to the curve's maximum crosses `thr` of its range."""
    v = ((s - s.min()) / (s.max() - s.min() + 1e-12)).values
    peak = int(np.argmax(v))
    below = np.where(v[:peak + 1] < thr)[0]
    return (below[-1] + 1 if len(below) else 0) / (len(v) - 1)


def zc(x):
    return (x - x.mean(0)) / (x.std(0) + 1e-9)


def main():
    rng = np.random.default_rng(SEED)
    rna = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    code = rna.obs["author_subtype"].astype(str).str.split(":").str[0]
    keep = (rna.obs["author_ct"] == "myeloid").values & ~code.isin(DC_LIKE).values
    m = rna[keep].copy()
    del rna
    m.obs["code"] = pd.Categorical(code[keep].values)
    atac = sc.read_h5ad(INTERIM / "gse243917_atac.h5ad")
    a = atac[m.obs_names].copy()
    peaks = atac.var[["chrom", "start", "end"]].copy()
    del atac
    print(f"Myeloid nuclei (no DCs): {m.n_obs}; per state: {m.obs['code'].value_counts().to_dict()}")

    hv = sc.AnnData(m.layers["counts"].copy(), obs=m.obs[["donor"]].copy(), var=pd.DataFrame(index=m.var_names))
    sc.pp.normalize_total(hv, target_sum=1e4); sc.pp.log1p(hv)
    sc.pp.highly_variable_genes(hv, n_top_genes=2000, flavor="seurat", batch_key="donor")
    hv = hv[:, hv.var["highly_variable"]].copy()
    sc.pp.scale(hv, max_value=10)
    sc.tl.pca(hv, n_comps=30, random_state=SEED)
    m.obsm["X_h"] = run_harmony(hv.obsm["X_pca"], m.obs, "donor")
    a.obsm["X_lsi_sub"] = a.obsm["X_lsi_harmony"]
    a.obs["code"] = m.obs["code"].values

    # T1: PAGA on RNA, leave-one-sample-out, and on ATAC
    conn = paga_conn(m, "X_h")
    sc.tl.umap(m, random_state=SEED)
    m7, res_max, res_arg, top = t1_stat(conn)
    print(f"T1 RNA PAGA: M-6 to M-7 {m7:.3f}; to best resident ({res_arg}) {res_max:.3f}; strongest partner {top}")
    lodo = []
    for d in sorted(m.obs["donor"].unique()):
        sub = m[m.obs["donor"] != d].copy()
        sub.obs["code"] = sub.obs["code"].cat.remove_unused_categories()
        if (sub.obs["code"] == "M-6").sum() < 20:
            continue
        c = paga_conn(sub, "X_h")
        x7, xr, xa, xt = t1_stat(c)
        lodo.append(dict(left_out=d, M6_to_M7=x7, M6_to_best_resident=xr, best_resident=xa, strongest_partner=xt))
    lodo = pd.DataFrame(lodo)
    ca = paga_conn(a, "X_lsi_sub", metric="cosine")
    a7, ar, aa, at = t1_stat(ca)
    print(f"T1 ATAC PAGA: M-6 to M-7 {a7:.3f}; to best resident ({aa}) {ar:.3f}; strongest partner {at}")
    out = conn.loc[["M-6"]].T.rename(columns={"M-6": "rna_paga_to_M6"}).join(ca.loc[["M-6"]].T.rename(columns={"M-6": "atac_paga_to_M6"}))
    out.round(4).to_csv(TAB / "myeloid_paga.tsv", sep="\t")

    # programme genes and chromatin (gene activity)
    genes = step08.read_genes()
    genes = genes[genes.index.isin(m.var_names)]
    P2G = step08.peak_gene_matrix(peaks, genes)
    n_peaks = pd.Series(np.asarray(P2G.sum(axis=0)).ravel(), index=genes.index)
    C = a.layers["counts"].tocsr()
    depth = np.asarray(C.sum(axis=1)).ravel()
    GA = sp.csr_matrix(C @ P2G)
    GA = sp.diags(1e4 / depth) @ GA
    GA.data = np.log1p(GA.data)
    X = m[:, genes.index].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    is6 = (m.obs["code"] == "M-6").values
    lfc = np.log2((np.expm1(X[is6]).mean(0) + 1e-3) / (np.expm1(X[~is6]).mean(0) + 1e-3))
    det = (X[is6] > 0).mean(0)
    cand = pd.DataFrame({"lfc": lfc, "det": det, "n_peaks": n_peaks.values, "mean": X.mean(0)}, index=genes.index)
    prog = cand[(cand.lfc > 1) & (cand.det >= 0.2) & (cand.n_peaks > 0)].sort_values("lfc", ascending=False).index[:50].tolist()
    print(f"M-6 programme: {len(prog)} genes, e.g. {', '.join(prog[:15])}")

    W = m.obsp["connectivities"] + sp.eye(m.n_obs)
    W = sp.diags(1 / np.asarray(W.sum(axis=1)).ravel()) @ W
    gidx = {g: i for i, g in enumerate(genes.index)}
    GAd = GA.tocsc()
    Zr = zc(X)

    def scores(gl):
        ii = [gidx[g] for g in gl]
        r = Zr[:, ii].mean(1)
        g = GAd[:, ii].toarray()
        g = zc(g).mean(1)
        rs, gs = W @ r, W @ g
        return rs, gs, zc(gs) - zc(rs)

    def t2_stat(prim):
        df = pd.DataFrame({"p": prim, "donor": m.obs["donor"].values, "code": m.obs["code"].astype(str).values})
        df["grp"] = np.where(df.code == "M-7", "M-7", np.where(df.code.isin(RESIDENT), "resident", "other"))
        g = df[df.grp != "other"].groupby(["donor", "grp"], observed=True).p.agg(["mean", "size"]).unstack()
        g = g[(g["size"] >= 10).all(axis=1)]
        return g["mean"]["M-7"] - g["mean"]["resident"]

    rs, gs, prim = scores(prog)
    diff = t2_stat(prim)
    w = wilcoxon(diff)
    # matched random sets: RNA mean decile x peak-count tercile among genes with >= 1 peak, programme genes excluded
    allp = cand[cand.n_peaks > 0].copy()
    allp["bin"] = pd.qcut(allp["mean"].rank(method="first"), 10, labels=False) * 3 + \
        pd.qcut(allp["n_peaks"].rank(method="first"), 3, labels=False)
    want = allp.loc[prog, "bin"].values
    pool = allp.drop(index=prog)
    members = {b: pool.index[pool.bin == b].values for b in np.unique(want)}

    # T3: path M-7 -> M-6 in the PAGA spanning tree, DPT from the most monocyte-like M-7 nucleus
    cats = list(conn.index)
    tree = m.uns["paga"]["connectivities_tree"]
    tree = ((tree + tree.T) > 0).astype(float)
    _, pred = shortest_path(tree, method="D", unweighted=True, return_predecessors=True, indices=cats.index("M-7"))
    path, j = [], cats.index("M-6")
    while j != -9999 and j != cats.index("M-7"):
        path.append(cats[j]); j = pred[j]
    path = ["M-7"] + path[::-1]
    mono = zc(X[:, [gidx[g] for g in MONO if g in gidx]]).mean(1)
    homeo = zc(X[:, [gidx[g] for g in MAC_HOMEOSTATIC if g in gidx]]).mean(1)
    sc.tl.diffmap(m, random_state=SEED)
    # roots from neighbourhood-smoothed scores (a raw maximum picked an outlier nucleus far from its own state)
    mono_s, homeo_s = W @ mono, W @ homeo
    roots = {"M-7 monocyte-like": np.where(m.obs["code"] == "M-7")[0][np.argmax(mono_s[m.obs["code"] == "M-7"])],
             "M-0 homeostatic": np.where(m.obs["code"] == "M-0")[0][np.argmax(homeo_s[m.obs["code"] == "M-0"])]}
    dpt, sanity = {}, {}
    for name, r in roots.items():
        m.uns["iroot"] = int(r)
        sc.tl.dpt(m)
        dpt[name] = m.obs["dpt_pseudotime"].values.copy()
        med = pd.Series(dpt[name]).groupby(m.obs["code"].astype(str).values).median().sort_values()
        print(f"DPT from {name}: state order {', '.join(f'{k} {v:.2f}' for k, v in med.items())}")
        sanity[name] = med.index[0] == name.split()[0]
    on_path = m.obs["code"].astype(str).isin(path).values

    def lag(rs_, gs_, t, nb=20, thr=0.5):
        q = pd.qcut(pd.Series(t[on_path]).rank(method="first"), nb, labels=False).values
        cr = pd.Series(rs_[on_path]).groupby(q).mean(); cg = pd.Series(gs_[on_path]).groupby(q).mean()
        return half(cr, thr) - half(cg, thr), cr, cg
    t_main = dpt["M-7 monocyte-like"]
    lag_obs, cr, cg = lag(rs, gs, t_main)

    rand_t2, rand_lag, rand_sc = [], [], []
    for _ in range(N_RAND):
        gl = [rng.choice(members[b]) for b in want]
        r_, g_, p_ = scores(gl)
        rand_t2.append(t2_stat(p_).mean())
        rand_lag.append(lag(r_, g_, t_main)[0])
        rand_sc.append((r_, g_))
    rand_t2, rand_lag = np.array(rand_t2), np.array(rand_lag)
    # T3 sensitivity (added after the half-max definition was corrected): bins x crossing threshold
    grid = []
    for nb in (10, 20, 30):
        for thr in (0.4, 0.5, 0.6):
            o = lag(rs, gs, t_main, nb, thr)[0]
            rl = np.array([lag(r_, g_, t_main, nb, thr)[0] for r_, g_ in rand_sc])
            grid.append(dict(n_bins=nb, threshold=thr, lag=o, random_median=float(np.median(rl)),
                             p_vs_random=(1 + (rl >= o).sum()) / (1 + len(rl))))
    grid = pd.DataFrame(grid)
    grid.to_csv(TAB / "myeloid_lag_sensitivity.tsv", sep="\t", index=False, float_format="%.4g")
    print("T3 sensitivity grid:\n", grid.round(3).to_string(index=False))

    # T3 controls (added after seeing T3): RNA-selected markers rise late by construction, so (a) an a priori IFN list
    # not selected on these data, and (b) the same marker selection for other end states of the tree
    def tree_path(target):
        _, pr = shortest_path(tree, method="D", unweighted=True, return_predecessors=True, indices=cats.index("M-7"))
        pth, k = [], cats.index(target)
        while k != -9999 and k != cats.index("M-7"):
            pth.append(cats[k]); k = pr[k]
        return ["M-7"] + pth[::-1]

    def lag_on(rs_, gs_, mask, nb=20):
        q = pd.qcut(pd.Series(t_main[mask]).rank(method="first"), nb, labels=False).values
        cr = pd.Series(rs_[mask]).groupby(q).mean(); cg = pd.Series(gs_[mask]).groupby(q).mean()
        return half(cr) - half(cg)

    ctrl = []
    sets_ctrl = {"M-6 a priori IFN list": ("M-6", [g for g in ["CXCL9", "CXCL10", "CXCL11", "STAT1", "GBP1", "GBP5", "IDO1",
                                                                 "IRF1", "CD274", "APOL6"] if g in gidx and cand.loc[g, "n_peaks"] > 0])}
    for tgt in ("M-4", "M-0", "M-5"):
        it = (m.obs["code"] == tgt).values
        lf = np.log2((np.expm1(X[it]).mean(0) + 1e-3) / (np.expm1(X[~it & (m.obs["code"] != "M-6").values]).mean(0) + 1e-3))
        dt = (X[it] > 0).mean(0)
        c2 = pd.DataFrame({"lfc": lf, "det": dt, "n_peaks": n_peaks.values}, index=genes.index)
        sets_ctrl[f"{tgt} RNA markers (same selection)"] = (tgt, c2[(c2.lfc > 1) & (c2.det >= 0.2) & (c2.n_peaks > 0)]
                                                           .sort_values("lfc", ascending=False).index[:50].tolist())
    for name, (tgt, gl) in sets_ctrl.items():
        pth = tree_path(tgt)
        mask = m.obs["code"].astype(str).isin(pth).values
        r_, g_, _ = scores(gl)
        o = lag_on(r_, g_, mask)
        w_ = allp.loc[gl, "bin"].values
        mem = {b: allp.drop(index=gl).index[allp.drop(index=gl).bin == b].values for b in np.unique(w_)}
        rl = []
        for _ in range(200):
            rr, gg, _ = scores([rng.choice(mem[b]) for b in w_])
            rl.append(lag_on(rr, gg, mask))
        rl = np.array(rl)
        ctrl.append(dict(gene_set=name, n_genes=len(gl), path=" > ".join(pth), lag=o, random_median=float(np.median(rl)),
                         p_vs_random=(1 + (rl >= o).sum()) / (1 + len(rl))))
    ctrl = pd.DataFrame(ctrl)
    ctrl.to_csv(TAB / "myeloid_lag_controls.tsv", sep="\t", index=False, float_format="%.4g")
    print("T3 controls:\n", ctrl.round(3).to_string(index=False))

    by_state = pd.DataFrame({"code": m.obs["code"].astype(str).values, "donor": m.obs["donor"].values, "priming": prim,
                             "rna_score": rs, "chromatin_score": gs, "dpt_from_M7": t_main, "dpt_from_M0": dpt["M-0 homeostatic"]})
    st = by_state.groupby("code").agg(n=("priming", "size"), median_priming=("priming", "median"),
                                      median_rna=("rna_score", "median"), median_chromatin=("chromatin_score", "median"),
                                      median_dpt_from_M7=("dpt_from_M7", "median"), median_dpt_from_M0=("dpt_from_M0", "median"))
    st.round(4).to_csv(TAB / "myeloid_priming_by_state.tsv", sep="\t")
    print(st.round(3).to_string())

    rows = [dict(test="T1 PAGA (RNA): M-6 to M-7 minus M-6 to best resident", value=m7 - res_max,
                 note=f"M-7 {m7:.3f}; {res_arg} {res_max:.3f}; strongest partner {top}"),
            dict(test="T1 leave-one-sample-out: share of runs with M-7 > best resident", value=float((lodo.M6_to_M7 > lodo.M6_to_best_resident).mean()),
                 note=f"{len(lodo)} runs; strongest partners {lodo.strongest_partner.value_counts().to_dict()}"),
            dict(test="T1 PAGA (ATAC): M-6 to M-7 minus M-6 to best resident", value=a7 - ar,
                 note=f"M-7 {a7:.3f}; {aa} {ar:.3f}; strongest partner {at}"),
            dict(test="T2 priming M-7 minus resident, per-sample mean (Wilcoxon over samples)", value=float(diff.mean()),
                 p=w.pvalue, note=f"{len(diff)} samples; {int((diff > 0).sum())} positive"),
            dict(test="T2 vs random gene sets (share with mean difference >= observed)", value=float(np.median(rand_t2)),
                 p=(1 + (rand_t2 >= diff.mean()).sum()) / (1 + N_RAND), note="value = random median"),
            dict(test="T3 lag (RNA half-max minus chromatin half-max, pseudotime quantile) on path " + " > ".join(path),
                 value=lag_obs, p=(1 + (rand_lag >= lag_obs).sum()) / (1 + N_RAND), note=f"random median {np.median(rand_lag):.3f}"),
            dict(test="T3 sensitivity: grid cells (bins x threshold) with lag > 0 and p vs random < 0.05",
                 value=int(((grid.lag > 0) & (grid.p_vs_random < 0.05)).sum()), note=f"of {len(grid)}; lags {grid.lag.min():.2f} to {grid.lag.max():.2f}"),
            dict(test="T3 control: a priori IFN list (not selected on these data), lag on the M-6 path", value=ctrl.lag.iloc[0],
                 p=ctrl.p_vs_random.iloc[0], note=f"{ctrl.n_genes.iloc[0]} genes"),
            dict(test="T3 control: other end states with the same marker selection, lag > 0 and p < 0.05",
                 value=int(((ctrl.lag.iloc[1:] > 0) & (ctrl.p_vs_random.iloc[1:] < 0.05)).sum()),
                 note="; ".join(f"{r.gene_set.split()[0]} lag {r.lag:.2f} p {r.p_vs_random:.2g}" for r in ctrl.iloc[1:].itertuples())),
            dict(test="Sanity: root state comes first in pseudotime (M-7 root; M-0 root)", value=int(sanity["M-7 monocyte-like"]),
                 note=f"M-7 root {'passes' if sanity['M-7 monocyte-like'] else 'FAILS (T3 not interpretable)'}; "
                      f"M-0 root {'passes' if sanity['M-0 homeostatic'] else 'fails'}"),
            dict(test="Sensitivity: median DPT rank of M-6 (root M-7 vs root M-0)",
                 value=float(st.median_dpt_from_M7.rank().loc["M-6"]), note=f"from M-0: rank {st.median_dpt_from_M0.rank().loc['M-6']:.0f} of {len(st)}")]
    res = pd.DataFrame(rows)
    res.to_csv(TAB / "myeloid_trajectory_tests.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
    print(res.round(4).to_string(index=False))
    print("LOSO:\n", lodo.round(3).to_string(index=False))

    fig, axes = plt.subplots(1, 4, figsize=(24, 5.5))
    ax = axes[0]
    um = m.obsm["X_umap"]
    for i, c in enumerate(m.obs["code"].cat.categories):
        s = (m.obs["code"] == c).values
        ax.scatter(um[s, 0], um[s, 1], s=3 if c != "M-6" else 14, label=c, color="k" if c == "M-6" else f"C{i % 10}", alpha=0.7)
    ax.legend(markerscale=3, fontsize=7); ax.set_xticks([]); ax.set_yticks([]); ax.set_title("Myeloid nuclei (M-6 in black)", fontsize=10)
    ax = axes[1]
    v = conn.loc["M-6"].drop("M-6").sort_values()
    va = ca.loc["M-6"].reindex(v.index)
    y = np.arange(len(v))
    ax.barh(y + 0.2, v.values, height=0.4, label="RNA graph"); ax.barh(y - 0.2, va.values, height=0.4, label="ATAC graph")
    ax.set_yticks(y, v.index); ax.legend(fontsize=8); ax.set_xlabel("PAGA connectivity to M-6"); ax.set_title("T1: neighbours of M-6", fontsize=10)
    ax = axes[2]
    order = st.sort_values("median_priming").index
    for i, c in enumerate(order):
        d = by_state[by_state.code == c].groupby("donor").priming.mean()
        ax.scatter(np.full(len(d), i) + rng.uniform(-.1, .1, len(d)), d, s=15, color="C3" if c in ("M-7",) else "C0" if c in RESIDENT else "C7")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(order)), order); ax.set_ylabel("M-6 programme priming (chromatin z - RNA z), sample means")
    ax.set_title(f"T2: M-7 - resident {diff.mean():.2f} (p {w.pvalue:.2g}; vs random p {res.p.iloc[4]:.2g})", fontsize=10)
    ax = axes[3]
    xq = np.arange(len(cr)) / 19
    ax.plot(xq, (cr - cr.min()) / (cr.max() - cr.min()), label="RNA", color="C0")
    ax.plot(xq, (cg - cg.min()) / (cg.max() - cg.min()), label="chromatin (gene activity)", color="C3")
    ax.axhline(0.5, color="k", lw=0.6, ls="--"); ax.legend(fontsize=8)
    ax.set_xlabel("pseudotime quantile on path " + " > ".join(path)); ax.set_ylabel("scaled M-6 programme score")
    ax.set_title(f"T3: lag {lag_obs:.2f} (vs random p {res.p.iloc[5]:.2g})" +
                 ("" if sanity["M-7 monocyte-like"] else "\nM-7 root fails sanity check: not interpreted"), fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "myeloid_trajectory.png", dpi=150)


if __name__ == "__main__":
    main()
