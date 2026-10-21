"""Spatial test of the IFNG -> IFN-activated macrophage hypothesis (Xenium, GSE347958; Sisó et al. 2026 Front Immunol).

Hypothesis fixed before analysis (H1, primary, 5K-panel RA sections): a macrophage's IFN programme score rises with
the number of IFNG+ T/NK cells within 30 um, beyond what the number of IFNG- T/NK cells within 30 um explains.
  Model per section: IFN score ~ n(IFNG+ lymphocytes within 30 um) + n(IFNG- lymphocytes within 30 um) (OLS).
  Null: IFNG labels permuted among the section's lymphocytes (positions fixed), 1,000 times -> empirical p for the
  IFNG+ coefficient. Section effects and empirical P values are reported separately; no normal-tail extrapolation of effect z scores.
Secondary (all 17 sections, 480-gene panel lacks IFNG/STAT1): share of CXCL9/CXCL10-high macrophages, RA vs NH/OA.
Cell types: per-section Leiden clusters named by marker scores (as in step 02).

Input : data/raw/GSE347958/* (per-sample matrix/features/barcodes/cells.parquet)
Output: results/tables/xenium_sections.tsv, xenium_ifng_neighbourhood.tsv, xenium_cxcl9_macrophages.tsv
        results/figures/xenium_spatial.png
"""
import glob
import gzip
import io
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import statsmodels.api as sm
from scipy.spatial import cKDTree
from scipy.stats import norm, mannwhitneyu

from config import RAW, PROCESSED, FIG, TAB, SEED

D = RAW / "GSE347958"
RADIUS = 30.0
N_PERM = 1000
MARKERS = {
    "macrophage": ["CD68", "CD163", "CSF1R", "CD14", "MERTK", "FOLR2", "LYVE1", "C1QA", "MS4A7"],
    "T": ["CD3E", "CD2", "CD4", "CD8A", "CD3D", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "PRF1", "KLRF1"],
    "B": ["MS4A1", "CD79A", "BANK1"],
    "plasma": ["MZB1", "JCHAIN", "TNFRSF17"],
    "fibroblast": ["PDGFRA", "PRG4", "COL1A1", "DCN", "FAP", "THY1"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "PLVAP"],
    "mural": ["ACTA2", "MYH11", "RGS5", "NOTCH3"],
}
IFN_GENES = ["CXCL9", "CXCL10", "CXCL11", "STAT1", "GBP1", "GBP5", "IDO1", "IRF1", "CD274", "APOL6"]


def sections():
    out = []
    for m in sorted(glob.glob(str(D / "*_matrix.mtx.gz"))):
        stem = m[:-len("_matrix.mtx.gz")]
        name = re.search(r"(CT_[\w-]+?)_(?:HU_)?(OA|RA|NORMAL)_SYNOVIUM", stem)
        out.append(dict(stem=stem, gsm=re.search(r"(GSM\d+)", stem).group(1), block=name.group(1),
                        condition={"NORMAL": "NH"}.get(name.group(2), name.group(2)),
                        panel="5K" if "_HU_" in stem else "480"))
    return out


def load(sec):
    X = scipy.io.mmread(sec["stem"] + "_matrix.mtx.gz").T.tocsr()
    feats = pd.read_csv(sec["stem"] + "_features.tsv.gz", sep="\t", header=None)
    bcs = pd.read_csv(sec["stem"] + "_barcodes.tsv.gz", header=None)[0].values
    keep = (feats[2] == "Gene Expression").values
    a = sc.AnnData(X=X[:, keep].astype(np.float32), obs=pd.DataFrame(index=bcs),
                   var=pd.DataFrame(index=feats.loc[keep, 1].values))
    a.var_names_make_unique()
    with gzip.open(sec["stem"] + "_cells.parquet.gz", "rb") as fh:
        cells = pd.read_parquet(io.BytesIO(fh.read())).set_index("cell_id")
    a = a[a.obs_names.isin(cells.index)].copy()
    a.obs[["x", "y"]] = cells.loc[a.obs_names, ["x_centroid", "y_centroid"]].values
    return a


def annotate(a, panel):
    a.layers["counts"] = a.X.copy()
    sc.pp.filter_cells(a, min_counts=20 if panel == "5K" else 10)
    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    sc.pp.pca(a, n_comps=30, random_state=SEED)
    sc.pp.neighbors(a, n_neighbors=15, random_state=SEED)
    sc.tl.leiden(a, resolution=1.0, flavor="igraph", n_iterations=2, random_state=SEED)
    scores = {}
    for ct, genes in MARKERS.items():
        g = [x for x in genes if x in a.var_names]
        if g:
            sc.tl.score_genes(a, g, score_name=f"s_{ct}", random_state=SEED)
            scores[ct] = f"s_{ct}"
    cl = a.obs.groupby("leiden", observed=True)[list(scores.values())].mean()
    z = (cl - cl.mean()) / cl.std()
    a.obs["cell_type"] = a.obs["leiden"].map(z.idxmax(axis=1).str.replace("s_", "")).astype(str)
    return a


def ifng_test(a, rng):
    mac = a.obs["cell_type"].eq("macrophage").values
    lym = a.obs["cell_type"].isin(["T", "NK"]).values
    ifng = np.asarray(a[:, "IFNG"].layers["counts"].todense()).ravel() > 0
    genes = [g for g in IFN_GENES if g in a.var_names]
    Xm = a[mac][:, genes].X
    Xm = np.asarray(Xm.todense()) if hasattr(Xm, "todense") else Xm
    score = ((Xm - Xm.mean(0)) / (Xm.std(0) + 1e-9)).mean(1)
    xy = a.obs[["x", "y"]].values
    tree = cKDTree(xy[lym])
    nb = tree.query_ball_point(xy[mac], r=RADIUS)
    lym_ifng = ifng[lym]

    def coef(labels):
        pos = np.array([labels[i].sum() for i in nb])
        neg = np.array([len(i) - labels[i].sum() for i in nb])
        if not pos.any():
            # No positive-exposure contrast is estimable. Keep the negative-neighbour
            # descriptive coefficient and consume the same RNG permutations below.
            fit = sm.OLS(score, sm.add_constant(neg, has_constant='add')).fit()
            return np.nan, fit.params[1], pos, neg
        fit = sm.OLS(score, sm.add_constant(np.column_stack([pos, neg]))).fit()
        return fit.params[1], fit.params[2], pos, neg

    b_pos, b_neg, pos, neg = coef(lym_ifng)
    null = np.array([coef(rng.permutation(lym_ifng))[0] for _ in range(N_PERM)])
    p = (1 + np.sum(null >= b_pos)) / (1 + N_PERM)
    zval = (b_pos - null.mean()) / (null.std() + 1e-12)
    estimable=bool(pos.any())
    if not estimable:p=np.nan;zval=np.nan
    near = pos > 0
    return dict(test_status='estimable' if estimable else 'not_estimable_no_IFNGpos_neighbours', n_macrophages=int(mac.sum()), n_lymphocytes=int(lym.sum()), pct_lymph_IFNG_pos=100 * lym_ifng.mean(),
                coef_IFNGpos_neighbour=b_pos, coef_IFNGneg_neighbour=b_neg, perm_null_mean=null.mean(), z=zval,
                p_perm_one_sided=p, mean_score_with_IFNGpos_neighbour=score[near].mean() if near.any() else np.nan,
                mean_score_without=score[~near].mean(), pct_mac_with_IFNGpos_neighbour=100 * near.mean(),
                ifn_genes_used=",".join(genes))


def main():
    rng = np.random.default_rng(SEED)
    secs = sections()
    rows, nb_rows, cx_rows = [], [], []
    example = None
    for sec in secs:
        a = annotate(load(sec), sec["panel"])
        if sec['block'] == 'CT_21045-A-15' and sec['panel'] == '5K':
            pos = np.asarray(a[:, 'CXCL13'].layers['counts'].todense()).ravel() > 0
            pd.DataFrame({'x': a.obs.x.to_numpy(), 'y': a.obs.y.to_numpy(),
                          'cell_type': a.obs.cell_type.to_numpy(), 'CXCL13': pos}).to_csv(
                PROCESSED / 'xenium_plot_CT_21045-A-15.tsv.gz', sep='\t', index=False,
                compression={'method': 'gzip', 'mtime': 0})
        comp = a.obs["cell_type"].value_counts(normalize=True)
        rows.append(dict(**{k: v for k, v in sec.items() if k != "stem"}, n_cells=a.n_obs,
                         **{f"frac_{k}": round(float(v), 4) for k, v in comp.items()}))
        mac = a[a.obs["cell_type"] == "macrophage"]
        cx = [g for g in ("CXCL9", "CXCL10") if g in mac.var_names]
        hi = (np.asarray(mac[:, cx].layers["counts"].sum(axis=1)).ravel() >= 2) if cx else np.array([])
        cx_rows.append(dict(gsm=sec["gsm"], block=sec["block"], condition=sec["condition"], panel=sec["panel"],
                            n_macrophages=mac.n_obs, pct_CXCL9_10_high_macrophages=100 * hi.mean() if len(hi) else np.nan))
        if sec["panel"] == "5K" and "IFNG" in a.var_names:
            r = ifng_test(a, rng)
            nb_rows.append(dict(gsm=sec["gsm"], block=sec["block"], condition=sec["condition"], **r))
            print(f"{sec['block']} {sec['condition']}: macrophages {r['n_macrophages']}, IFNG+ lymph {r['pct_lymph_IFNG_pos']:.1f}%, "
                  f"coef IFNG+ {r['coef_IFNGpos_neighbour']:.3f} vs IFNG- {r['coef_IFNGneg_neighbour']:.3f}, perm p {r['p_perm_one_sided']:.3f}")
            if sec["condition"] == "RA" and (example is None or r["n_macrophages"] > example[1]):
                example = (a, r["n_macrophages"], sec["block"])
        else:
            print(f"{sec['block']} {sec['condition']} ({sec['panel']}): {a.n_obs} cells, macrophages {mac.n_obs}")

    pd.DataFrame(rows).to_csv(TAB / "xenium_sections.tsv", sep="\t", index=False)
    nbh = pd.DataFrame(nb_rows)
    nbh.round(5).to_csv(TAB / "xenium_ifng_neighbourhood.tsv", sep="\t", index=False)
    cxt = pd.DataFrame(cx_rows)
    cxt.round(3).to_csv(TAB / "xenium_cxcl9_macrophages.tsv", sep="\t", index=False)

    ra = nbh[nbh["condition"] == "RA"]
    n_support=int((ra.p_perm_one_sided < 0.05).sum())
    print(f"\nH1 (RA, 5K panel, n = {len(ra)} sections): positive IFNG+ neighbour coefficient in "
          f"{int((ra.coef_IFNGpos_neighbour > 0).sum())}/{len(ra)}; p < 0.05 in {int((ra.p_perm_one_sided < 0.05).sum())}/{len(ra)}; "
          "no combined tail P (effect z scores are descriptive)")
    c480 = cxt[cxt.panel == "480"]
    u = mannwhitneyu(c480.loc[c480.condition == "RA", "pct_CXCL9_10_high_macrophages"],
                     c480.loc[c480.condition != "RA", "pct_CXCL9_10_high_macrophages"])
    print(f"Secondary (480 panel): CXCL9/10-high macrophages RA median "
          f"{c480.loc[c480.condition == 'RA', 'pct_CXCL9_10_high_macrophages'].median():.1f}% vs NH/OA "
          f"{c480.loc[c480.condition != 'RA', 'pct_CXCL9_10_high_macrophages'].median():.1f}% (Mann-Whitney p = {u.pvalue:.3f})")
    pd.DataFrame([dict(test="H1 IFNG+ neighbour effect on macrophage IFN score (RA 5K sections)", n_sections=len(ra),
                       n_positive=int((ra.coef_IFNGpos_neighbour > 0).sum()), n_perm_p_lt_0_05=int((ra.p_perm_one_sided < 0.05).sum()),
                       summary_method="section effects and empirical P; no combined tail P"),
                  dict(test="Secondary CXCL9/10-high macrophage share RA vs NH/OA (480 panel)", n_sections=len(c480),
                       summary_method="Mann-Whitney", mann_whitney_p=u.pvalue)]) \
        .to_csv(TAB / "xenium_tests.tsv", sep="\t", index=False, float_format="%.4g")

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    if example is not None:
        a, _, blk = example
        ax = axes[0]
        colors = {"macrophage": "C1", "T": "C0", "NK": "C2", "fibroblast": "#cccccc", "endothelial": "C4"}
        for ct, col in colors.items():
            s = a.obs[a.obs.cell_type == ct]
            ax.scatter(s.x, s.y, s=0.3, c=col, label=ct, rasterized=True)
        ifn = a.obs.cell_type.isin(["T", "NK"]).values & (np.asarray(a[:, "IFNG"].layers["counts"].todense()).ravel() > 0)
        ax.scatter(a.obs.x[ifn], a.obs.y[ifn], s=4, c="red", label="IFNG+ T/NK", rasterized=True)
        ax.set_aspect("equal"); ax.legend(markerscale=6, fontsize=7, loc="upper right"); ax.set_title(f"{blk} (RA, Xenium 5K)", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    ax = axes[1]
    x = np.arange(len(nbh))
    ax.bar(x - 0.2, nbh["coef_IFNGpos_neighbour"], 0.4, label="per IFNG+ T/NK within 30 um", color="C3")
    ax.bar(x + 0.2, nbh["coef_IFNGneg_neighbour"], 0.4, label="per IFNG- T/NK within 30 um", color="C7")
    for i, r in nbh.iterrows():
        ax.text(i - 0.2, r["coef_IFNGpos_neighbour"], f"p={r['p_perm_one_sided']:.3f}", ha="center", va="bottom", fontsize=7)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x, [f"{b}\n{c}" for b, c in zip(nbh.block, nbh.condition)], fontsize=7)
    ax.set_ylabel("Change in macrophage IFN score (OLS coefficient)")
    ax.set_title(f"H1: IFNG+ vs IFNG- lymphocyte neighbours\n{n_support}/{len(ra)} RA sections: empirical P < 0.05", fontsize=10)
    ax.legend(fontsize=7)
    ax = axes[2]
    for i, cond in enumerate(["NH", "OA", "RA"]):
        v = c480.loc[c480.condition == cond, "pct_CXCL9_10_high_macrophages"]
        ax.scatter(np.full(len(v), i) + rng.uniform(-.1, .1, len(v)), v, s=40)
    ax.set_xticks(range(3), ["Healthy", "OA", "RA"]); ax.set_ylabel("% macrophages with >= 2 CXCL9/CXCL10 counts")
    ax.set_title(f"Secondary (480-gene panel): RA vs NH/OA p = {u.pvalue:.3f}", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "xenium_spatial.png", dpi=150)


if __name__ == "__main__":
    main()
