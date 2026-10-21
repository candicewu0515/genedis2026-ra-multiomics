"""GSE243917: per-donor cell-composition and macrophage-state features vs synovial histology.

Statistics are donor-level (n = 11 RA; the OA sample is shown but excluded from correlations).
Composition is CLR-transformed. Every feature is also tested against the inferred nuclei-prep
batch, because the low-mito batch lost most fibroblasts.

Input : data/interim/gse243917_annotated.h5ad
Output: data/processed/gse243917_donor_features.tsv
        results/tables/gse243917_feature_vs_histology.tsv
        results/figures/gse243917_composition.png, gse243917_macrophage_index.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from config import INTERIM, PROCESSED, FIG, TAB, MAC_HOMEOSTATIC, MAC_INFLAMMATORY

IMMUNE = ["Tcell", "NK", "Bcell", "plasma", "myeloid"]
MERTK_POS = ["M-0", "M-1", "M-2", "M-3"]          # author subsets: MERTK+ SELENOP+/S100A8+/HBEGF+
MAC_INFL = ["M-4", "M-6", "M-7"]                   # SPP1+, STAT1+ CXCL10+, IL1B+ FCN1+
DC_LIKE = ["M-9", "M-10", "M-11", "M-12", "M-13", "M-14"]


def clr(frac, pseudo=0.5):
    x = np.log(frac + pseudo / frac.sum(axis=1).max() / 100)
    return x.sub(x.mean(axis=1), axis=0)


def pseudobulk_score(adata, genes, donors):
    genes = [g for g in genes if g in adata.raw.var_names]
    X = adata.raw[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else X
    df = pd.DataFrame(X, columns=genes, index=adata.obs_names)
    df["donor"] = adata.obs["donor"].values
    m = df.groupby("donor").mean().reindex(donors)
    return ((m - m.mean()) / m.std()).mean(axis=1)


def main():
    adata = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    obs = adata.obs
    labels = obs.groupby("donor", observed=True)[["disease", "lymphocytic_inflammation_grade",
                                                  "lining_hyperplasia_grade", "prep_batch"]].first()
    donors = labels.index.tolist()

    counts = pd.crosstab(obs["donor"], obs["cell_type"]).reindex(donors)
    frac_all = counts.div(counts.sum(axis=1), axis=0)
    imm = counts[IMMUNE]
    frac_imm = imm.div(imm.sum(axis=1), axis=0)

    feats = pd.DataFrame(index=donors)
    feats["n_nuclei"] = counts.sum(axis=1)
    for c in frac_all.columns:
        feats[f"frac_all_{c}"] = frac_all[c]
    for c, v in clr(frac_imm).items():
        feats[f"clr_immune_{c}"] = v

    # Macrophage state from author subtypes (macrophages only, DCs excluded)
    mye = obs[obs["author_ct"] == "myeloid"].copy()
    mye["subset"] = mye["author_subtype"].astype(str).str.split(":").str[0]
    mac = mye[~mye["subset"].isin(DC_LIKE)]
    tab = pd.crosstab(mac["donor"], mac["subset"]).reindex(donors).fillna(0)
    feats["n_macrophages"] = tab.sum(axis=1)
    feats["mac_frac_MERTKpos"] = tab[MERTK_POS].sum(axis=1) / tab.sum(axis=1)
    feats["mac_frac_inflammatory"] = tab[MAC_INFL].sum(axis=1) / tab.sum(axis=1)
    feats["mac_homeostasis_index"] = np.log2((tab[MERTK_POS].sum(axis=1) + 1) / (tab[MAC_INFL].sum(axis=1) + 1))

    mac_cells = adata[adata.obs_names.isin(mac.index)]
    feats["mac_score_homeostatic"] = pseudobulk_score(mac_cells, MAC_HOMEOSTATIC, donors)
    feats["mac_score_inflammatory"] = pseudobulk_score(mac_cells, MAC_INFLAMMATORY, donors)

    # Fibroblast lining vs sublining (author subtypes); only useful where fibroblasts were captured
    fib = obs[obs["author_ct"] == "fibroblast"].copy()
    fib["subset"] = fib["author_subtype"].astype(str).str.split(":").str[0]
    ft = pd.crosstab(fib["donor"], fib["subset"]).reindex(donors).fillna(0)
    feats["n_fibroblasts"] = ft.sum(axis=1)
    feats["fib_frac_lining"] = ft[["F-0", "F-1"]].sum(axis=1) / ft.sum(axis=1)
    feats["fib_frac_NOTCH3_sublining"] = ft.get("F-7", 0) / ft.sum(axis=1)

    out = labels.join(feats)
    out.round(4).to_csv(PROCESSED / "gse243917_donor_features.tsv", sep="\t")

    # Donor-level association tests, RA only
    ra = out[out["disease"] == "RA"]
    test_cols = [c for c in feats.columns if not c.startswith("n_")]
    rows = []
    for c in test_cols:
        sub = ra[ra[c].notna()]
        if c.startswith("fib_"):
            sub = sub[sub["n_fibroblasts"] >= 50]
        r_inf, p_inf = spearmanr(sub[c], sub["lymphocytic_inflammation_grade"])
        r_lin, p_lin = spearmanr(sub[c], sub["lining_hyperplasia_grade"])
        lo, hi = sub.loc[sub["prep_batch"] == "lowMT", c], sub.loc[sub["prep_batch"] == "highMT", c]
        p_b = mannwhitneyu(lo, hi).pvalue if len(lo) > 1 and len(hi) > 1 else np.nan
        rows.append(dict(feature=c, n=len(sub), rho_inflammation=r_inf, p_inflammation=p_inf,
                         rho_lining=r_lin, p_lining=p_lin, p_batch=p_b))
    res = pd.DataFrame(rows)
    for col in ("p_inflammation", "p_lining"):
        res[col.replace("p_", "fdr_")] = multipletests(res[col].fillna(1), method="fdr_bh")[1]
    res["batch_confounded"] = res["p_batch"] < 0.05
    res = res.sort_values("p_inflammation")
    res.round(4).to_csv(TAB / "gse243917_feature_vs_histology.tsv", sep="\t", index=False)
    print(res.round(3).to_string(index=False))

    # Figure 1: composition per donor
    order = out.sort_values(["disease", "lymphocytic_inflammation_grade", "prep_batch"]).index
    fig, axes = plt.subplots(1, 2, figsize=(16, 5))
    frac_all.loc[order].plot.bar(stacked=True, ax=axes[0], width=0.85, colormap="tab10")
    axes[0].set_title("All nuclei")
    frac_imm.loc[order].plot.bar(stacked=True, ax=axes[1], width=0.85, colormap="Set2")
    axes[1].set_title("Immune compartment only")
    for ax in axes:
        ax.set_xticklabels([f"{d.replace('BRI-', '')}\n{out.loc[d, 'disease']} g{out.loc[d, 'lymphocytic_inflammation_grade']}\n"
                            f"{out.loc[d, 'prep_batch']}" for d in order], rotation=0, fontsize=7)
        ax.set_ylabel("Fraction"); ax.legend(bbox_to_anchor=(1, 1), fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "gse243917_composition.png", dpi=150)

    # Figure 2: macrophage homeostasis index vs inflammation grade
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, c in zip(axes, ["mac_homeostasis_index", "mac_score_inflammatory"]):
        for b, mk in (("highMT", "o"), ("lowMT", "s")):
            s = out[out["prep_batch"] == b]
            ax.scatter(s["lymphocytic_inflammation_grade"] + np.random.default_rng(0).uniform(-.08, .08, len(s)),
                       s[c], marker=mk, label=b, s=50,
                       c=["grey" if d == "OA" else "C3" for d in s["disease"]])
        r = res.set_index("feature").loc[c]
        ax.set_title(f"{c}\nRA: rho={r.rho_inflammation:.2f}, p={r.p_inflammation:.3f}")
        ax.set_xlabel("Lymphocytic inflammation grade (OA in grey)")
        ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "gse243917_macrophage_index.png", dpi=150)


if __name__ == "__main__":
    main()
