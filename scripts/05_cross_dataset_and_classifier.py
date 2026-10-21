"""Cross-dataset check and a leakage-safe supervised baseline.

A. Signature transfer: genes DE between active RA and remission in E-MTAB-8322 macrophages
   (scRNA) are scored in GSE243917 macrophage pseudobulk (snRNA, independent donors) and
   correlated with synovial lymphocytic-inflammation grade. Nothing from GSE243917 is used
   to build the signature.
B. Classifier: active RA vs remission in E-MTAB-8322 at donor level (n = 19), leave-one-donor-out,
   using only a priori macrophage features (literature gene programs, label-free clustering),
   compared with a majority-class baseline and a single-feature model; label-permutation p-value.

Input : results/tables/emtab8322_DE_active_vs_remission.tsv, data/processed/emtab8322_donor_features.tsv,
        data/interim/gse243917_annotated.h5ad, data/processed/gse243917_donor_features.tsv
Output: results/tables/signature_genes.tsv, cross_dataset_signature_scores.tsv, classifier_loocv.tsv
        results/figures/cross_dataset_signature.png, classifier_loocv.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, balanced_accuracy_score, roc_curve
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from config import INTERIM, PROCESSED, FIG, TAB, SEED

N_SIG = 50
MIN_MAC = 50
DC_LIKE = ["M-9", "M-10", "M-11", "M-12", "M-13", "M-14"]


def build_signatures():
    """Primary: genes up in active RA at FDR<0.05 (used one-sided if the other direction has <10 FDR hits).
    Sensitivity: top-50 per direction by Wald statistic (nominal), as originally pre-specified."""
    de = pd.read_csv(TAB / "emtab8322_DE_active_vs_remission.tsv", sep="\t", index_col=0).dropna(subset=["pvalue"])
    de = de[~de.index.str.match(r"^(MT-|RP[SL]\d|ENSG)")]
    fdr_up = de[(de["padj"] < 0.05) & (de["log2FoldChange"] > 0.5)].sort_values("stat", ascending=False)
    fdr_dn = de[(de["padj"] < 0.05) & (de["log2FoldChange"] < -0.5)].sort_values("stat")
    sigs = {"primary_FDR": (fdr_up.index.tolist() if len(fdr_up) >= 10 else [],
                            fdr_dn.index.tolist() if len(fdr_dn) >= 10 else [])}
    nom_up = de[de["log2FoldChange"] > 0.5].sort_values("stat", ascending=False).head(N_SIG)
    nom_dn = de[de["log2FoldChange"] < -0.5].sort_values("stat").head(N_SIG)
    sigs["sensitivity_nominal50x50"] = (nom_up.index.tolist(), nom_dn.index.tolist())
    rows = [dict(signature=k, gene=g, direction=d) for k, (u, dn) in sigs.items()
            for d, gl in (("up_in_active", u), ("up_in_remission", dn)) for g in gl]
    pd.DataFrame(rows).merge(de[["log2FoldChange", "stat", "pvalue", "padj"]], left_on="gene", right_index=True) \
        .to_csv(TAB / "signature_genes.tsv", sep="\t", index=False)
    for k, (u, dn) in sigs.items():
        print(f"{k}: {len(u)} up in active, {len(dn)} up in remission")
        print("  up in active   :", ", ".join(u[:20]))
        if dn:
            print("  up in remission:", ", ".join(dn[:20]))
    return sigs


def score_gse243917(sigs):
    adata = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    sub = adata.obs["author_subtype"].astype(str).str.split(":").str[0]
    mac = adata[(adata.obs["author_ct"] == "myeloid") & ~sub.isin(DC_LIKE)]
    n_mac = mac.obs["donor"].value_counts()
    donors = sorted(n_mac[n_mac >= MIN_MAC].index)     # pseudobulk needs enough cells per donor
    C = mac.layers["counts"]
    pb = np.vstack([np.asarray(C[(mac.obs["donor"] == d).values].sum(axis=0)).ravel() for d in donors])
    pb = pd.DataFrame(pb, index=donors, columns=mac.var_names)
    cpm = np.log2(pb.div(pb.sum(axis=1), axis=0) * 1e6 + 1)
    z = (cpm - cpm.mean()) / cpm.std().replace(0, np.nan)
    feats = pd.read_csv(PROCESSED / "gse243917_donor_features.tsv", sep="\t", index_col=0)
    out = feats[["disease", "lymphocytic_inflammation_grade", "lining_hyperplasia_grade", "prep_batch", "n_macrophages"]].copy()
    stats = {}
    for k, (up, dn) in sigs.items():
        up_p, dn_p = [g for g in up if g in z.columns], [g for g in dn if g in z.columns]
        score = z[up_p].mean(axis=1) - (z[dn_p].mean(axis=1) if dn_p else 0)
        out[f"score_{k}"] = score.reindex(out.index)
        ra = out[(out["disease"] == "RA") & out[f"score_{k}"].notna()]
        rho, p = spearmanr(ra[f"score_{k}"], ra["lymphocytic_inflammation_grade"])
        rho_l, p_l = spearmanr(ra[f"score_{k}"], ra["lining_hyperplasia_grade"])
        p_b = mannwhitneyu(ra.loc[ra.prep_batch == "lowMT", f"score_{k}"],
                           ra.loc[ra.prep_batch == "highMT", f"score_{k}"]).pvalue
        stats[k] = dict(signature=k, n_RA=len(ra), genes_up=f"{len(up_p)}/{len(up)}", genes_down=f"{len(dn_p)}/{len(dn)}",
                        rho_inflammation=rho, p_inflammation=p, rho_lining=rho_l, p_lining=p_l, p_batch=p_b,
                        OA_score=out.loc[out["disease"] == "OA", f"score_{k}"].iloc[0], RA_median=ra[f"score_{k}"].median())
    out.to_csv(TAB / "cross_dataset_signature_scores.tsv", sep="\t")
    st = pd.DataFrame(stats.values())
    st.round(4).to_csv(TAB / "cross_dataset_signature_stats.tsv", sep="\t", index=False)
    print(f"GSE243917 transfer (donors with >= {MIN_MAC} macrophages):")
    print(st.round(3).to_string(index=False))
    return out, st.set_index("signature")


def loocv(X, y, model_fn):
    prob = np.zeros(len(y))
    for i in range(len(y)):
        tr = np.arange(len(y)) != i
        m = model_fn().fit(X[tr], y[tr])
        prob[i] = m.predict_proba(X[i:i + 1])[0, 1]
    return prob


def classifier():
    f = pd.read_csv(PROCESSED / "emtab8322_donor_features.tsv", sep="\t", index_col=0)
    f = f[f["group"].isin(["RA_naive_active", "RA_resistant_active", "RA_remission"])]
    y = (f["group"] != "RA_remission").astype(int).values   # 1 = active
    model = lambda: make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=1000))
    sets = {
        "single: mean_score_homeostatic": ["mean_score_homeostatic"],
        "a priori programs (3 features)": ["mean_score_homeostatic", "mean_score_inflammatory",
                                           "mac_frac_inflammatory"],
    }
    rng = np.random.default_rng(SEED)
    rows, curves = [], {}
    maj = np.full(len(y), y.mean())
    rows.append(dict(model="baseline: majority class", n=len(y), auc=0.5,
                     balanced_accuracy=0.5, accuracy=max(y.mean(), 1 - y.mean()), perm_p=np.nan))
    for name, cols in sets.items():
        X = f[cols].values
        prob = loocv(X, y, model)
        auc = roc_auc_score(y, prob)
        null = []
        for _ in range(500):
            yp = rng.permutation(y)
            null.append(roc_auc_score(yp, loocv(X, yp, model)))
        perm_p = (1 + np.sum(np.array(null) >= auc)) / (1 + len(null))
        pred = (prob >= 0.5).astype(int)
        rows.append(dict(model=name, n=len(y), auc=auc, balanced_accuracy=balanced_accuracy_score(y, pred),
                         accuracy=(pred == y).mean(), perm_p=perm_p))
        curves[name] = roc_curve(y, prob)
    res = pd.DataFrame(rows)
    res.round(3).to_csv(TAB / "classifier_loocv.tsv", sep="\t", index=False)
    print(res.round(3).to_string(index=False))
    return res, curves


def main():
    sigs = build_signatures()
    scores, stats = score_gse243917(sigs)
    res, curves = classifier()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    rng = np.random.default_rng(SEED)
    col = "score_primary_FDR"
    r = stats.loc["primary_FDR"]
    for b, mk in (("highMT", "o"), ("lowMT", "s")):
        s = scores[(scores["prep_batch"] == b) & scores[col].notna()]
        ax.scatter(s["lymphocytic_inflammation_grade"] + rng.uniform(-.08, .08, len(s)), s[col],
                   marker=mk, s=55, label=b, c=["grey" if d == "OA" else "C3" for d in s["disease"]])
    ax.set_xlabel("Lymphocytic inflammation grade (GSE243917; OA grey)")
    ax.set_ylabel("Active-RA macrophage score\n(FDR genes from E-MTAB-8322)")
    ax.set_title(f"Signature transfer, RA n={int(r.n_RA)}: rho={r.rho_inflammation:.2f}, p={r.p_inflammation:.3f}", fontsize=10)
    ax.legend(fontsize=8)
    ax = axes[1]
    for name, (fpr, tpr, _) in curves.items():
        a = res.set_index("model").loc[name, "auc"]
        ax.plot(fpr, tpr, label=f"{name} (AUC {a:.2f})")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8, label="majority baseline (AUC 0.50)")
    ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
    ax.set_title("Active vs remission, leave-one-donor-out (n=19)", fontsize=10)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "cross_dataset_and_classifier.png", dpi=150)


if __name__ == "__main__":
    main()
