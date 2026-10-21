"""Clinical-response test: does the active-RA macrophage signature relate to anti-TNF response?

Cohort: GSE198520, synovial RNA-seq from 46 RA patients before and 12 weeks after anti-TNF
(certolizumab n = 27, etanercept n = 19); EULAR response: good (n = 19), moderate (n = 13), non (n = 14).
(R4RA, the first choice, is not usable: processed expression and treatment-arm labels are not public.)

Signature (fixed before this cohort was opened): FDR genes up in active-RA macrophages, E-MTAB-8322 (step 05);
13 of 15 genes present. Score = mean per-gene z of log2(CPM + 1); z parameters from baseline samples only.

Hypotheses, written before outcome data were analysed (all two-sided):
  H1  baseline score vs disease activity: Spearman with DAS28 and histological inflammatory score;
      Kruskal-Wallis across pathotypes.
  H2  baseline score vs response: EULAR good vs non-responders (moderate excluded), Mann-Whitney + AUC
      (bootstrap 95% CI, label-permutation p); ordinal Spearman non < moderate < good (all 46).
      Comparators: majority class, baseline DAS28 alone, 1,000 random gene sets matched on expression decile.
  H3  change with treatment: paired pre->post change (Wilcoxon) overall and per response group; change in good vs
      non-responders (Mann-Whitney); Spearman of score change vs DAS28 change.
  Confounders: RIN (RNA integrity, 11 missing), biologic; logistic models response ~ score + biologic (+ RIN).

Added after the first run (disclosed in the report): the first principal component (~57 % of variance) tracked library
size, mitochondrial-read fraction and RIN, and also response; random gene sets gave AUC ~0.15. So the analysis also
reports (a) log library size and mito fraction as stand-alone "technical" predictors, (b) the signature score
residualised on log library size + mito fraction (fit label-free on all samples), (c) logistic models with both
covariates, and (d) whether the signature genes are more library-size-dependent than expression-matched random genes.
CPM uses nuclear genes only (mitochondrial reads excluded from the library size).

Input : data/raw/GSE198520/*, results/tables/signature_genes.tsv
Output: data/processed/gse198520_sample_scores.tsv
        results/tables/gse198520_tests.tsv, gse198520_auc.tsv
        results/figures/gse198520_response.png
"""
import gzip

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr, mannwhitneyu, kruskal, wilcoxon
from sklearn.metrics import roc_auc_score, roc_curve

from config import RAW, PROCESSED, FIG, TAB, SEED, MAC_HOMEOSTATIC, MAC_INFLAMMATORY

D = RAW / "GSE198520"
N_RAND = 1000
N_BOOT = 2000


def metadata():
    rows = {}
    with gzip.open(D / "GSE198520_series_matrix.txt.gz", "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                rows["title"] = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
                rows[vals[0].split(":")[0]] = [v.split(": ", 1)[1] if ": " in v else np.nan for v in vals]
    m = pd.DataFrame(rows)
    m["sample"] = m["title"] + "_bx"
    m["patient"] = m["title"].str.rsplit("_", n=1).str[0]
    m["response"] = m["status"].map({"Non-Responder": "non", "Mod-Responder": "moderate", "Responder": "good"})
    m["response_ord"] = m["response"].map({"non": 0, "moderate": 1, "good": 2})
    for c in ("das28 score", "delta das28", "inflammatory score", "rin", "crp", "age"):
        m[c] = pd.to_numeric(m[c].replace({"N/A": np.nan, "NA": np.nan}), errors="coerce")
    m["histology"] = m["aggregational histological score"].str.replace("Pauci immune", "Pauci-immune")
    return m.set_index("sample")


def auc_ci(y, s, rng):
    a = roc_auc_score(y, s)
    boots = []
    idx = np.arange(len(y))
    while len(boots) < N_BOOT:
        b = rng.choice(idx, len(idx), replace=True)
        if len(np.unique(y[b])) == 2:
            boots.append(roc_auc_score(y[b], s[b]))
    perm = np.array([roc_auc_score(rng.permutation(y), s) for _ in range(N_BOOT)])
    p_perm = (1 + np.sum(np.abs(perm - 0.5) >= abs(a - 0.5))) / (1 + len(perm))
    return a, np.quantile(boots, 0.025), np.quantile(boots, 0.975), p_perm


def main():
    rng = np.random.default_rng(SEED)
    counts = pd.read_csv(D / "GSE198520_Raw_gene_count_matrix.txt.gz", sep="\t", index_col=0)
    meta = metadata()
    counts = counts[meta.index]
    mito = counts.index.str.startswith("MT-")
    meta["log_libsize"] = np.log(counts.sum(axis=0))
    meta["pct_mito"] = counts[mito].sum(axis=0) / counts.sum(axis=0) * 100
    counts = counts[~mito]
    logcpm = np.log2(counts / counts.sum(axis=0) * 1e6 + 1)
    pre = meta.index[meta["timepoint"] == "Pre"]
    mu, sd = logcpm[pre].mean(axis=1), logcpm[pre].std(axis=1).replace(0, np.nan)
    z = logcpm.sub(mu, axis=0).div(sd, axis=0)

    sg = pd.read_csv(TAB / "signature_genes.tsv", sep="\t")
    sig = sg.loc[(sg["signature"] == "primary_FDR") & (sg["direction"] == "up_in_active"), "gene"].tolist()
    sig = [g for g in sig if g in z.index and np.isfinite(z.loc[g, pre]).all()]
    print(f"Signature genes available: {len(sig)} -> {', '.join(sig)}")
    meta["signature_score"] = z.loc[sig].mean(axis=0)
    for name, genes in (("homeostatic_score", MAC_HOMEOSTATIC), ("inflammatory_program_score", MAC_INFLAMMATORY)):
        gg = [g for g in genes if g in z.index and np.isfinite(z.loc[g, pre]).all()]
        meta[name] = z.loc[gg].mean(axis=0)
    # technical residual: score adjusted for log library size + mito fraction (label-free, all 92 samples)
    Xt = sm.add_constant(meta[["log_libsize", "pct_mito"]])
    tech = sm.OLS(meta["signature_score"], Xt).fit()
    meta["signature_score_tech_resid"] = tech.resid
    print(f"R2 of signature score ~ log library size + mito fraction (all samples): {tech.rsquared:.2f}")
    meta.to_csv(PROCESSED / "gse198520_sample_scores.tsv", sep="\t")

    b = meta.loc[pre].copy()
    rows = []

    def add(h, test, stat, p, n, note=""):
        rows.append(dict(hypothesis=h, test=test, statistic=stat, p=p, n=n, note=note))

    # H1 disease activity
    for c in ("das28 score", "inflammatory score", "crp"):
        ok = b[c].notna()
        r, p = spearmanr(b.loc[ok, "signature_score"], b.loc[ok, c])
        add("H1", f"Spearman score vs baseline {c}", r, p, int(ok.sum()))
    groups = [g["signature_score"].values for _, g in b.groupby("pathotype")]
    add("H1", "Kruskal score ~ pathotype", kruskal(*groups).statistic, kruskal(*groups).pvalue, len(b),
        "; ".join(f"{k} median {v:.2f}" for k, v in b.groupby("pathotype")["signature_score"].median().items()))

    # H2 response
    gn = b[b["response"].isin(["good", "non"])]
    y = (gn["response"] == "good").astype(int).values
    u = mannwhitneyu(gn.loc[gn.response == "good", "signature_score"], gn.loc[gn.response == "non", "signature_score"])
    add("H2", "Mann-Whitney score good vs non", u.statistic, u.pvalue, len(gn),
        f"median good {gn.loc[gn.response == 'good', 'signature_score'].median():.2f}, "
        f"non {gn.loc[gn.response == 'non', 'signature_score'].median():.2f}")
    r, p = spearmanr(b["signature_score"], b["response_ord"])
    add("H2", "Spearman score vs ordinal response (non<moderate<good)", r, p, len(b))
    r, p = spearmanr(b["signature_score"], b["delta das28"], nan_policy="omit")
    add("H2", "Spearman baseline score vs delta DAS28 (negative = improvement)", r, p, int(b["delta das28"].notna().sum()))

    auc_rows = []
    for name, col in (("signature score", "signature_score"), ("baseline DAS28", "das28 score"),
                      ("homeostatic program", "homeostatic_score"), ("inflammatory program", "inflammatory_program_score"),
                      ("log library size (technical)", "log_libsize"), ("mito fraction (technical)", "pct_mito"),
                      ("signature score, technical residual", "signature_score_tech_resid")):
        ok = gn[col].notna().values
        a, lo, hi, pp = auc_ci(y[ok], gn[col].values[ok], rng)
        auc_rows.append(dict(predictor=name, n=int(ok.sum()), auc_good_vs_non=a, ci_low=lo, ci_high=hi, perm_p=pp))
    auc_rows.append(dict(predictor="majority class", n=len(y), auc_good_vs_non=0.5, ci_low=np.nan, ci_high=np.nan, perm_p=np.nan))
    # random expression-matched gene sets
    pool = z.index[np.isfinite(z[pre]).all(axis=1) & (counts[pre].median(axis=1) >= 5)]
    expr = logcpm.loc[pool, pre].mean(axis=1)
    dec = pd.qcut(expr, 10, labels=False)
    sig_dec = [dec[g] if g in dec.index else int(pd.cut([logcpm.loc[g, pre].mean()], np.r_[-np.inf, expr.quantile(np.linspace(.1, .9, 9)).values, np.inf], labels=False)[0]) for g in sig]
    null = []
    for _ in range(N_RAND):
        pick = [rng.choice(dec.index[dec == q]) for q in sig_dec]
        null.append(roc_auc_score(y, z.loc[pick, gn.index].mean(axis=0).values))
    null = np.array(null)
    obs = auc_rows[0]["auc_good_vs_non"]
    auc_rows[0]["random_geneset_median_auc"] = np.median(null)
    auc_rows[0]["random_geneset_pct_ge_obs"] = float(np.mean(null >= obs))
    aucs = pd.DataFrame(auc_rows)
    aucs.round(4).to_csv(TAB / "gse198520_auc.tsv", sep="\t", index=False)
    print(aucs.round(3).to_string(index=False))

    # Confounders: RIN and biologic
    ok = b["rin"].notna()
    r, p = spearmanr(b.loc[ok, "signature_score"], b.loc[ok, "rin"])
    add("confound", "Spearman baseline score vs RIN", r, p, int(ok.sum()))
    for covs, label in ((["biologic"], "score + biologic"), (["biologic", "rin"], "score + biologic + RIN")):
        d = gn.dropna(subset=[c for c in covs if c == "rin"])
        X = pd.get_dummies(d[["signature_score"] + covs], columns=["biologic"], drop_first=True, dtype=float)
        X = sm.add_constant(X)
        fit = sm.Logit((d["response"] == "good").astype(int), X).fit(disp=0)
        add("H2", f"Logistic good vs non ~ {label}: OR per 1 SD of score",
            float(np.exp(fit.params["signature_score"] * d["signature_score"].std())), fit.pvalues["signature_score"], len(d))

    for v in ("log_libsize", "pct_mito"):
        r, p = spearmanr(b[v], b["response_ord"])
        add("confound", f"Spearman baseline {v} vs ordinal response", r, p, len(b))
        r, p = spearmanr(b[v], b["signature_score"])
        add("confound", f"Spearman baseline {v} vs signature score", r, p, len(b))
    d = gn.copy()
    fit = sm.Logit((d["response"] == "good").astype(int),
                   sm.add_constant(d[["signature_score", "log_libsize", "pct_mito"]])).fit(disp=0)
    add("H2-adjusted", "Logistic good vs non ~ score + log library size + mito: OR per 1 SD of score",
        float(np.exp(fit.params["signature_score"] * d["signature_score"].std())), fit.pvalues["signature_score"], len(d))
    r, p = spearmanr(b["signature_score_tech_resid"], b["response_ord"])
    add("H2-adjusted", "Spearman technical-residual score vs ordinal response", r, p, len(b))
    ok = b["inflammatory score"].notna()
    r, p = spearmanr(b.loc[ok, "inflammatory score"], b.loc[ok, "log_libsize"])
    add("confound", "Spearman histological inflammatory score vs log library size", r, p, int(ok.sum()),
        "library size may partly reflect tissue cellularity, so the technical adjustment is conservative")
    r, p = spearmanr(b.loc[ok, "signature_score_tech_resid"], b.loc[ok, "inflammatory score"])
    add("H1-adjusted", "Spearman technical-residual score vs histological inflammatory score", r, p, int(ok.sum()))
    # dominant expression axis at baseline (label-free PCA of well-measured genes)
    keep = (counts[pre] >= 10).mean(axis=1) >= 0.5
    Xp = logcpm.loc[keep, pre].T
    Xp = Xp - Xp.mean()
    U, S, _ = np.linalg.svd(Xp.values, full_matrices=False)
    b["PC1"] = U[:, 0] * S[0]
    var1 = S[0] ** 2 / np.sum(S ** 2)
    for v in ("log_libsize", "pct_mito", "rin", "response_ord", "signature_score"):
        ok = b[v].notna()
        r, p = spearmanr(b.loc[ok, "PC1"], b.loc[ok, v])
        add("confound", f"Spearman baseline PC1 ({var1:.0%} of variance) vs {v}", r, p, int(ok.sum()))
    # are the signature genes themselves library-size dependent?
    ok_g = z.index[np.isfinite(z[pre]).all(axis=1)]
    lib_pre = meta.loc[pre, "log_libsize"]
    grho = pd.Series({g: spearmanr(logcpm.loc[g, pre], lib_pre)[0] for g in ok_g})
    gbin = pd.qcut(logcpm.loc[ok_g, pre].mean(axis=1), 10, labels=False)
    obs_r = grho[sig].mean()
    null_r = np.array([grho[[rng.choice(gbin.index[gbin == gbin[g]]) for g in sig]].mean() for _ in range(2000)])
    add("confound", "Mean per-gene rho with log library size: signature vs expression-matched random sets",
        obs_r, (1 + np.sum(null_r >= obs_r)) / (1 + len(null_r)), len(sig), f"random median {np.median(null_r):.2f}")

    # H3 longitudinal
    wide = meta.pivot_table(index="patient", columns="timepoint", values="signature_score")
    resp = meta[meta["timepoint"] == "Pre"].set_index("patient")
    wide = wide.join(resp[["response", "delta das28", "biologic"]]).dropna(subset=["Pre", "Post"])
    wide["delta_score"] = wide["Post"] - wide["Pre"]
    w = wilcoxon(wide["Post"], wide["Pre"])
    add("H3", "Wilcoxon post vs pre (all)", w.statistic, w.pvalue, len(wide), f"median change {wide['delta_score'].median():.2f}")
    for g in ("good", "moderate", "non"):
        s = wide[wide["response"] == g]
        w = wilcoxon(s["Post"], s["Pre"])
        add("H3", f"Wilcoxon post vs pre ({g})", w.statistic, w.pvalue, len(s), f"median change {s['delta_score'].median():.2f}")
    u = mannwhitneyu(wide.loc[wide.response == "good", "delta_score"], wide.loc[wide.response == "non", "delta_score"])
    add("H3", "Mann-Whitney score change good vs non", u.statistic, u.pvalue,
        int(wide["response"].isin(["good", "non"]).sum()))
    r, p = spearmanr(wide["delta_score"], wide["delta das28"], nan_policy="omit")
    add("H3", "Spearman score change vs DAS28 change", r, p, int(wide["delta das28"].notna().sum()))
    lib_w = meta.pivot_table(index="patient", columns="timepoint", values="log_libsize")
    wide["delta_loglib"] = (lib_w["Post"] - lib_w["Pre"]).reindex(wide.index)
    for g in ("good", "non"):
        s_ = wide[wide["response"] == g]
        w = wilcoxon(s_["delta_loglib"])
        add("H3-adjusted", f"Wilcoxon log library size change ({g})", w.statistic, w.pvalue, len(s_),
            f"median change {s_['delta_loglib'].median():.2f}")
    res_w = meta.pivot_table(index="patient", columns="timepoint", values="signature_score_tech_resid")
    wide["delta_resid"] = (res_w["Post"] - res_w["Pre"]).reindex(wide.index)
    for g in ("good", "moderate", "non"):
        s_ = wide[wide["response"] == g]
        w = wilcoxon(s_["delta_resid"])
        add("H3-adjusted", f"Wilcoxon technical-residual score change ({g})", w.statistic, w.pvalue, len(s_),
            f"median change {s_['delta_resid'].median():.2f}")
    u = mannwhitneyu(wide.loc[wide.response == "good", "delta_resid"], wide.loc[wide.response == "non", "delta_resid"])
    add("H3-adjusted", "Mann-Whitney technical-residual change good vs non", u.statistic, u.pvalue,
        int(wide["response"].isin(["good", "non"]).sum()))

    tests = pd.DataFrame(rows)
    tests.round(4).to_csv(TAB / "gse198520_tests.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(tests.round(4).to_string(index=False))

    # Figure
    fig, axes = plt.subplots(1, 4, figsize=(22, 5))
    order = ["non", "moderate", "good"]
    ax = axes[0]
    for i, g in enumerate(order):
        v = b.loc[b.response == g, "signature_score"]
        ax.boxplot(v, positions=[i], widths=0.5, showfliers=False)
        ax.scatter(i + rng.uniform(-.12, .12, len(v)), v, s=22, zorder=3,
                   c=["C0" if x == "Certolizumab" else "C1" for x in b.loc[v.index, "biologic"]])
    t = tests.set_index("test")
    ax.set_xticks(range(3), ["EULAR non", "moderate", "good"])
    ax.set_ylabel("Baseline active-RA macrophage score")
    ax.set_title(f"Baseline score by response\ngood vs non p={t.loc['Mann-Whitney score good vs non', 'p']:.3f}", fontsize=10)
    ax = axes[1]
    ax.scatter(b["log_libsize"], b["signature_score"], s=22, c=[{"non": "C3", "moderate": "C7", "good": "C2"}[r] for r in b.response])
    r = t.loc["Spearman baseline log_libsize vs signature score"]
    ax.set_xlabel("log library size (red non, grey moderate, green good)"); ax.set_ylabel("Baseline score")
    ax.set_title(f"Technical axis: score vs library size rho={r.statistic:.2f}", fontsize=10)
    ax = axes[2]
    for name, col in (("signature score", "signature_score"), ("baseline DAS28", "das28 score"),
                      ("log library size (technical)", "log_libsize"),
                      ("signature score, technical residual", "signature_score_tech_resid")):
        ok = gn[col].notna().values
        fpr, tpr, _ = roc_curve(y[ok], gn[col].values[ok])
        ax.plot(fpr, tpr, label=f"{name} (AUC {aucs.set_index('predictor').loc[name, 'auc_good_vs_non']:.2f})")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8, label="majority class (0.50)")
    ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
    ax.set_title(f"Good vs non-responders (n={len(y)}); random gene sets median AUC {np.median(null):.2f}", fontsize=10)
    ax.legend(fontsize=8, loc="lower right")
    ax = axes[3]
    for g, col in (("non", "C3"), ("moderate", "C7"), ("good", "C2")):
        s = wide[wide.response == g]
        for _, rr in s.iterrows():
            ax.plot([0, 1], [rr["Pre"], rr["Post"]], color=col, alpha=0.35, lw=0.8)
        ax.plot([0, 1], [s["Pre"].median(), s["Post"].median()], color=col, lw=3, label=f"{g} (median)")
    ax.set_xticks([0, 1], ["Pre", "Week 12"]); ax.set_ylabel("Score")
    ax.set_title("H3: change with anti-TNF", fontsize=10)
    ax.legend(fontsize=8)
    fig.suptitle("GSE198520 anti-TNF cohort: active-RA macrophage signature (from E-MTAB-8322) "
                 "(blue = certolizumab, orange = etanercept)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "gse198520_response.png", dpi=150)


if __name__ == "__main__":
    main()
