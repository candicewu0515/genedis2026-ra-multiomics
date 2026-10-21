"""Bulk deconvolution of GSE198520 (anti-TNF cohort) with a GSE243917 single-nucleus reference.

Reference states: broad cell types plus myeloid split into M-6 IFN-activated (STAT1+ CXCL10+), SPP1+ (M-4),
MERTK+ resident (M-0..M-3), monocyte-derived IL1B+ FCN1+ (M-7), other macrophages and DCs.
Signature: mean CPM per state; top 50 marker genes per state by log fold change vs the other states, kept if detected
in the bulk panel; each gene scaled by its reference mean (applied to bulk too); non-negative least squares per bulk
sample, proportions rescaled to sum 1. (A first version without gene scaling failed the sanity checks: one state
absorbed ~50 % and fibroblasts were ~3 %.)

Checks and tests (H fixed before outcome analysis):
  sanity  estimated lymphoid share (T + B + plasma) higher in lymphoid than fibroid pathotype; correlates with
          histological inflammatory score
  D1      baseline M-6 share, good vs non-responders (Mann-Whitney; logistic adjusted for log library size + mito)
  D2      change of M-6 share week 12 vs baseline, good vs non-responders
Caveat: nuclear (snRNA) reference vs whole-cell bulk RNA; M-6 is defined by 65 nuclei.

Output: data/processed/gse198520_deconvolution.tsv, results/tables/deconvolution_tests.tsv,
        results/figures/deconvolution.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm
from scipy.optimize import nnls
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

from config import RAW, INTERIM, PROCESSED, FIG, TAB, SEED

N_MARK = 50


def state_of(ct, sub):
    code = sub.split(":")[0]
    if ct != "myeloid":
        return ct
    if code == "M-6":
        return "mac_M6_IFN"
    if code == "M-4":
        return "mac_SPP1"
    if code in ("M-0", "M-1", "M-2", "M-3"):
        return "mac_MERTK_resident"
    if code == "M-7":
        return "mono_derived_IL1B"
    if code in ("M-9", "M-10", "M-11", "M-12", "M-13", "M-14"):
        return "DC"
    return "mac_other"


def main():
    rna = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    states = np.array([state_of(c, s) for c, s in zip(rna.obs["author_ct"].astype(str), rna.obs["author_subtype"].astype(str))])
    C = rna.layers["counts"].tocsr()
    lib = np.asarray(C.sum(axis=1)).ravel()
    cpm_cells = C.multiply(1e6 / lib[:, None]).tocsr()
    order = sorted(set(states))
    ref = pd.DataFrame({s: np.asarray(cpm_cells[states == s].mean(axis=0)).ravel() for s in order}, index=rna.var_names)
    ref = ref[~ref.index.str.startswith("MT-")]

    bulk = pd.read_csv(RAW / "GSE198520" / "GSE198520_Raw_gene_count_matrix.txt.gz", sep="\t", index_col=0)
    bulk = bulk[~bulk.index.str.startswith("MT-")]
    common = ref.index.intersection(bulk.index)
    ref, bulk = ref.loc[common], bulk.loc[common]
    lref = np.log2(ref + 1)
    markers = set()
    for s in order:
        lfc = lref[s] - lref.drop(columns=s).max(axis=1)
        markers |= set(lfc[(ref[s] > 10)].sort_values(ascending=False).index[:N_MARK])
    markers = sorted(markers)
    # keep markers measured in the bulk panel; scale each gene by its mean in the reference so that highly expressed
    # genes do not dominate the fit (CIBERSORT-style standardisation; same factor applied to bulk)
    detected = bulk.index[(bulk >= 5).mean(axis=1) >= 0.5]
    markers = [g for g in markers if g in detected]
    S = ref.loc[markers]
    B = bulk.loc[markers] / bulk.sum(axis=0) * 1e6
    scale = S.mean(axis=1)
    S = S.div(scale, axis=0)
    B = B.div(scale, axis=0)
    print(f"Reference states: {[str(o) for o in order]}; marker genes used (detected in bulk): {len(markers)}")

    props = {}
    for smp in B.columns:
        w, _ = nnls(S.values, B[smp].values)
        props[smp] = w / w.sum() if w.sum() > 0 else w
    P = pd.DataFrame(props, index=order).T
    meta = pd.read_csv(PROCESSED / "gse198520_sample_scores.tsv", sep="\t", index_col=0)
    P = P.join(meta[["patient", "timepoint", "response", "response_ord", "pathotype", "inflammatory score", "biologic",
                     "log_libsize", "pct_mito"]])
    P["lymphoid_share"] = P[[c for c in ("Tcell", "Bcell", "plasma", "NK") if c in P]].sum(axis=1)
    P.round(5).to_csv(PROCESSED / "gse198520_deconvolution.tsv", sep="\t")

    rows = []
    b = P[P.timepoint == "Pre"]
    lf = b[b.pathotype.isin(["Lymphoid", "Fibroid"])]
    u = mannwhitneyu(lf.loc[lf.pathotype == "Lymphoid", "lymphoid_share"], lf.loc[lf.pathotype == "Fibroid", "lymphoid_share"])
    rows.append(dict(test="sanity: lymphoid share, lymphoid vs fibroid pathotype", n=len(lf),
                     statistic=lf.loc[lf.pathotype == "Lymphoid", "lymphoid_share"].median() - lf.loc[lf.pathotype == "Fibroid", "lymphoid_share"].median(),
                     p=u.pvalue))
    ok = b["inflammatory score"].notna()
    r, p = spearmanr(b.loc[ok, "lymphoid_share"], b.loc[ok, "inflammatory score"])
    rows.append(dict(test="sanity: lymphoid share vs histological inflammatory score", n=int(ok.sum()), statistic=r, p=p))
    gn = b[b.response.isin(["good", "non"])]
    for st in ("mac_M6_IFN", "mac_SPP1", "mac_MERTK_resident", "fibroblast", "lymphoid_share"):
        u = mannwhitneyu(gn.loc[gn.response == "good", st], gn.loc[gn.response == "non", st])
        r_lib, _ = spearmanr(b[st], b["log_libsize"])
        X = sm.add_constant(pd.DataFrame({"z": (gn[st] - gn[st].mean()) / (gn[st].std() + 1e-12),
                                          "log_libsize": gn.log_libsize, "pct_mito": gn.pct_mito}))
        try:
            fit = sm.Logit((gn.response == "good").astype(int), X).fit(disp=0)
            or_adj, p_adj = float(np.exp(fit.params["z"])), float(fit.pvalues["z"])
        except Exception:
            or_adj, p_adj = np.nan, np.nan
        rows.append(dict(test=f"D1 baseline {st}: good vs non", n=len(gn),
                         statistic=gn.loc[gn.response == "good", st].median() - gn.loc[gn.response == "non", st].median(),
                         p=u.pvalue, rho_with_log_libsize=r_lib, OR_adjusted_per_SD=or_adj, p_adjusted=p_adj))
    wide = P.pivot_table(index="patient", columns="timepoint", values="mac_M6_IFN")
    resp = P[P.timepoint == "Pre"].set_index("patient")["response"]
    wide = wide.join(resp).dropna(subset=["Pre", "Post"])
    wide["delta"] = wide["Post"] - wide["Pre"]
    for g in ("good", "non"):
        s = wide[wide.response == g]
        rows.append(dict(test=f"D2 M-6 share change week 12 - baseline ({g})", n=len(s), statistic=s.delta.median(),
                         p=wilcoxon(s.delta).pvalue if (s.delta != 0).any() else 1.0))
    u = mannwhitneyu(wide.loc[wide.response == "good", "delta"], wide.loc[wide.response == "non", "delta"])
    rows.append(dict(test="D2 M-6 share change, good vs non", n=int(wide.response.isin(["good", "non"]).sum()),
                     statistic=np.nan, p=u.pvalue))
    res = pd.DataFrame(rows)
    res.to_csv(TAB / "deconvolution_tests.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 220)
    print(res.round(4).to_string(index=False))
    print("\nMedian baseline proportions:"); print(b[order].median().round(3).to_string())

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    ax = axes[0]
    b.sort_values("lymphoid_share")[order].plot.bar(stacked=True, ax=ax, width=0.9, legend=False, colormap="tab20")
    ax.set_xticks([]); ax.set_ylabel("estimated proportion"); ax.set_title("GSE198520 baseline biopsies (NNLS)", fontsize=10)
    ax.legend(order, fontsize=6, bbox_to_anchor=(1, 1))
    ax = axes[1]
    for i, g in enumerate(["non", "moderate", "good"]):
        v = b.loc[b.response == g, "mac_M6_IFN"]
        ax.boxplot(v, positions=[i], widths=0.5, showfliers=False)
        ax.scatter(i + np.random.default_rng(SEED).uniform(-.1, .1, len(v)), v, s=18, zorder=3)
    t = res.set_index("test")
    ax.set_xticks(range(3), ["non", "moderate", "good"]); ax.set_ylabel("estimated M-6 share")
    ax.set_title(f"Baseline M-6 share by EULAR response\ngood vs non p = {t.loc['D1 baseline mac_M6_IFN: good vs non', 'p']:.3f}", fontsize=10)
    ax = axes[2]
    for g, col in (("non", "C3"), ("good", "C2")):
        s = wide[wide.response == g]
        for _, r in s.iterrows():
            ax.plot([0, 1], [r.Pre, r.Post], color=col, alpha=0.4)
    ax.set_xticks([0, 1], ["baseline", "week 12"]); ax.set_ylabel("estimated M-6 share")
    ax.set_title("M-6 share with anti-TNF (green good, red non)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "deconvolution.png", dpi=150)


if __name__ == "__main__":
    main()
