"""Blood: is the synovial active-RA macrophage / IFN-gamma programme visible in circulating monocytes?

Cohort: CELLxGENE "Single-cell RNA-seq analysis reveals cell subsets and gene signatures associated with RA disease
activity" (JCI Insight 2024, doi:10.1172/jci.insight.178499): PBMC 10x 3' v3, 18 RA patients each paired with a
control (pair_index_CW), CDAI / DAS28-CRP, three sequencing batches.

Hypotheses fixed before outcome data were analysed (two-sided):
  H1  monocyte pseudobulk signature score higher in RA than in the matched control (paired Wilcoxon);
  H2  within RA, score correlates with DAS28-CRP and CDAI (Spearman);
  same tests for the Hallmark IFN-gamma response score (minus IFNG) and for the share of IFITM3+ (IFN-like) monocytes.
Scores: mean per-gene z of log2(CPM + 1) of donor monocyte pseudobulks (classical, non-classical, IL1B+, IFITM3+).
Batch is checked (Kruskal-Wallis of scores across batches).

Input : data/raw/cellxgene_blood/RA_disease_activity_blood.h5ad, results/tables/signature_genes.tsv, Hallmark GMT
Output: data/processed/blood_donor_scores.tsv, results/tables/blood_tests.tsv, results/figures/blood_validation.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import wilcoxon, spearmanr, kruskal, mannwhitneyu

from config import RAW, PROCESSED, FIG, TAB, SEED

MONO = ["Classical Monocytes", "Non-classical Monocytes", "IL1b-Monocytes", "IFITM3 Monocytes"]


def main():
    a = sc.read_h5ad(RAW / "cellxgene_blood" / "RA_disease_activity_blood.h5ad")
    sym = a.raw.var["feature_name"].astype(str).values if "feature_name" in a.raw.var else a.raw.var_names.values
    obs = a.obs
    mono = obs["fine_annot"].isin(MONO).values
    donors = sorted(obs["donor_id"].unique())
    C = a.raw.X
    rows, mats = [], []
    for d in donors:
        m = mono & (obs["donor_id"] == d).values
        if m.sum() < 100:
            continue
        v = np.asarray(C[m].sum(axis=0)).ravel()
        mats.append(v)
        o = obs.loc[obs["donor_id"] == d].iloc[0]
        rows.append(dict(donor=d, disease=o["disease"], pair=o["pair_index_CW"], batch=o["batch"],
                         das28crp=pd.to_numeric(o["MY_das28crp4"], errors="coerce"),
                         cdai=pd.to_numeric(o["MY_cdai"], errors="coerce"), mtx=o["MY_MTX"],
                         n_monocytes=int(m.sum()),
                         frac_IFITM3_mono=float((obs.loc[m, "fine_annot"] == "IFITM3 Monocytes").mean())))
    meta = pd.DataFrame(rows).set_index("donor")
    pb = pd.DataFrame(np.vstack(mats), index=meta.index, columns=sym)
    pb = pb.T.groupby(level=0).sum().T                                      # collapse duplicate symbols
    logcpm = np.log2(pb.div(pb.sum(axis=1), axis=0) * 1e6 + 1)
    expressed = logcpm.columns[(pb >= 10).mean(axis=0) >= 0.5]
    z = (logcpm[expressed] - logcpm[expressed].mean()) / logcpm[expressed].std()

    sg = pd.read_csv(TAB / "signature_genes.tsv", sep="\t")
    sig = [g for g in sg.loc[(sg.signature == "primary_FDR") & (sg.direction == "up_in_active"), "gene"] if g in z.columns]
    hall = {}
    for line in open(next((RAW / "genesets").glob("h.all*.gmt"))):
        f = line.rstrip("\n").split("\t"); hall[f[0]] = f[2:]
    ifn = [g for g in hall["HALLMARK_INTERFERON_GAMMA_RESPONSE"] if g in z.columns and g != "IFNG"]
    meta["signature_score"] = z[sig].mean(axis=1)
    meta["ifng_response_score"] = z[ifn].mean(axis=1)
    meta.to_csv(PROCESSED / "blood_donor_scores.tsv", sep="\t")
    print(f"Donors with >= 100 monocytes: {len(meta)} ({(meta.disease == 'rheumatoid arthritis').sum()} RA); "
          f"signature genes measured: {len(sig)}/15 ({', '.join(sig)}); IFN-gamma genes: {len(ifn)}")

    tests = []
    ra = meta[meta["disease"] == "rheumatoid arthritis"]
    ctl = meta[meta["disease"] == "normal"]
    for col in ("signature_score", "ifng_response_score", "frac_IFITM3_mono"):
        pairs = ra[["pair", col]].merge(ctl[["pair", col]], on="pair", suffixes=("_ra", "_ctl"))
        w = wilcoxon(pairs[f"{col}_ra"], pairs[f"{col}_ctl"])
        tests.append(dict(hypothesis="H1", test=f"paired Wilcoxon RA vs matched control: {col}", n=len(pairs),
                          statistic=float((pairs[f"{col}_ra"] - pairs[f"{col}_ctl"]).median()), p=w.pvalue,
                          note="statistic = median within-pair difference (RA - control)"))
        for clin in ("das28crp", "cdai"):
            ok = ra[clin].notna()
            r, p = spearmanr(ra.loc[ok, col], ra.loc[ok, clin])
            tests.append(dict(hypothesis="H2", test=f"Spearman {col} vs {clin} (RA)", n=int(ok.sum()), statistic=r, p=p))
        k = kruskal(*[g[col].values for _, g in meta.groupby("batch")])
        tests.append(dict(hypothesis="batch", test=f"Kruskal {col} across sequencing batches", n=len(meta),
                          statistic=k.statistic, p=k.pvalue))
    mtx = ra.dropna(subset=["mtx"])
    if mtx["mtx"].nunique() == 2:
        u = mannwhitneyu(mtx.loc[mtx.mtx == 1, "signature_score"], mtx.loc[mtx.mtx == 0, "signature_score"])
        tests.append(dict(hypothesis="confound", test="Mann-Whitney signature score, methotrexate yes vs no (RA)",
                          n=len(mtx), statistic=u.statistic, p=u.pvalue))
    tests = pd.DataFrame(tests)
    tests.round(4).to_csv(TAB / "blood_tests.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print(tests.round(3).to_string(index=False))

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for ax, col, title in ((axes[0], "signature_score", "Active-RA macrophage signature"),
                           (axes[1], "ifng_response_score", "Hallmark IFN-gamma response")):
        pairs = ra[["pair", col]].merge(ctl[["pair", col]], on="pair", suffixes=("_ra", "_ctl"))
        for _, r in pairs.iterrows():
            ax.plot([0, 1], [r[f"{col}_ctl"], r[f"{col}_ra"]], color="grey", alpha=0.5)
        ax.scatter(np.zeros(len(pairs)), pairs[f"{col}_ctl"], color="C0"); ax.scatter(np.ones(len(pairs)), pairs[f"{col}_ra"], color="C3")
        p = tests.set_index("test").loc[f"paired Wilcoxon RA vs matched control: {col}", "p"]
        ax.set_xticks([0, 1], ["matched control", "RA"]); ax.set_ylabel("monocyte pseudobulk score (z)")
        ax.set_title(f"{title}\npaired p = {p:.3f}", fontsize=10)
    ax = axes[2]
    ax.scatter(ra["das28crp"], ra["signature_score"], color="C3")
    r = tests.set_index("test").loc["Spearman signature_score vs das28crp (RA)"]
    ax.set_xlabel("DAS28-CRP"); ax.set_ylabel("monocyte signature score")
    ax.set_title(f"Within RA: rho = {r.statistic:.2f}, p = {r.p:.3f}", fontsize=10)
    fig.suptitle("Blood monocytes (CELLxGENE RA disease-activity PBMC cohort)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "blood_validation.png", dpi=150)


if __name__ == "__main__":
    main()
