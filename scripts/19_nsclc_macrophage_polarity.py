"""Cross-indication: the synovial CXCL9 (IFN) vs SPP1 macrophage axis in NSCLC (TCGA LUAD / LUSC).

Rationale: in RA synovium the IFN-activated M-6 state is CXCL9/CXCL10-high and the M-4 state is SPP1-high (checked
below in GSE243917). In tumours a CXCL9:SPP1 macrophage polarity has been linked to prognosis. We test the same axis,
and the RA-derived active-macrophage signature, in TCGA lung adenocarcinoma and squamous carcinoma.

Hypotheses fixed before survival data were analysed:
  H1  primary: polarity = z(CXCL9) - z(SPP1) in primary tumours; higher polarity -> longer overall survival
      (Cox PH per SD, adjusted for age and stage I-II vs III-IV), tested in LUAD and LUSC separately (two-sided p).
  H2  secondary: RA active-macrophage signature score (same model).
  H3  descriptive: polarity in tumour vs adjacent normal (paired Wilcoxon where both exist).
Data: UCSC Xena TCGA HiSeqV2 (log2(RSEM + 1)); TCGA-CDR curated survival (Liu et al. 2018 Cell).

Output: results/tables/nsclc_ra_axis_check.tsv, nsclc_cox.tsv, nsclc_tumour_vs_normal.tsv,
        data/processed/nsclc_sample_scores.tsv, results/figures/nsclc_polarity.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import wilcoxon
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import SurvfuncRight, survdiff

from config import RAW, INTERIM, PROCESSED, FIG, TAB

T = RAW / "tcga"


def ra_axis_check():
    rna = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    my = rna[rna.obs["author_ct"] == "myeloid"]
    df = sc.get.obs_df(my, ["CXCL9", "CXCL10", "SPP1"], use_raw=True)
    df["subtype"] = my.obs["author_subtype"].astype(str).values
    tab = df.groupby("subtype").agg(n=("CXCL9", "size"), CXCL9=("CXCL9", "mean"), CXCL10=("CXCL10", "mean"), SPP1=("SPP1", "mean"))
    tab = tab[tab.n >= 30].sort_values("CXCL9", ascending=False)
    tab.round(3).to_csv(TAB / "nsclc_ra_axis_check.tsv", sep="\t")
    print("RA synovium (GSE243917) myeloid subtypes, mean log-normalised expression:")
    print(tab.round(2).to_string())


def load(cancer):
    x = pd.read_csv(T / f"{cancer}_HiSeqV2.gz", sep="\t", index_col=0)
    clin = pd.read_csv(T / f"{cancer}_clinicalMatrix.tsv", sep="\t", index_col=0, low_memory=False)
    return x, clin


def stage_group(s):
    s = str(s).upper()
    if "IV" in s or "III" in s:
        return 1
    if "I" in s:
        return 0
    return np.nan


def main():
    ra_axis_check()
    surv = pd.read_csv(T / "Survival_SupplementalTable_S1_20171025_xena_sp.tsv", sep="\t", index_col=0)
    sg = pd.read_csv(TAB / "signature_genes.tsv", sep="\t")
    sig = sg.loc[(sg.signature == "primary_FDR") & (sg.direction == "up_in_active"), "gene"].tolist()

    cox_rows, tn_rows, all_scores = [], [], []
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    for k, cancer in enumerate(("LUAD", "LUSC")):
        x, clin = load(cancer)
        stype = x.columns.str[13:15]
        tum = x.columns[stype == "01"]
        z = x.sub(x[tum].mean(axis=1), axis=0).div(x[tum].std(axis=1).replace(0, np.nan), axis=0)
        sc_df = pd.DataFrame(index=x.columns)
        sc_df["polarity_CXCL9_minus_SPP1"] = z.loc["CXCL9"] - z.loc["SPP1"]
        g = [s for s in sig if s in z.index]
        sc_df["ra_signature"] = z.loc[g].mean(axis=0)
        sc_df["cancer"] = cancer
        sc_df["sample_type"] = stype
        all_scores.append(sc_df)

        # H3 tumour vs adjacent normal (paired by patient)
        norm = x.columns[stype == "11"]
        pairs = [(t, n) for n in norm for t in tum if t[:12] == n[:12]]
        if pairs:
            d = np.array([sc_df.loc[t, "polarity_CXCL9_minus_SPP1"] - sc_df.loc[n, "polarity_CXCL9_minus_SPP1"] for t, n in pairs])
            tn_rows.append(dict(cancer=cancer, n_pairs=len(pairs), median_tumour_minus_normal=np.median(d),
                                p_wilcoxon=wilcoxon(d).pvalue))

        # H1/H2 survival in primary tumours (one sample per patient)
        t = sc_df.loc[tum].copy()
        t["patient"] = t.index.str[:15]
        t = t.join(surv[["OS", "OS.time"]], on="patient")
        t = t.join(clin[["age_at_initial_pathologic_diagnosis", "pathologic_stage"]], how="left")
        t["age"] = pd.to_numeric(t["age_at_initial_pathologic_diagnosis"], errors="coerce")
        t["stage_III_IV"] = t["pathologic_stage"].map(stage_group)
        t = t.dropna(subset=["OS", "OS.time", "age", "stage_III_IV"])
        t = t[t["OS.time"] > 0]
        t = t[~t.index.str[:12].duplicated()]                      # one primary tumour per patient
        for score in ("polarity_CXCL9_minus_SPP1", "ra_signature"):
            X = pd.DataFrame({"score_per_SD": (t[score] - t[score].mean()) / t[score].std(),
                              "age_per_10y": t["age"] / 10, "stage_III_IV": t["stage_III_IV"]})
            res = PHReg(t["OS.time"].values, X.values, status=t["OS"].values.astype(int)).fit()
            hr = np.exp(res.params); ci = np.exp(res.conf_int())
            hi = t[score] >= t[score].median()
            lr = survdiff(t["OS.time"].values, t["OS"].values.astype(int), hi.astype(int).values)
            cox_rows.append(dict(cancer=cancer, score=score, n=len(t), events=int(t["OS"].sum()),
                                 HR_per_SD=hr[0], CI_low=ci[0, 0], CI_high=ci[0, 1], p=res.pvalues[0],
                                 HR_stage=hr[2], p_stage=res.pvalues[2], logrank_median_split_p=lr[1]))
            if score == "polarity_CXCL9_minus_SPP1":
                ax = axes[k]
                for grp, lab, col in ((True, "CXCL9-high polarity", "C3"), (False, "SPP1-high polarity", "C0")):
                    sf = SurvfuncRight(t.loc[hi == grp, "OS.time"].values / 365.25, t.loc[hi == grp, "OS"].values.astype(int))
                    ax.step(sf.surv_times, sf.surv_prob, where="post", color=col, label=f"{lab} (n={int((hi == grp).sum())})")
                ax.set_xlabel("Years"); ax.set_ylabel("Overall survival")
                ax.set_title(f"TCGA {cancer}: CXCL9 - SPP1 polarity (median split)\n"
                             f"adjusted HR per SD = {hr[0]:.2f} ({ci[0, 0]:.2f}-{ci[0, 1]:.2f}), p = {res.pvalues[0]:.3g}", fontsize=9)
                ax.legend(fontsize=8); ax.set_xlim(0, 10)

    scores = pd.concat(all_scores)
    scores.round(4).to_csv(PROCESSED / "nsclc_sample_scores.tsv", sep="\t")
    cox = pd.DataFrame(cox_rows)
    cox.to_csv(TAB / "nsclc_cox.tsv", sep="\t", index=False, float_format="%.4g")
    tn = pd.DataFrame(tn_rows)
    tn.to_csv(TAB / "nsclc_tumour_vs_normal.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 220)
    print("\nCox PH (overall survival; adjusted for age and stage):"); print(cox.round(3).to_string(index=False))
    print("\nTumour vs adjacent normal, polarity:"); print(tn.round(4).to_string(index=False))

    ax = axes[2]
    s = scores[scores.sample_type.isin(["01", "11"])]
    data = [s.loc[(s.cancer == c) & (s.sample_type == st), "polarity_CXCL9_minus_SPP1"] for c in ("LUAD", "LUSC") for st in ("11", "01")]
    ax.boxplot(data, showfliers=False)
    ax.set_xticks(range(1, 5), ["LUAD normal", "LUAD tumour", "LUSC normal", "LUSC tumour"], fontsize=8)
    ax.set_ylabel("CXCL9 - SPP1 polarity (z, tumour-referenced)")
    ax.set_title("Tumour vs adjacent normal lung", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "nsclc_polarity.png", dpi=150)


if __name__ == "__main__":
    main()
