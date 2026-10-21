"""Disease-stage validation in bulk synovium (GSE89408; Guo et al.; 218 biopsies, total RNA-seq):
healthy 28, arthralgia 10, undifferentiated arthritis 6, early RA 57, established RA 95, OA 22.

Normalisation: DESeq2 median-of-ratios size factors (CPM was dropped: random genes came out higher in RA than OA,
a composition effect of a few very abundant transcripts). Genes with CPM >= 1 in >= 20 % of biopsies.
Scores (mean z of log2 normalised counts across all biopsies):
  S1 primary   active-RA macrophage signature (15 FDR genes up in active RA, step 05)
  S2           a priori IFN macrophage programme (CXCL9/10/11, STAT1, GBP1/5, IDO1, IRF1, CD274, APOL6; as step 24)
  S3           homeostatic MERTK+ macrophage programme (expected lower in RA)
  macrophage   pan-macrophage abundance (CD68, CSF1R, C1QA, C1QB, CD14, CD163, AIF1), used as a covariate
Hypotheses fixed before analysis (OLS with age, sex and log library size; HC3 standard errors):
  H1 S1 higher in RA (early + established) than healthy
  H2 S1 higher in RA than OA (age-matched comparison; the cleaner test)
  H3 S1 rises with stage healthy (0) < arthralgia (1) < early RA (2)  (ordinal term)
  H4 H2 holds after adjusting for macrophage abundance (state, not just more macrophages)
  Specificity: every test is repeated on 1,000 random expression-matched gene sets of the same size (empirical p).
  Added after seeing S3 rise in RA vs healthy (healthy synovium holds far fewer macrophages): every test is also run
  with macrophage abundance as a covariate (H*m).
Exploratory: early vs established RA; ACPA+ vs ACPA- RA; healthy biopsies aged >= 40 only (the healthy group is younger).
Undifferentiated arthritis has no age recorded and enters the figure only.

Output: data/processed/gse89408_sample_scores.tsv, results/tables/gse89408_tests.tsv,
        results/figures/gse89408_disease_stage.png
"""
import gzip

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from config import RAW, PROCESSED, FIG, TAB, SEED, MAC_HOMEOSTATIC

D = RAW / "GSE89408"
IFN = ["CXCL9", "CXCL10", "CXCL11", "STAT1", "GBP1", "GBP5", "IDO1", "IRF1", "CD274", "APOL6"]
MAC = ["CD68", "CSF1R", "C1QA", "C1QB", "CD14", "CD163", "AIF1"]
PREFIX = {"healthy tissue": "normal_tissue", "arthralgia tissue": "AG_tissue", "undifferentiated arthritis tissue": "undiff_tissue",
          "osteoarthritis tissue": "OA_tissue", "rheumatoid arthritis tissue": "RA_tissue"}
GROUPS = {"Normal": "healthy", "Arthralgia": "arthralgia", "Undifferentiated arthritis": "UA",
          "Rheumatoid arthritis (early)": "early RA", "Rheumatoid arthritis (established)": "established RA",
          "Osteoarthritis": "OA"}
ORDER = ["healthy", "arthralgia", "UA", "early RA", "established RA", "OA"]


def metadata():
    rows = {}
    for line in gzip.open(D / "GSE89408_series_matrix.txt.gz", "rt"):
        if not line.startswith("!Sample_"):
            continue
        f = line.rstrip("\n").split("\t"); v = [x.strip('"') for x in f[1:]]
        if f[0] == "!Sample_characteristics_ch1":
            rows[v[0].split(":")[0]] = [x.split(": ", 1)[1] if ": " in x else np.nan for x in v]
        elif f[0] in ("!Sample_title", "!Sample_geo_accession"):
            rows[f[0][8:]] = v
    m = pd.DataFrame(rows)
    # count-matrix columns follow the GEO titles with a shorter prefix and the same number
    num = m["title"].str.extract(r"(\d+)$")[0]
    base = m["title"].str.replace(r" \d+$", "", regex=True)
    m["column"] = base.map(PREFIX) + "_" + num
    m["group"] = m["disease"].map(GROUPS)
    m["age"] = pd.to_numeric(m["age"], errors="coerce")
    m["female"] = (m["Sex"] == "F").astype(int)
    return m.set_index("column")


def score(z, genes):
    g = [x for x in genes if x in z.index]
    return z.loc[g].mean(axis=0), g


def main():
    rng = np.random.default_rng(SEED)
    meta = metadata()
    X = pd.read_csv(D / "GSE89408_GEO_count_matrix_rename.txt.gz", sep="\t", index_col=0)
    assert set(X.columns) == set(meta.index), "count columns do not match GEO titles"
    X = X[meta.index].round()
    X = X[~X.index.duplicated()]
    lib = X.sum(axis=0)
    cpm = X / lib * 1e6
    keep = (cpm >= 1).mean(axis=1) >= 0.2
    pos = X[(X > 0).all(axis=1)]
    lg = np.log(pos)
    sf = np.exp((lg.sub(lg.mean(axis=1), axis=0)).median(axis=0))       # median-of-ratios size factors
    L = np.log2(X[keep] / sf + 1)
    z = L.sub(L.mean(axis=1), axis=0).div(L.std(axis=1) + 1e-9, axis=0)
    print(f"Biopsies {X.shape[1]}; genes kept {keep.sum()}; size factors from {len(pos)} genes "
          f"(range {sf.min():.2f}-{sf.max():.2f}); groups: {meta.group.value_counts().to_dict()}")

    sig = pd.read_csv(TAB / "signature_genes.tsv", sep="\t")
    s1_genes = sig.loc[(sig.signature == "primary_FDR") & (sig.direction == "up_in_active"), "gene"].tolist()
    meta["S1"], g1 = score(z, s1_genes)
    meta["S2"], g2 = score(z, IFN)
    meta["S3"], g3 = score(z, MAC_HOMEOSTATIC)
    meta["macrophage"], gm = score(z, MAC)
    meta["log_lib"] = np.log(lib[meta.index])
    print(f"S1 genes found {len(g1)}/{len(s1_genes)}; S2 {len(g2)}/{len(IFN)}; S3 {len(g3)}/{len(MAC_HOMEOSTATIC)}; macrophage {len(gm)}")
    meta.drop(columns=["title"]).round(4).to_csv(PROCESSED / "gse89408_sample_scores.tsv", sep="\t")
    print(meta.groupby("group")[["S1", "S2", "S3", "macrophage", "age"]].median().loc[ORDER].round(2).to_string())

    def fit(df, formula, term):
        f = smf.ols(formula, data=df).fit(cov_type="HC3")
        return f.params[term], f.bse[term], f.tvalues[term], f.pvalues[term], int(f.nobs)

    adj = "age + female + log_lib"
    ra = meta.group.isin(["early RA", "established RA"])
    d1 = meta[ra | (meta.group == "healthy")].assign(RA=lambda d: d.group.ne("healthy").astype(int)).dropna(subset=["age"])
    d2 = meta[ra | (meta.group == "OA")].assign(RA=lambda d: d.group.ne("OA").astype(int)).dropna(subset=["age"])
    d3 = meta[meta.group.isin(["healthy", "arthralgia", "early RA"])].assign(
        stage=lambda d: d.group.map({"healthy": 0, "arthralgia": 1, "early RA": 2})).dropna(subset=["age"])
    d5 = meta[ra].assign(early=lambda d: (d.group == "early RA").astype(int)).dropna(subset=["age"])
    d6 = meta[ra & meta["acpa positivity"].isin(["Positive", "Negative"])].assign(
        acpa=lambda d: (d["acpa positivity"] == "Positive").astype(int)).dropna(subset=["age"])
    d7 = d1[(d1.group != "healthy") | (d1.age >= 40)]
    tests = [("H1 RA vs healthy", d1, "RA", ""), ("H2 RA vs OA", d2, "RA", ""),
             ("H3 stage trend healthy < arthralgia < early RA", d3, "stage", ""),
             ("H4 RA vs OA, adjusted for macrophage abundance", d2, "RA", " + macrophage"),
             ("H1m RA vs healthy, adjusted for macrophage abundance", d1, "RA", " + macrophage"),
             ("H3m stage trend, adjusted for macrophage abundance", d3, "stage", " + macrophage"),
             ("E1 early vs established RA", d5, "early", ""), ("E2 ACPA+ vs ACPA- RA", d6, "acpa", ""),
             ("E3 RA vs healthy aged >= 40", d7, "RA", "")]
    # random expression-matched gene sets (same draws for every test of a score)
    mean_expr = L.mean(axis=1)
    bins = pd.qcut(mean_expr.rank(method="first"), 20, labels=False)
    members = {q: bins.index[bins == q].values for q in range(20)}
    rows = []
    for s, genes in (("S1", g1), ("S2", g2), ("S3", g3)):
        want = bins.loc[genes].values
        rand = np.vstack([z.loc[[rng.choice(members[q]) for q in want]].mean(axis=0).values for _ in range(1000)])
        rand = pd.DataFrame(rand.T, index=z.columns)
        for test, df, term, extra in tests:
            if test.startswith("E") and s == "S3":
                continue
            b, se, t, p, n = fit(df, f"{s} ~ {term} + {adj}{extra}", term)
            tn = np.array([fit(df.assign(R=rand.loc[df.index, k].values), f"R ~ {term} + {adj}{extra}", term)[2]
                           for k in range(rand.shape[1])])
            side = np.sign(t) if t != 0 else 1
            rows.append(dict(score=s, test=test, n=n, beta=b, se=se, t=t, p=p, random_median_t=float(np.median(tn)),
                             p_vs_random=(1 + (side * tn >= side * t).sum()) / (1 + len(tn))))
    res = pd.DataFrame(rows)
    res.to_csv(TAB / "gse89408_tests.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 200)
    print(res.round(4).to_string(index=False))
    r_age = np.corrcoef(meta.loc[meta.group == "healthy", "S1"], meta.loc[meta.group == "healthy", "age"])[0, 1]
    print(f"Healthy biopsies: S1 vs age r = {r_age:.2f}")

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    for ax, s, title in zip(axes, ("S1", "S2", "S3"), ("S1 active-RA macrophage signature (primary)", "S2 IFN macrophage programme",
                                                         "S3 homeostatic MERTK+ programme")):
        for i, g in enumerate(ORDER):
            v = meta.loc[meta.group == g, s]
            ax.boxplot(v, positions=[i], widths=0.55, showfliers=False)
            ax.scatter(i + rng.uniform(-.15, .15, len(v)), v, s=8, alpha=0.6, color="C3" if "RA" in g else "C0")
        ax.set_xticks(range(len(ORDER)), [f"{g}\n(n={int((meta.group == g).sum())})" for g in ORDER], fontsize=8)
        t = res[res.score == s].set_index("test")
        ax.set_title(f"{title}\nRA vs healthy p={t.loc['H1 RA vs healthy', 'p']:.1g}; RA vs OA p={t.loc['H2 RA vs OA', 'p']:.1g}; "
                     f"vs OA beyond random genes p={t.loc['H2 RA vs OA', 'p_vs_random']:.2g}", fontsize=9)
        ax.set_ylabel("score (mean z)")
    fig.suptitle("GSE89408 synovial biopsies by disease stage (OLS adjusted for age, sex, library size)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "gse89408_disease_stage.png", dpi=150)


if __name__ == "__main__":
    main()
