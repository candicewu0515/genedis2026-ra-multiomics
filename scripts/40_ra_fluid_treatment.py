"""Second RA single-cell cohort: synovial fluid before and after adalimumab (anti-TNF) or tofacitinib (JAK inhibitor).

Cohort: GSE296117 (scRNA-seq, 10x), 9 RA patients, each sampled treatment-naive (BT) and on treatment (AT); 5 on
adalimumab, 4 on tofacitinib; author cell types used as deposited. Synovial fluid, not tissue; no response labels.
Programmes defined in RA synovial tissue and not re-fitted (scripts/state_markers.py): M6_programme (30 genes),
Tph_programme (30 genes); plus the 15-gene active-RA macrophage signature (step 03) and Hallmark IFN-gamma response
(IFNG removed).

Hypotheses fixed before the GSE296117 data were scored (two-sided):
  V1  an M-6-like state exists among treatment-naive fluid macrophages: the best-M6 cluster stands out more than the best
      cluster of 200 expression-matched random gene sets (same test as N1 in step 41);
  V2  primary: the patient-level macrophage M6 score falls on treatment (AT - BT; paired Wilcoxon, n = 9), and by more
      than the same change for 200 matched random gene sets; drugs reported separately (descriptive, n = 5 and 4);
      same tests for the active-RA signature and Hallmark IFN-gamma (secondary);
  V3  the share of M-6-like macrophages falls on treatment (cells above the 90th percentile of treatment-naive
      macrophage M6 scores; paired Wilcoxon);
  V4  the Tph score of fluid T cells changes on treatment (paired Wilcoxon; secondary).
Units: patient x time point with >= 50 cells of the type. Scores: per-gene z of log2(CPM + 1) across units, averaged.

Input : data/raw/gse296117/converted/ (from scripts/convert_gse296117.R), results/tables/signature_genes.tsv,
        results/tables/state_markers_gse243917.tsv, data/raw/genesets/h.all.v2024.1.Hs.symbols.gmt
Output: results/tables/ra_fluid_tests.tsv, ra_fluid_m6_cluster_markers.tsv; data/processed/ra_fluid_unit_scores.tsv;
        results/figures/ra_fluid_treatment.png
"""
import gzip

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
from scipy.stats import wilcoxon

from config import RAW, PROCESSED, FIG, TAB, SEED
from state_markers import markers, matched_sets, unit_scores, logcpm, cluster_best_state

D = RAW / "gse296117" / "converted"
N_RANDOM = 200
MIN_CELLS = 50


def load():
    X = scipy.io.mmread(D / "counts.mtx.gz").tocsr().T.tocsr()          # cells x genes
    genes = pd.read_csv(D / "genes.tsv.gz", header=None)[0].astype(str)
    meta = pd.read_csv(D / "meta.tsv.gz", sep="\t", low_memory=False)
    assert X.shape == (len(meta), len(genes))
    return X, genes, meta


def hallmark_ifng():
    for line in open(RAW / "genesets" / "h.all.v2024.1.Hs.symbols.gmt"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "HALLMARK_INTERFERON_GAMMA_RESPONSE":
            return [g for g in f[2:] if g != "IFNG"]


def pseudobulk(X, meta, genes, celltype):
    m = (meta.celltype == celltype).to_numpy()
    key = (meta.patient + "|" + meta.treatment)[m]
    n = key.value_counts()
    units = n[n >= MIN_CELLS].index
    Xm = X[m]
    pb = pd.DataFrame(np.vstack([np.asarray(Xm[(key == u).to_numpy()].sum(axis=0)).ravel() for u in units]),
                      index=units, columns=genes)
    pb = pb.T.groupby(level=0).sum().T
    info = meta[m].assign(unit=key.values).drop_duplicates("unit").set_index("unit").loc[units, ["patient", "treatment", "drug"]]
    info["n_cells"] = n.loc[units].values
    return pb, info


def paired_change(S, info, col):
    w = S[[col]].join(info).pivot(index="patient", columns="treatment", values=col).dropna()
    return (w["AT"] - w["BT"])


def main():
    rng = np.random.default_rng(SEED)
    mk = markers()
    m6 = mk.loc[mk.programme == "M6_programme", "gene"].tolist()
    tph = mk.loc[mk.programme == "Tph_programme", "gene"].tolist()
    sg = pd.read_csv(TAB / "signature_genes.tsv", sep="\t")
    act = sg.loc[(sg.signature == "primary_FDR") & (sg.direction == "up_in_active"), "gene"].tolist()
    ifng = hallmark_ifng()

    X, genes, meta = load()
    drug = meta.drop_duplicates("patient").set_index("patient")["drug"]
    print(f"GSE296117: {X.shape[0]:,} cells, {meta.patient.nunique()} patients; drugs: {drug.value_counts().to_dict()}")
    rows = []

    # V1 M-6-like state among treatment-naive macrophages
    sel = np.where(((meta.celltype == "Macrophage") & (meta.treatment == "BT")).to_numpy())[0]
    row, mkc, cellres = cluster_best_state(X[sel], meta.iloc[sel].reset_index(drop=True), genes, m6, "patient",
                                           "patient", rng, N_RANDOM, seed=SEED)
    rows.append(dict(test="V1 M-6-like cluster among treatment-naive fluid macrophages", **row))
    mkc.to_csv(TAB / "ra_fluid_m6_cluster_markers.tsv", sep="\t", index=False, float_format="%.4g")
    print(rows[-1])

    # V2 macrophage programmes, AT vs BT
    pb, info = pseudobulk(X, meta, genes, "Macrophage")
    L = logcpm(pb); L = L[L.columns[(L > 1).mean() >= 0.5]]
    rand = matched_sets([g for g in m6 if g in L.columns], L.mean(), N_RANDOM, rng)
    S = unit_scores(L, {"M6": m6, "active_RA_signature": act, "Hallmark_IFNG": ifng,
                        **{f"r{i}": r for i, r in enumerate(rand)}})
    null = np.array([np.median(paired_change(S, info, f"r{i}")) for i in range(N_RANDOM)])
    for col in ("M6", "active_RA_signature", "Hallmark_IFNG"):
        d = paired_change(S, info, col)
        p_rand = (1 + (np.abs(null - np.median(null)) >= abs(np.median(d) - np.median(null))).sum()) / (1 + N_RANDOM)
        rows.append(dict(test=f"V2 macrophage {col}: on treatment - naive (paired)", n=int(len(d)),
                         statistic=float(np.median(d)), null_median=float(np.median(null)), p=float(wilcoxon(d).pvalue),
                         p_vs_random=p_rand if col == "M6" else np.nan,
                         detail=f"fell in {(d < 0).sum()}/{len(d)}; adalimumab median {d[drug.loc[d.index] == 'Adalimumab'].median():.2f} "
                                f"(n={(drug.loc[d.index] == 'Adalimumab').sum()}), tofacitinib median "
                                f"{d[drug.loc[d.index] == 'Tofacitinib'].median():.2f} (n={(drug.loc[d.index] == 'Tofacitinib').sum()})"))

    # V3 share of M-6-like macrophages (cell level; threshold from treatment-naive macrophages)
    from state_markers import cell_scores
    mac = np.where((meta.celltype == "Macrophage").to_numpy())[0]
    cs = cell_scores(X[mac], genes, m6, rng)
    naive = (meta.treatment.iloc[mac] == "BT").to_numpy()
    thr = np.quantile(cs[naive], 0.90)
    share = pd.Series(cs > thr).groupby([meta.patient.iloc[mac].values, meta.treatment.iloc[mac].values]).mean().unstack()
    share = share.dropna()
    d3 = share["AT"] - share["BT"]
    rows.append(dict(test="V3 share of M-6-like macrophages: on treatment - naive (paired)", n=int(len(d3)),
                     statistic=float(np.median(d3)), p=float(wilcoxon(d3).pvalue),
                     detail=f"fell in {(d3 < 0).sum()}/{len(d3)}; naive median {share['BT'].median():.1%}, on treatment {share['AT'].median():.1%}"))

    # V4 Tph programme in T cells
    pbt, infot = pseudobulk(X, meta, genes, "T cell")
    Lt = logcpm(pbt); Lt = Lt[Lt.columns[(Lt > 1).mean() >= 0.5]]
    St = unit_scores(Lt, {"Tph": tph})
    d4 = paired_change(St, infot, "Tph")
    rows.append(dict(test="V4 T-cell Tph score: on treatment - naive (paired)", n=int(len(d4)),
                     statistic=float(np.median(d4)), p=float(wilcoxon(d4).pvalue), detail=f"fell in {(d4 < 0).sum()}/{len(d4)}"))

    tests = pd.DataFrame(rows)
    tests.to_csv(TAB / "ra_fluid_tests.tsv", sep="\t", index=False, float_format="%.4g")
    sh = share.stack().rename("m6_like_share")
    sh.index = [f"{p}|{t}" for p, t in sh.index]
    units = S[["M6", "active_RA_signature", "Hallmark_IFNG"]].join(info).join(
        St.rename(columns={"Tph": "tcell_Tph"}), how="outer").join(sh, how="left")
    units["patient"] = units.index.str.split("|").str[0]; units["treatment"] = units.index.str.split("|").str[1]
    units["drug"] = units["patient"].map(drug)
    units.round(4).to_csv(PROCESSED / "ra_fluid_unit_scores.tsv", sep="\t")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
    print(tests.round(4).to_string(index=False))

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    for k, col in enumerate(("M6", "Hallmark_IFNG")):
        w = S[[col]].join(info).pivot(index="patient", columns="treatment", values=col).dropna()
        for p, r in w.iterrows():
            ax[k].plot([0, 1], [r["BT"], r["AT"]], color="C3" if drug[p] == "Tofacitinib" else "C0", alpha=0.8)
        ax[k].set_xticks([0, 1], ["naive", "on treatment"]); ax[k].set_title(f"Fluid macrophages: {col}\n(blue adalimumab, red tofacitinib)", fontsize=9)
    ax[2].hist(null, bins=30, color="0.75"); ax[2].axvline(np.median(paired_change(S, info, "M6")), color="C3")
    ax[2].set_title("V2 median change: M6 vs 200 random sets", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "ra_fluid_treatment.png", dpi=150)


if __name__ == "__main__":
    main()
