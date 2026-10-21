"""Proteomics layer: is the IFN-gamma / active-RA macrophage programme visible at protein level in RA synovium?

Data: PRIDE PXD062696 (Wang et al. 2025 Proteomes 13:17; same tissue blocks as the Xenium/Visium data):
laser-capture microdissected regions, timsTOF DIA, DIA-NN gene-level matrix. Regions: healthy (NH) lining (SL) and
sublining (SSL), 3 donors; RA lining, fibrous sublining, immune sublining (4 of 6 donors), 6 donors; OA 1 donor.

Hypotheses fixed before outcome analysis (two-sided Mann-Whitney, donor as the unit):
  P1  Hallmark IFN-gamma response protein score: RA sublining (donor mean of immune + fibrous regions) vs NH sublining
      (6 vs 3); and RA lining vs NH lining (6 vs 3).
  P2  same for the active-RA macrophage signature proteins that are quantified.
  P3  within RA donors that have both: immune vs fibrous sublining (paired; descriptive, n = 4).
  Control: 1,000 random protein sets matched on abundance decile -> percentile of the observed RA - NH difference.

Input : data/raw/PXD062696/*, results/tables/signature_genes.tsv, Hallmark GMT
Output: data/processed/proteomics_region_scores.tsv, results/tables/proteomics_tests.tsv, proteomics_detected_genes.tsv
        results/figures/proteomics.png
"""
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

from config import RAW, PROCESSED, FIG, TAB, SEED, MAC_HOMEOSTATIC, MAC_INFLAMMATORY

P = RAW / "PXD062696"


def main():
    rng = np.random.default_rng(SEED)
    design = pd.read_excel(P / "experiment_design.xlsx").set_index("Proteomics Sample ID")
    gg = pd.read_csv(P / "report-first-pass.gg_matrix.tsv", sep="\t")
    gg = gg[gg["Genes"].notna() & (gg["Genes"] != "-")]
    gg["gene"] = gg["Genes"].str.split(";").str[0]
    gg = gg.drop(columns="Genes").groupby("gene").max()
    gg.columns = [re.search(r"(P\d+-[LSI])_", c).group(1) for c in gg.columns]
    meta = design.loc[gg.columns].copy()
    meta["region"] = np.where(meta["Tissue description"].str.contains("Lining"), "lining",
                     np.where(meta["Tissue description"].str.contains("Immune"), "immune_sublining",
                              np.where(meta["Tissue description"].str.contains("Fibrous"), "fibrous_sublining", "sublining")))
    meta["donor"] = meta["CT#"].astype(str) + "_" + meta["LabLion Block ID"].astype(str)
    use = meta.index[meta["Disease"].isin(["NH", "RA"])]
    X = np.log2(gg[use])
    X = X[X.notna().mean(axis=1) >= 0.7]
    X = X - X.median(axis=0)                                           # per-region median normalisation
    Z = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1), axis=0)
    print(f"Regions: {len(use)} (NH {int((meta.loc[use, 'Disease'] == 'NH').sum())}, RA {int((meta.loc[use, 'Disease'] == 'RA').sum())}); "
          f"proteins (genes) quantified in >= 70 % of regions: {len(X)}")

    sg = pd.read_csv(TAB / "signature_genes.tsv", sep="\t")
    sig = [g for g in sg.loc[(sg.signature == "primary_FDR") & (sg.direction == "up_in_active"), "gene"] if g in Z.index]
    hall = {}
    for line in open(next((RAW / "genesets").glob("h.all*.gmt"))):
        f = line.rstrip("\n").split("\t"); hall[f[0]] = f[2:]
    ifn = [g for g in hall["HALLMARK_INTERFERON_GAMMA_RESPONSE"] if g in Z.index and g != "IFNG"]
    sets = {"ifng_response": ifn, "ra_signature": sig,
            "mac_homeostatic": [g for g in MAC_HOMEOSTATIC if g in Z.index],
            "mac_inflammatory": [g for g in MAC_INFLAMMATORY if g in Z.index]}
    pd.DataFrame([(k, len(v), ",".join(v)) for k, v in sets.items()], columns=["set", "n_detected", "genes"]) \
        .to_csv(TAB / "proteomics_detected_genes.tsv", sep="\t", index=False)
    for k, v in sets.items():
        print(f"  {k}: {len(v)} proteins quantified" + (f" ({', '.join(v)})" if len(v) <= 15 else ""))

    sc = pd.DataFrame({k: Z.loc[v].mean(axis=0, skipna=True) for k, v in sets.items() if len(v) >= 2})
    sc = meta.loc[use, ["Disease", "region", "donor", "Sample group"]].join(sc)
    sc.round(4).to_csv(PROCESSED / "proteomics_region_scores.tsv", sep="\t")

    def donor_level(df, col, disease, regions):
        s = df[(df["Disease"] == disease) & df["region"].isin(regions)]
        return s.groupby("donor")[col].mean().dropna()

    expr_bins = pd.qcut(X.mean(axis=1), 10, labels=False)
    rows = []
    for col in [c for c in sets if c in sc.columns]:
        for label, ra_regions, nh_regions in (("sublining", ["immune_sublining", "fibrous_sublining", "sublining"], ["sublining"]),
                                              ("lining", ["lining"], ["lining"])):
            ra = donor_level(sc, col, "RA", ra_regions); nh = donor_level(sc, col, "NH", nh_regions)
            u = mannwhitneyu(ra, nh)
            obs = ra.median() - nh.median()
            # random protein-set control (same size, matched abundance decile)
            want = expr_bins.loc[sets[col]].values
            null = []
            for _ in range(1000):
                pick = [rng.choice(expr_bins.index[expr_bins == q]) for q in want]
                s2 = sc[["Disease", "region", "donor"]].assign(rnd=Z.loc[pick].mean(axis=0).reindex(sc.index))
                null.append(donor_level(s2, "rnd", "RA", ra_regions).median() - donor_level(s2, "rnd", "NH", nh_regions).median())
            null = np.array(null)
            rows.append(dict(set=col, comparison=f"RA vs NH {label} (donor means)", n_RA=len(ra), n_NH=len(nh),
                             median_RA=ra.median(), median_NH=nh.median(), diff=obs, p_mannwhitney=u.pvalue,
                             random_set_pct_ge_obs=float(np.mean(null >= obs))))
        r = sc[sc["Disease"] == "RA"].pivot_table(index="donor", columns="region", values=col)
        if {"immune_sublining", "fibrous_sublining"} <= set(r.columns):
            pr = r[["immune_sublining", "fibrous_sublining"]].dropna()
            d = pr["immune_sublining"] - pr["fibrous_sublining"]
            rows.append(dict(set=col, comparison="RA immune vs fibrous sublining (paired)", n_RA=len(pr), n_NH=0,
                             median_RA=pr["immune_sublining"].median(), median_NH=pr["fibrous_sublining"].median(),
                             diff=d.median(), p_mannwhitney=wilcoxon(d).pvalue if len(d) >= 3 else np.nan,
                             random_set_pct_ge_obs=np.nan, n_positive=int((d > 0).sum())))
    res = pd.DataFrame(rows)
    res.to_csv(TAB / "proteomics_tests.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 220)
    print(res.round(3).to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    order = [("NH", "lining"), ("RA", "lining"), ("NH", "sublining"), ("RA", "fibrous_sublining"), ("RA", "immune_sublining")]
    for ax, col, title in ((axes[0], "ifng_response", f"Hallmark IFN-gamma response ({len(sets['ifng_response'])} proteins)"),
                           (axes[1], "ra_signature" if "ra_signature" in sc.columns else "mac_inflammatory",
                            f"Active-RA macrophage signature ({len(sets['ra_signature'])} proteins)")):
        for i, (dz, rg) in enumerate(order):
            v = sc.loc[(sc.Disease == dz) & (sc.region.isin([rg] + (["fibrous_sublining"] if rg == "fibrous_sublining" else []))), col]
            if rg == "fibrous_sublining":
                v = sc.loc[(sc.Disease == dz) & sc.region.isin(["fibrous_sublining", "sublining"]), col]
            ax.boxplot(v, positions=[i], widths=0.5, showfliers=False)
            ax.scatter(i + rng.uniform(-.1, .1, len(v)), v, s=25, zorder=3, color="C3" if dz == "RA" else "C0")
        ax.set_xticks(range(len(order)), ["NH lining", "RA lining", "NH sublining", "RA fibrous\nsublining", "RA immune\nsublining"], fontsize=8)
        ax.set_ylabel("protein score (mean z)"); ax.set_title(title, fontsize=10)
    fig.suptitle("LCM spatial proteomics of synovium (PRIDE PXD062696)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "proteomics.png", dpi=150)


if __name__ == "__main__":
    main()
