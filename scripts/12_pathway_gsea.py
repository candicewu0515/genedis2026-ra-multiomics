"""Pathway layer: pre-ranked GSEA (MSigDB Hallmark + Reactome).

Contrast A (primary): active RA vs remission, synovial macrophages, E-MTAB-8322 pseudobulk DESeq2 (step 04),
          ranked by Wald statistic.
Contrast B (cross-check, independent cohort): week 12 vs baseline in EULAR good responders to anti-TNF, GSE198520,
          paired PyDESeq2 (~ patient + timepoint; library size did not change in this group, step 11).
Question: are pathways up in active-RA macrophages the ones that fall after effective anti-TNF?

Input : results/tables/emtab8322_DE_active_vs_remission.tsv, data/raw/GSE198520/*, data/processed/gse198520_sample_scores.tsv
Output: data/raw/genesets/*.gmt (MSigDB download, versioned), results/tables/gsea_*.tsv, gse198520_DE_post_vs_pre_good.tsv
        results/figures/gsea_hallmark.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gseapy as gp
from scipy.stats import spearmanr

from config import RAW, PROCESSED, FIG, TAB, SEED

GS = RAW / "genesets"
GS.mkdir(exist_ok=True)
MSIGDB_VER = "2024.1.Hs"
LIBS = {"hallmark": "h.all", "reactome": "c2.cp.reactome"}


def load_gmt(name, cat):
    path = GS / f"{cat}.v{MSIGDB_VER}.symbols.gmt"
    if not path.exists():
        gmt = gp.Msigdb().get_gmt(category=cat, dbver=MSIGDB_VER)
        with open(path, "w") as fh:
            for k, genes in gmt.items():
                fh.write("\t".join([k, "msigdb"] + list(genes)) + "\n")
    sets = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        sets[f[0]] = f[2:]
    return sets


def prerank(rnk, sets, tag):
    res = gp.prerank(rnk=rnk, gene_sets=sets, min_size=15, max_size=500, permutation_num=1000,
                     seed=SEED, threads=8, outdir=None, verbose=False).res2d
    res = res.rename(columns={"Term": "pathway", "NES": "nes", "NOM p-val": "p", "FDR q-val": "fdr"})
    res["nes"] = res["nes"].astype(float); res["p"] = res["p"].astype(float); res["fdr"] = res["fdr"].astype(float)
    res["contrast"] = tag
    return res[["contrast", "pathway", "nes", "p", "fdr", "Lead_genes"]].sort_values("nes", ascending=False)


def paired_de_good_responders():
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.default_inference import DefaultInference
    meta = pd.read_csv(PROCESSED / "gse198520_sample_scores.tsv", sep="\t", index_col=0)
    meta = meta[meta["response"] == "good"]
    counts = pd.read_csv(RAW / "GSE198520" / "GSE198520_Raw_gene_count_matrix.txt.gz", sep="\t", index_col=0)
    counts = counts.loc[~counts.index.str.startswith("MT-"), meta.index].T
    counts = counts.loc[:, (counts >= 10).sum(axis=0) >= 0.5 * len(counts)]
    md = meta[["patient", "timepoint"]].copy()
    md["timepoint"] = md["timepoint"].map({"Pre": "pre", "Post": "post"})
    inf = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(counts=counts.astype(int), metadata=md, design="~patient + timepoint",
                       refit_cooks=True, inference=inf, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["timepoint", "post", "pre"], inference=inf, quiet=True)
    st.summary()
    de = st.results_df.sort_values("pvalue")
    de.round(5).to_csv(TAB / "gse198520_DE_post_vs_pre_good.tsv", sep="\t")
    print(f"GSE198520 good responders post vs pre (n = {md.patient.nunique()} pairs): {(de.padj < 0.05).sum()} genes FDR < 0.05")
    return de


def main():
    de_a = pd.read_csv(TAB / "emtab8322_DE_active_vs_remission.tsv", sep="\t", index_col=0).dropna(subset=["stat"])
    de_a = de_a[~de_a.index.str.match(r"^(ENSG|MT-)")]
    rnk_a = de_a["stat"].sort_values(ascending=False)
    de_b = paired_de_good_responders().dropna(subset=["stat"])
    rnk_b = de_b["stat"].sort_values(ascending=False)

    all_res = []
    for lib, cat in LIBS.items():
        sets = load_gmt(lib, cat)
        for tag, rnk in (("A_active_vs_remission_macrophages", rnk_a), ("B_antiTNF_good_post_vs_pre", rnk_b)):
            r = prerank(rnk.reset_index().rename(columns={rnk.index.name or "index": "gene"}), sets, tag)
            r["library"] = lib
            all_res.append(r)
    res = pd.concat(all_res)
    res.to_csv(TAB / "gsea_results.tsv", sep="\t", index=False)

    for lib in LIBS:
        a = res[(res.library == lib) & (res.contrast.str.startswith("A"))].set_index("pathway")
        b = res[(res.library == lib) & (res.contrast.str.startswith("B"))].set_index("pathway")
        common = a.index.intersection(b.index)
        r, p = spearmanr(a.loc[common, "nes"], b.loc[common, "nes"])
        up_a = a[(a.fdr < 0.05) & (a.nes > 0)].index
        down_b = b[(b.fdr < 0.05) & (b.nes < 0)].index
        both = up_a.intersection(down_b)
        print(f"\n[{lib}] A: {len(up_a)} up in active (FDR<0.05), {int(((a.fdr < 0.05) & (a.nes < 0)).sum())} down; "
              f"B: {len(down_b)} down after anti-TNF (FDR<0.05), {int(((b.fdr < 0.05) & (b.nes > 0)).sum())} up; "
              f"overlap up-in-active & down-after-treatment: {len(both)}; NES correlation across {len(common)} pathways "
              f"rho = {r:.2f} (p = {p:.1e})")
        if lib == "hallmark":
            show = a.loc[common].join(b.loc[common], lsuffix="_A", rsuffix="_B")[["nes_A", "fdr_A", "nes_B", "fdr_B"]]
            print(show.sort_values("nes_A", ascending=False).round(3).head(15).to_string())
        else:
            print("  top overlap:", ", ".join(sorted(both, key=lambda k: -a.loc[k, "nes"])[:12]))
        pd.DataFrame({"pathway": list(both)}).to_csv(TAB / f"gsea_overlap_{lib}.tsv", sep="\t", index=False)

    # Figure: Hallmark NES, both contrasts
    h = res[res.library == "hallmark"].pivot_table(index="pathway", columns="contrast", values="nes")
    f = res[res.library == "hallmark"].pivot_table(index="pathway", columns="contrast", values="fdr")
    cA, cB = "A_active_vs_remission_macrophages", "B_antiTNF_good_post_vs_pre"
    keep = f[(f[cA] < 0.05) | (f[cB] < 0.05)].index
    h = h.loc[keep].sort_values(cA)
    fig, ax = plt.subplots(figsize=(8, max(4, 0.28 * len(h))))
    y = np.arange(len(h))
    ax.barh(y + 0.2, h[cA], height=0.4, color="C3", label="active vs remission (macrophages, E-MTAB-8322)")
    ax.barh(y - 0.2, h[cB], height=0.4, color="C0", label="anti-TNF good responders, week 12 vs baseline (GSE198520)")
    for i, pw in enumerate(h.index):
        for off, c in ((0.2, cA), (-0.2, cB)):
            if f.loc[pw, c] < 0.05:
                ax.text(h.loc[pw, c] + (0.05 if h.loc[pw, c] >= 0 else -0.15), i + off, "*", va="center", fontsize=9)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y, [p.replace("HALLMARK_", "") for p in h.index], fontsize=7)
    ax.set_xlabel("Normalised enrichment score (* FDR < 0.05)")
    ax.legend(fontsize=7, loc="lower right")
    ax.set_title("Hallmark pathways significant in either contrast", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "gsea_hallmark.png", dpi=150)


if __name__ == "__main__":
    main()
