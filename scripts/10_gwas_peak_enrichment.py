"""Genetics layer: are RA risk loci enriched near cell-state-specific accessible chromatin?

Design (locus-window enrichment, exploratory alternative to S-LDSC, which needs full summary statistics):
- GWAS: GWAS Catalog associations whose mapped traits are all RA terms (RA, ACPA+/RF+ RA), p <= 5e-8,
  GRCh38; MHC (chr6:28.51-33.48 Mb) excluded. Lead SNPs +/- W are merged into loci. Negative control: LDL.
- Peak sets: top N specific peaks per state. Broad types: each vs all other nuclei. Subtypes: each vs the other
  subtypes of the same lineage (so the test is about the state, not the lineage).
- Statistic: number of loci with >= 1 state-specific peak. Null: 1,000 random peak sets of the same size matched
  on accessibility decile x promoter/distal. Enrichment = observed / null mean; empirical p; BH within level x trait.
- Windows: 25 kb primary; 10 kb and 100 kb sensitivity.

Input : data/interim/gse243917_atac.h5ad, data/raw/gwas/*_associations.tsv, Ensembl 98 GTF,
        results/tables/peak_gene_links.tsv.gz
Output: results/tables/gwas_loci_{RA,LDL}.tsv, gwas_enrichment.tsv, gwas_candidate_variants.tsv
        results/figures/gwas_enrichment.png
"""
import gzip
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

from config import RAW, INTERIM, FIG, TAB, SEED

N_TOP = 4000
N_PERM = 1000
WINDOWS = (10_000, 25_000, 100_000)
PRIMARY_W = 25_000
MIN_STATE = 100
MHC = ("chr6", 28_510_120, 33_480_577)
RA_TERMS = re.compile(r"rheumatoid|citrullinated|rheumatoid factor", re.I)


def lead_snps(key):
    d = pd.read_csv(RAW / "gwas" / f"{key}_associations.tsv", sep="\t")
    d = d[(d["p_value"] <= 5e-8) & d["location"].notna() & d["location"].astype(str).str.contains(":")]
    if key == "RA":   # keep associations whose traits are all RA-related (drop cross-disease pleiotropy scans)
        d = d[d["efo_traits"].fillna("").apply(lambda s: all(RA_TERMS.search(t) for t in s.split(";")))]
    d["chrom"] = "chr" + d["location"].str.split(":").str[0]
    d["pos"] = pd.to_numeric(d["location"].str.split(":").str[1], errors="coerce")
    d = d.dropna(subset=["pos"])
    d["pos"] = d["pos"].astype(int)
    d = d[~((d["chrom"] == MHC[0]) & d["pos"].between(MHC[1], MHC[2]))]
    return d.sort_values("p_value").drop_duplicates("rs_id")


def merge_loci(snps, w):
    loci = []
    for chrom, g in snps.sort_values("pos").groupby("chrom"):
        s0, e0, members = None, None, []
        for _, r in g.iterrows():
            s, e = r["pos"] - w, r["pos"] + w
            if s0 is not None and s <= e0:
                e0 = max(e0, e); members.append(r["rs_id"])
            else:
                if s0 is not None:
                    loci.append((chrom, s0, e0, members))
                s0, e0, members = s, e, [r["rs_id"]]
        loci.append((chrom, s0, e0, members))
    return pd.DataFrame(loci, columns=["chrom", "start", "end", "snps"])


def peak_locus_ids(peaks, loci):
    mid = ((peaks["start"] + peaks["end"]) // 2).values
    lid = np.full(len(peaks), -1)
    for k, r in loci.reset_index(drop=True).iterrows():
        hit = (peaks["chrom"].values == r["chrom"]) & (mid >= r["start"]) & (mid <= r["end"])
        lid[hit] = k
    return lid


def tss_table():
    rows = []
    with gzip.open(RAW / "annotation" / "Homo_sapiens.GRCh38.98.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene" or 'gene_biotype "protein_coding"' not in f[8]:
                continue
            name = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
            rows.append(("chr" + f[0], int(f[3]) if f[6] == "+" else int(f[4]), name))
    return pd.DataFrame(rows, columns=["chrom", "tss", "gene"])


def pseudobulk_cpm(C, labels, order):
    codes = pd.Categorical(labels, categories=order).codes
    keep = codes >= 0
    M = sp.csr_matrix((np.ones(keep.sum()), (codes[keep], np.where(keep)[0])), shape=(len(order), C.shape[0]))
    pb = np.asarray((M @ C).todense())
    return pb / pb.sum(axis=1, keepdims=True) * 1e6


def specific_sets(atac, ok_peak):
    C = atac.layers["counts"].tocsr()[:, ok_peak]
    ct = atac.obs["author_ct"].astype(str).values
    st = atac.obs["author_subtype"].astype(str).values
    sets = {}
    cts = [c for c in pd.unique(ct) if (ct == c).sum() >= MIN_STATE]
    cpm = pseudobulk_cpm(C, ct, cts)
    for i, c in enumerate(cts):
        others = np.delete(cpm, i, axis=0).mean(axis=0)
        score = np.log2((cpm[i] + 1) / (others + 1))
        score[cpm[i] < np.quantile(cpm[i], 0.5)] = -np.inf          # must be reasonably open in the state
        sets[("cell type", c)] = np.argsort(score)[::-1][:N_TOP]
    for parent in cts:
        subs = [s for s in pd.unique(st[ct == parent]) if (st == s).sum() >= MIN_STATE]
        if len(subs) < 2:
            continue
        m = ct == parent
        cpm_s = pseudobulk_cpm(C[m], st[m], subs)
        for i, s in enumerate(subs):
            others = np.delete(cpm_s, i, axis=0).mean(axis=0)
            score = np.log2((cpm_s[i] + 1) / (others + 1))
            score[cpm_s[i] < np.quantile(cpm_s[i], 0.5)] = -np.inf
            sets[("subtype", s)] = np.argsort(score)[::-1][:N_TOP]
    return sets


def main():
    rng = np.random.default_rng(SEED)
    atac = sc.read_h5ad(INTERIM / "gse243917_atac.h5ad")
    peaks = atac.var[["chrom", "start", "end", "frac_cells"]].copy()
    ok_peak = ~((peaks["chrom"] == MHC[0]) & (peaks["end"] >= MHC[1]) & (peaks["start"] <= MHC[2])).values
    peaks = peaks[ok_peak]
    mid = ((peaks["start"] + peaks["end"]) // 2).values

    tss = tss_table()
    near_tss = np.zeros(len(peaks), bool)
    nearest_gene = np.empty(len(peaks), dtype=object)
    for chrom, g in tss.groupby("chrom"):
        idx = np.where(peaks["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        t = np.sort(g["tss"].values); names = g.sort_values("tss")["gene"].values
        j = np.clip(np.searchsorted(t, mid[idx]), 1, len(t) - 1)
        left_closer = np.abs(mid[idx] - t[j - 1]) <= np.abs(mid[idx] - t[j])
        jj = np.where(left_closer, j - 1, j)
        near_tss[idx] = np.abs(mid[idx] - t[jj]) <= 2000
        nearest_gene[idx] = names[jj]
    strata = pd.qcut(peaks["frac_cells"].rank(method="first"), 10, labels=False).values * 2 + near_tss.astype(int)
    pools = {s: np.where(strata == s)[0] for s in np.unique(strata)}

    sets = specific_sets(atac, ok_peak)
    print(f"Peak sets: {sum(k[0] == 'cell type' for k in sets)} cell types, {sum(k[0] == 'subtype' for k in sets)} subtypes; "
          f"{N_TOP} peaks each; {len(peaks):,} non-MHC peaks in universe")

    snps = {k: lead_snps(k) for k in ("RA", "LDL")}
    rows = []
    for trait, sn in snps.items():
        for w in WINDOWS:
            loci = merge_loci(sn, w)
            if w == PRIMARY_W:
                loci.assign(n_snps=loci["snps"].str.len(), snps=loci["snps"].str.join(",")) \
                    .to_csv(TAB / f"gwas_loci_{trait}.tsv", sep="\t", index=False)
            lid = peak_locus_ids(peaks, loci)
            base = len(np.unique(lid[lid >= 0]))
            print(f"{trait}: {len(sn)} lead SNPs -> {len(loci)} loci at +/-{w // 1000} kb; loci with any peak: {base}")
            for (level, state), idx in sets.items():
                obs_hit = len(np.unique(lid[idx][lid[idx] >= 0]))
                need = pd.Series(strata[idx]).value_counts()
                null = np.empty(N_PERM)
                for b in range(N_PERM):
                    pick = np.concatenate([rng.choice(pools[s], n, replace=False) for s, n in need.items()])
                    h = lid[pick]
                    null[b] = len(np.unique(h[h >= 0]))
                rows.append(dict(trait=trait, window_kb=w // 1000, level=level, state=state, n_loci=len(loci),
                                 loci_hit=obs_hit, null_mean=null.mean(), enrichment=obs_hit / max(null.mean(), 1e-9),
                                 p_empirical=(1 + np.sum(null >= obs_hit)) / (1 + N_PERM)))
    res = pd.DataFrame(rows)
    res["fdr"] = np.nan
    for _, g in res.groupby(["trait", "window_kb", "level"]):
        res.loc[g.index, "fdr"] = multipletests(g["p_empirical"], method="fdr_bh")[1]
    res = res.sort_values(["trait", "window_kb", "level", "p_empirical"])
    res.round(4).to_csv(TAB / "gwas_enrichment.tsv", sep="\t", index=False)
    prim = res[res["window_kb"] == PRIMARY_W // 1000]
    for trait in ("RA", "LDL"):
        print(f"\n== {trait}, +/-{PRIMARY_W // 1000} kb")
        print(prim[prim["trait"] == trait].head(12)[["level", "state", "loci_hit", "null_mean", "enrichment",
                                                     "p_empirical", "fdr"]].round(3).to_string(index=False))
    rob = res[(res["trait"] == "RA")].pivot_table(index=["level", "state"], columns="window_kb", values="fdr")
    print("\nRA FDR by window (states with FDR < 0.1 at any window):")
    print(rob[(rob < 0.1).any(axis=1)].round(3).to_string())

    # Candidate RA variants sitting inside (+/-250 bp of) a state-specific peak, with nearest / linked genes
    links = pd.read_csv(TAB / "peak_gene_links.tsv.gz", sep="\t")
    links = links[(links["r"] > 0.45) & (links["fdr"] < 0.05)]
    ra = snps["RA"]
    cand = []
    peak_names = peaks.index.values
    for (level, state), idx in sets.items():
        sub = peaks.iloc[idx]
        smid = mid[idx]
        for _, s in ra.iterrows():
            hit = np.where((sub["chrom"].values == s["chrom"]) & (np.abs(smid - s["pos"]) <= 250))[0]
            for h in hit:
                pk = peak_names[idx[h]]
                lg = links.loc[links["peak"] == pk, "gene"].tolist()
                cand.append(dict(rs_id=s["rs_id"], location=s["location"], p_value=s["p_value"],
                                 catalog_mapped_genes=s["mapped_genes"], level=level, state=state, peak=pk,
                                 nearest_protein_coding_gene=nearest_gene[idx[h]],
                                 linked_genes=";".join(lg), study=s["first_author"]))
    cand = pd.DataFrame(cand).sort_values(["level", "state", "p_value"])
    cand.to_csv(TAB / "gwas_candidate_variants.tsv", sep="\t", index=False)
    print(f"\nRA lead SNPs within 250 bp of a state-specific peak: {cand['rs_id'].nunique()} SNPs, {len(cand)} SNP-state pairs")
    print(cand.groupby(["level", "state"]).size().sort_values(ascending=False).head(12).to_string())

    # Figure: primary-window enrichment, RA vs LDL, cell types and subtypes
    fig, axes = plt.subplots(1, 2, figsize=(16, 9), gridspec_kw={"width_ratios": [1, 2.4]})
    for ax, level in zip(axes, ("cell type", "subtype")):
        p = prim[prim["level"] == level].pivot_table(index="state", columns="trait", values="enrichment")
        f = prim[prim["level"] == level].pivot_table(index="state", columns="trait", values="fdr")
        p = p.sort_values("RA")
        y = np.arange(len(p))
        ax.barh(y + 0.2, p["RA"], height=0.4, color="C3", label="RA")
        ax.barh(y - 0.2, p["LDL"], height=0.4, color="C7", label="LDL (control)")
        for i, st in enumerate(p.index):
            if f.loc[st, "RA"] < 0.05:
                ax.text(p.loc[st, "RA"] + 0.03, i + 0.2, "*", va="center", fontsize=11)
        ax.axvline(1, color="k", lw=0.8, ls="--")
        ax.set_yticks(y, p.index, fontsize=7 if level == "subtype" else 9)
        ax.set_xlabel(f"Loci hit / matched-null mean (+/-{PRIMARY_W // 1000} kb)")
        ax.set_title(f"{level}: specific peaks vs GWAS loci (* RA FDR < 0.05)", fontsize=10)
        ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "gwas_enrichment.png", dpi=150)


if __name__ == "__main__":
    main()
