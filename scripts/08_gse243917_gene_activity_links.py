"""GSE243917: gene activity, RNA-chromatin concordance, peak-gene links, and the active-RA
macrophage signature at the chromatin level.

- Gene activity = ATAC counts in peaks whose midpoint lies in the gene body or 2 kb upstream of the TSS
  (Ensembl 98 = 10x GRCh38-2020-A annotation).
- Pseudobulk units = author subtype x donor with >= 30 nuclei.
- Peak-gene links: Pearson r between peak and gene log-CPM across units, peaks within +/-250 kb of the TSS;
  empirical p from a trans (different-chromosome) null of random peak-gene pairs; BH FDR.
- Signature: the FDR genes up in active RA macrophages (E-MTAB-8322, step 05) are examined in GSE243917
  myeloid nuclei on both layers, and scored per donor on each layer.

Input : data/interim/gse243917_annotated.h5ad, gse243917_atac.h5ad, results/tables/signature_genes.tsv,
        data/raw/annotation/Homo_sapiens.GRCh38.98.gtf.gz, data/processed/gse243917_donor_features.tsv
Output: results/tables/rna_atac_concordance.tsv, peak_gene_links.tsv.gz, signature_chromatin.tsv,
                       signature_links.tsv, multiome_donor_signature.tsv
        data/processed/gse243917_multiome_donor_features.tsv
        results/figures/signature_rna_vs_atac.png
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
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from config import RAW, INTERIM, PROCESSED, FIG, TAB, SEED

UPSTREAM = 2000
WINDOW = 250_000
MIN_UNIT = 30
MIN_MAC = 50
DC_LIKE = ["M-9", "M-10", "M-11", "M-12", "M-13", "M-14"]


def read_genes():
    rows = []
    with gzip.open(RAW / "annotation" / "Homo_sapiens.GRCh38.98.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene":
                continue
            name = re.search(r'gene_name "([^"]+)"', f[8])
            btype = re.search(r'gene_biotype "([^"]+)"', f[8])
            if name is None or btype.group(1) not in ("protein_coding", "lncRNA"):
                continue
            rows.append(("chr" + f[0], int(f[3]), int(f[4]), f[6], name.group(1)))
    g = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "gene"]).drop_duplicates("gene")
    g["tss"] = np.where(g["strand"] == "+", g["start"], g["end"])
    g["reg_start"] = np.where(g["strand"] == "+", g["start"] - UPSTREAM, g["start"])
    g["reg_end"] = np.where(g["strand"] == "+", g["end"], g["end"] + UPSTREAM)
    return g.set_index("gene")


def peak_gene_matrix(peaks, genes):
    """Binary peaks x genes matrix: peak midpoint within gene body + upstream promoter."""
    mid = ((peaks["start"] + peaks["end"]) // 2).values
    rows, cols = [], []
    gi = {g: j for j, g in enumerate(genes.index)}
    for chrom, gsub in genes.groupby("chrom"):
        pidx = np.where(peaks["chrom"].values == chrom)[0]
        if len(pidx) == 0:
            continue
        order = np.argsort(mid[pidx]); pidx = pidx[order]; m = mid[pidx]
        lo = np.searchsorted(m, gsub["reg_start"].values, "left")
        hi = np.searchsorted(m, gsub["reg_end"].values, "right")
        for g, a, b in zip(gsub.index, lo, hi):
            if b > a:
                rows.extend(pidx[a:b]); cols.extend([gi[g]] * (b - a))
    return sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=(len(peaks), len(genes)))


def pseudobulk(X, units, order):
    """Sum rows of X by unit label, return dense units x features in the given order."""
    codes = pd.Categorical(units, categories=order).codes
    keep = codes >= 0
    M = sp.csr_matrix((np.ones(keep.sum()), (codes[keep], np.where(keep)[0])), shape=(len(order), X.shape[0]))
    out = M @ X
    return np.asarray(out.todense()) if sp.issparse(out) else out


def logcpm(m):
    return np.log2(m / m.sum(axis=1, keepdims=True) * 1e6 + 1)


def zcols(m):
    sd = m.std(axis=0, ddof=1)
    sd[sd == 0] = np.nan
    return (m - m.mean(axis=0)) / sd


def main():
    rng = np.random.default_rng(SEED)
    rna = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    atac = sc.read_h5ad(INTERIM / "gse243917_atac.h5ad")
    common = rna.obs_names.intersection(atac.obs_names)
    rna, atac = rna[common], atac[common]
    obs = rna.obs.copy()
    obs["unit"] = obs["author_subtype"].astype(str) + "|" + obs["donor"].astype(str)
    vc = obs["unit"].value_counts()
    units = sorted(vc[vc >= MIN_UNIT].index)
    print(f"Pseudobulk units (subtype x donor, >= {MIN_UNIT} nuclei): {len(units)}")

    genes = read_genes()
    genes = genes[genes.index.isin(rna.var_names)]
    peaks = atac.var[["chrom", "start", "end"]]
    P2G = peak_gene_matrix(peaks, genes)
    C_atac = atac.layers["counts"].tocsr()
    GA = C_atac @ P2G                                     # nuclei x genes gene activity
    print(f"Gene activity: {genes.shape[0]} genes, {(np.asarray(P2G.sum(axis=0)).ravel() > 0).sum()} with >= 1 peak")

    rna_pb = logcpm(pseudobulk(rna[:, genes.index].layers["counts"], obs["unit"], units) + 0.0)
    ga_pb = logcpm(pseudobulk(GA, obs["unit"], units) + 0.0)
    peak_pb = logcpm(pseudobulk(C_atac, obs["unit"], units) + 0.0)

    # 1. RNA vs gene-activity concordance across units, per gene (expressed, variable genes)
    expressed = (rna_pb.mean(axis=0) > 3) & (np.asarray(P2G.sum(axis=0)).ravel() > 0)
    zr, zg = zcols(rna_pb[:, expressed]), zcols(ga_pb[:, expressed])
    r_gene = np.nansum(zr * zg, axis=0) / (len(units) - 1)
    conc = pd.DataFrame({"gene": genes.index[expressed], "r_rna_vs_gene_activity": r_gene}).sort_values("r_rna_vs_gene_activity")
    conc.round(4).to_csv(TAB / "rna_atac_concordance.tsv", sep="\t", index=False)
    print(f"RNA vs gene activity across {len(units)} units: median per-gene r = {np.nanmedian(r_gene):.2f} "
          f"over {expressed.sum()} expressed genes; {np.mean(r_gene > 0.5):.1%} with r > 0.5")

    # 2. Peak-gene links within +/-250 kb of the TSS
    zp = zcols(peak_pb)
    zg_all = zcols(rna_pb)
    gene_ok = (rna_pb.mean(axis=0) > 3) & np.isfinite(zg_all).all(axis=0)
    mid = ((peaks["start"] + peaks["end"]) // 2).values
    pchrom = peaks["chrom"].values
    link_rows = []
    for j in np.where(gene_ok)[0]:
        g = genes.iloc[j]
        idx = np.where((pchrom == g.chrom) & (np.abs(mid - g.tss) <= WINDOW))[0]
        idx = idx[np.isfinite(zp[:, idx]).all(axis=0)]
        if len(idx) == 0:
            continue
        r = zp[:, idx].T @ zg_all[:, j] / (len(units) - 1)
        link_rows.append(pd.DataFrame({"gene": genes.index[j], "peak": peaks.index[idx], "r": r,
                                       "distance_to_tss": mid[idx] - g.tss}))
    links = pd.concat(link_rows, ignore_index=True)
    # trans null: random gene / peak pairs on different chromosomes
    gj = rng.choice(np.where(gene_ok)[0], 200_000)
    pk = rng.choice(np.where(np.isfinite(zp).all(axis=0))[0], 200_000)
    trans = pchrom[pk] != genes["chrom"].values[gj]
    null = np.sort(np.einsum("ij,ij->j", zp[:, pk[trans]], zg_all[:, gj[trans]]) / (len(units) - 1))
    links["p_empirical"] = 1 - np.searchsorted(null, links["r"].values, "left") / len(null)
    links["p_empirical"] = links["p_empirical"].clip(lower=1 / len(null))
    links["fdr"] = multipletests(links["p_empirical"], method="fdr_bh")[1]
    sig = links[(links["r"] > 0.45) & (links["fdr"] < 0.05)]
    links.round(4).to_csv(TAB / "peak_gene_links.tsv.gz", sep="\t", index=False)
    print(f"Peak-gene tests: {len(links):,} cis pairs for {links.gene.nunique():,} genes; trans-null 99th pct r = "
          f"{np.quantile(null, 0.99):.2f}; significant links (r > 0.45, FDR < 0.05): {len(sig):,} "
          f"covering {sig.gene.nunique():,} genes; median |distance| {sig.distance_to_tss.abs().median() / 1000:.0f} kb")

    # 3. Active-RA macrophage signature on both layers (myeloid subtypes, pooled over donors)
    sg = pd.read_csv(TAB / "signature_genes.tsv", sep="\t")
    sig_genes = sg.loc[(sg["signature"] == "primary_FDR") & (sg["direction"] == "up_in_active"), "gene"].tolist()
    sig_genes = [g for g in sig_genes if g in genes.index]
    gpos = [genes.index.get_loc(g) for g in sig_genes]
    my = obs["author_ct"] == "myeloid"
    subsets = obs.loc[my, "author_subtype"].astype(str)
    sub_order = [s for s in subsets.value_counts().index if subsets.value_counts()[s] >= 30]
    my_idx = np.where(my.values)[0]
    rna_sub = logcpm(pseudobulk(rna[my_idx][:, genes.index].layers["counts"], subsets.values, sub_order) + 0.0)[:, gpos]
    ga_sub = logcpm(pseudobulk(GA[my_idx], subsets.values, sub_order) + 0.0)[:, gpos]
    chrom_tab = []
    for k, g in enumerate(sig_genes):
        r, p = spearmanr(rna_sub[:, k], ga_sub[:, k])
        nl = sig[sig["gene"] == g]
        chrom_tab.append(dict(gene=g, n_peaks_in_gene_region=int(P2G[:, genes.index.get_loc(g)].sum()),
                              rho_rna_vs_activity_across_myeloid_subsets=r, p=p,
                              top_rna_subset=sub_order[int(np.argmax(rna_sub[:, k]))],
                              top_activity_subset=sub_order[int(np.argmax(ga_sub[:, k]))],
                              n_linked_peaks=len(nl), max_link_r=nl["r"].max() if len(nl) else np.nan))
    chrom_tab = pd.DataFrame(chrom_tab)
    chrom_tab.round(3).to_csv(TAB / "signature_chromatin.tsv", sep="\t", index=False)
    sig[sig["gene"].isin(sig_genes)].sort_values(["gene", "r"], ascending=[True, False]).round(4) \
        .to_csv(TAB / "signature_links.tsv", sep="\t", index=False)
    print(chrom_tab.round(2).to_string(index=False))

    # 4. Donor-level signature on each layer (macrophages, donors with >= 50)
    sub_code = obs["author_subtype"].astype(str).str.split(":").str[0]
    mac = (obs["author_ct"] == "myeloid") & ~sub_code.isin(DC_LIKE)
    n_mac = obs.loc[mac, "donor"].value_counts()
    donors = sorted(n_mac[n_mac >= MIN_MAC].index)
    mac_idx = np.where(mac.values)[0]
    dons = obs["donor"].astype(str).values[mac_idx]
    rna_d = zcols(logcpm(pseudobulk(rna[mac_idx][:, genes.index].layers["counts"], dons, donors) + 0.0))
    ga_d = zcols(logcpm(pseudobulk(GA[mac_idx], dons, donors) + 0.0))
    sig_peaks = [p for p in sig.loc[sig["gene"].isin(sig_genes), "peak"].unique()]
    pk_pos = [peaks.index.get_loc(p) for p in sig_peaks]
    pk_d = zcols(logcpm(pseudobulk(C_atac[mac_idx], dons, donors) + 0.0)[:, pk_pos]) if pk_pos else None
    feats = pd.read_csv(PROCESSED / "gse243917_donor_features.tsv", sep="\t", index_col=0)
    d = pd.DataFrame(index=donors)
    d["rna_signature_score"] = np.nanmean(rna_d[:, gpos], axis=1)
    d["atac_gene_activity_signature_score"] = np.nanmean(ga_d[:, gpos], axis=1)
    if pk_d is not None:
        d["atac_linked_peak_signature_score"] = np.nanmean(pk_d, axis=1)
    d = feats[["disease", "lymphocytic_inflammation_grade", "lining_hyperplasia_grade", "prep_batch"]].join(d, how="right")
    d.round(4).to_csv(TAB / "multiome_donor_signature.tsv", sep="\t")
    ra = d[d["disease"] == "RA"]
    stats = []
    for c in [c for c in d.columns if c.endswith("_score")]:
        r_i, p_i = spearmanr(ra[c], ra["lymphocytic_inflammation_grade"])
        r_r, p_r = spearmanr(ra[c], ra["rna_signature_score"])
        p_b = mannwhitneyu(ra.loc[ra.prep_batch == "lowMT", c], ra.loc[ra.prep_batch == "highMT", c]).pvalue
        stats.append(dict(score=c, n_RA=len(ra), rho_vs_inflammation=r_i, p_inflammation=p_i,
                          rho_vs_rna_score=r_r, p_vs_rna=p_r, p_batch=p_b))
    stats = pd.DataFrame(stats)

    # Control: same nuclei feed both layers, so donor-wide quality could couple any gene set's RNA and ATAC
    # scores. Null = random gene sets of the same size, matched on RNA expression decile.
    ra_rows = np.array([i for i, dn in enumerate(donors) if d.loc[dn, "disease"] == "RA"])
    expr_mean = logcpm(pseudobulk(rna[mac_idx][:, genes.index].layers["counts"], dons, donors) + 0.0).mean(axis=0)
    has_peak = np.asarray(P2G.sum(axis=0)).ravel() > 0
    pool_ok = np.isfinite(rna_d).all(axis=0) & np.isfinite(ga_d).all(axis=0) & has_peak & (expr_mean > 3)
    dec = pd.qcut(expr_mean[pool_ok], 10, labels=False)
    pool_idx = np.where(pool_ok)[0]
    sig_dec = [dec[np.searchsorted(pool_idx, j)] if j in set(pool_idx) else None for j in gpos]
    null_rho = []
    for _ in range(1000):
        pick = [rng.choice(pool_idx[dec == q]) for q in sig_dec if q is not None]
        a = np.nanmean(rna_d[np.ix_(ra_rows, pick)], axis=1)
        b = np.nanmean(ga_d[np.ix_(ra_rows, pick)], axis=1)
        null_rho.append(spearmanr(a, b)[0])
    null_rho = np.array(null_rho)
    obs_rho = stats.set_index("score").loc["atac_gene_activity_signature_score", "rho_vs_rna_score"]
    p_null = (1 + np.sum(null_rho >= obs_rho)) / (1 + len(null_rho))
    stats["null_median_rho_vs_rna"] = np.where(stats["score"] == "atac_gene_activity_signature_score",
                                               np.median(null_rho), np.nan)
    stats["p_vs_random_genesets"] = np.where(stats["score"] == "atac_gene_activity_signature_score", p_null, np.nan)
    stats.round(4).to_csv(TAB / "multiome_signature_stats.tsv", sep="\t", index=False)
    print(stats.round(3).to_string(index=False))
    print(f"Random expression-matched gene sets (n=1000): median rho(RNA, ATAC) = {np.median(null_rho):.2f}, "
          f"95th pct = {np.quantile(null_rho, 0.95):.2f}; observed {obs_rho:.2f}, empirical p = {p_null:.3f}")

    # Which myeloid subset carries the signature? mean z across signature genes, per layer
    zr_sub = (rna_sub - rna_sub.mean(axis=0)) / (rna_sub.std(axis=0) + 1e-9)
    zg_sub = (ga_sub - ga_sub.mean(axis=0)) / (ga_sub.std(axis=0) + 1e-9)
    rank = pd.DataFrame({"subset": sub_order, "rna_signature_z": zr_sub.mean(axis=1),
                         "atac_signature_z": zg_sub.mean(axis=1)}).sort_values("rna_signature_z", ascending=False)
    rank.round(3).to_csv(TAB / "signature_by_myeloid_subset.tsv", sep="\t", index=False)
    print(rank.round(2).to_string(index=False))

    # Integrated donor feature table (RNA composition + both signature layers)
    feats.join(d.drop(columns=["disease", "lymphocytic_inflammation_grade", "lining_hyperplasia_grade", "prep_batch"])) \
        .round(4).to_csv(PROCESSED / "gse243917_multiome_donor_features.tsv", sep="\t")

    # Figure: RNA vs gene activity across myeloid subsets for the signature genes; donor-level scores
    fig, axes = plt.subplots(1, 3, figsize=(20, 5.5), gridspec_kw={"width_ratios": [1, 1, 0.9]})
    for ax, m, title in ((axes[0], rna_sub, "RNA (log-CPM, z by gene)"), (axes[1], ga_sub, "Gene activity (ATAC, z by gene)")):
        z = (m - m.mean(axis=0)) / (m.std(axis=0) + 1e-9)
        im = ax.imshow(z.T, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
        ax.set_xticks(range(len(sub_order)), [s.split(":")[0] + ":" + s.split(":")[1][:14] for s in sub_order],
                      rotation=90, fontsize=7)
        ax.set_yticks(range(len(sig_genes)), sig_genes, fontsize=8)
        ax.set_title(f"Active-RA macrophage genes in GSE243917 myeloid subsets\n{title}", fontsize=9)
        plt.colorbar(im, ax=ax, shrink=0.6)
    ax = axes[2]
    col = "atac_gene_activity_signature_score"
    for b, mk in (("highMT", "o"), ("lowMT", "s")):
        s = d[d["prep_batch"] == b]
        ax.scatter(s["rna_signature_score"], s[col], marker=mk, s=55, label=b,
                   c=["grey" if x == "OA" else "C3" for x in s["disease"]])
        for name, row in s.iterrows():
            ax.annotate(f"{name.replace('BRI-', '')} g{row.lymphocytic_inflammation_grade}", (row.rna_signature_score, row[col]),
                        fontsize=6, xytext=(3, 3), textcoords="offset points")
    st = stats.set_index("score").loc[col]
    ax.set_xlabel("RNA signature score (macrophage pseudobulk)")
    ax.set_ylabel("ATAC gene-activity signature score")
    ax.set_title(f"Donor level, RA n={int(st.n_RA)}: rho(RNA, ATAC) = {st.rho_vs_rna_score:.2f}, p = {st.p_vs_rna:.3f}",
                 fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(FIG / "signature_rna_vs_atac.png", dpi=150)


if __name__ == "__main__":
    main()
