"""TF motif accessibility (chromVAR deviations, pychromvar + JASPAR CORE vertebrates) on pseudobulk units.

Single-nucleus chromVAR would need a dense 31k x 132k matrix (> 16 GB), so deviations are computed on units that
share one background model:
  - subtype x donor units with >= 10 nuclei (author subtypes)
  - one macrophage unit per donor (all non-DC myeloid nuclei)
Questions:
  A. sanity: lineage-defining motifs per broad cell type (expected SPI1/CEBP myeloid, EBF1/PAX5 B, RUNX/ETS T ...);
  B. M-6 STAT1+ CXCL10+ vs the other macrophage subtypes of the same donor (paired, donors with an M-6 unit);
  C. donor level: macrophage motif accessibility vs the RNA active-RA signature score of the same macrophages
     (step 08), and whether STAT/IRF motifs rank high among all motifs.

Input : data/interim/gse243917_atac.h5ad, data/interim/hg38.fa (UCSC), results/tables/multiome_donor_signature.tsv
Output: results/tables/chromvar_unit_deviations.tsv.gz, chromvar_celltype_top.tsv, chromvar_M6_vs_macrophages.tsv,
        chromvar_donor_signature_link.tsv; results/figures/chromvar.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import pysam
import pychromvar as pc
from pyjaspar import jaspardb
from scipy.stats import spearmanr, wilcoxon

from config import INTERIM, FIG, TAB, SEED

MIN_UNIT = 10
M6 = "M-6: STAT1+ CXCL10+"
DC_LIKE = ["M-9", "M-10", "M-11", "M-12", "M-13", "M-14"]
FOCUS = r"^(STAT|IRF|NFKB|REL|SPI1|CEBP|MAFB|NR4A|FOS|JUN|BATF|ATF3|PRDM1|MITF|TFEC)"


def build_units(atac):
    obs = atac.obs
    sub = obs["author_subtype"].astype(str)
    code = sub.str.split(":").str[0]
    lab = (sub + "|" + obs["donor"].astype(str)).values
    units = pd.Series(lab).value_counts()
    units = list(units[units >= MIN_UNIT].index)
    mac = (obs["author_ct"] == "myeloid").values & ~code.isin(DC_LIKE).values
    labels = np.where(np.isin(lab, units), lab, None).astype(object)
    groups = units + [f"MACROPHAGES|{d}" for d in sorted(obs["donor"].unique())]
    C = atac.layers["counts"].tocsr()
    rows = []
    for g in groups:
        if g.startswith("MACROPHAGES|"):
            idx = np.where(mac & (obs["donor"].values == g.split("|")[1]))[0]
        else:
            idx = np.where(labels == g)[0]
        rows.append((g, idx))
    rows = [(g, i) for g, i in rows if len(i) >= MIN_UNIT]
    M = sp.csr_matrix((np.ones(sum(len(i) for _, i in rows)),
                       (np.concatenate([[k] * len(i) for k, (_, i) in enumerate(rows)]), np.concatenate([i for _, i in rows]))),
                      shape=(len(rows), C.shape[0]))
    X = (M @ C).toarray().astype(np.float32)
    meta = pd.DataFrame({"unit": [g for g, _ in rows], "n_nuclei": [len(i) for _, i in rows]}).set_index("unit")
    meta["state"] = meta.index.str.split("|").str[0]
    meta["donor"] = meta.index.str.split("|").str[1]
    ct_of = obs.groupby(sub.values)["author_ct"].first()
    meta["cell_type"] = [("macrophage_pool" if s == "MACROPHAGES" else ct_of.get(s, "?")) for s in meta["state"]]
    ad = sc.AnnData(X=X, obs=meta, var=pd.DataFrame(index=atac.var_names))
    ad = ad[:, X.sum(axis=0) > 0].copy()
    return ad


def main():
    atac = sc.read_h5ad(INTERIM / "gse243917_atac.h5ad")
    ad = build_units(atac)
    del atac
    print(f"Units: {ad.n_obs} ({(ad.obs.state == 'MACROPHAGES').sum()} donor macrophage pools), peaks: {ad.n_vars:,}")

    fa = INTERIM / "hg38.fa"
    if not (INTERIM / "hg38.fa.fai").exists():
        pysam.faidx(str(fa))
    pc.add_peak_seq(ad, genome_file=str(fa), delimiter=":|-")
    pc.add_gc_bias(ad)
    np.random.seed(SEED)
    pc.get_bg_peaks(ad, niterations=50, n_jobs=8)
    motifs = jaspardb(release="JASPAR2026").fetch_motifs(collection="CORE", tax_group=["vertebrates"])
    pc.match_motif(ad, motifs=motifs, p_value=5e-5)
    dev = pc.compute_deviations(ad, n_jobs=8)
    dev.obs = ad.obs.copy()
    names = pd.Series(dev.var_names).str.split(".").str[-1].values     # "MA0080.7.SPI1" -> "SPI1"
    z = pd.DataFrame(dev.X, index=dev.obs_names, columns=[f"{n}|{i}" for n, i in zip(names, dev.var_names)])
    z.round(4).to_csv(TAB / "chromvar_unit_deviations.tsv.gz", sep="\t")
    print(f"Motifs: {z.shape[1]}")

    # A. sanity check per broad cell type
    meta = dev.obs
    rows = []
    for ct in sorted(meta["cell_type"].unique()):
        if ct in ("macrophage_pool", "?"):
            continue
        m = z[meta["cell_type"] == ct].mean() - z[(meta["cell_type"] != ct) & (meta["cell_type"] != "macrophage_pool")].mean()
        for k, v in m.sort_values(ascending=False).head(8).items():
            rows.append(dict(cell_type=ct, motif=k.split("|")[0], jaspar_id=k.split("|")[1], delta_z=v))
    top = pd.DataFrame(rows)
    top.round(3).to_csv(TAB / "chromvar_celltype_top.tsv", sep="\t", index=False)
    print("\nA. Top motifs per cell type:")
    print(top.groupby("cell_type")["motif"].apply(lambda s: ", ".join(s.head(6))).to_string())

    # B. M-6 vs other macrophage subtypes, paired within donor
    is_mac_sub = meta["cell_type"].eq("myeloid") & ~meta["state"].str.split(":").str[0].isin(DC_LIKE)
    diffs = []
    for d in sorted(meta["donor"].unique()):
        m6 = meta.index[(meta["state"] == M6) & (meta["donor"] == d)]
        oth = meta.index[is_mac_sub & (meta["state"] != M6) & (meta["donor"] == d)]
        if len(m6) and len(oth):
            diffs.append((z.loc[m6].mean() - z.loc[oth].mean()).rename(d))
    if diffs:
        D = pd.DataFrame(diffs)
        res_b = pd.DataFrame({"mean_diff_z": D.mean(), "n_donors": len(D), "n_positive": (D > 0).sum()})
        if len(D) >= 3:
            res_b["p_wilcoxon"] = [wilcoxon(D[c]).pvalue if (D[c] != 0).any() else 1.0 for c in D.columns]
        res_b["motif"] = [c.split("|")[0] for c in res_b.index]
        res_b = res_b.sort_values("mean_diff_z", ascending=False)
        res_b.round(4).to_csv(TAB / "chromvar_M6_vs_macrophages.tsv", sep="\t")
        print(f"\nB. M-6 vs other macrophages (paired within {len(D)} donors with an M-6 unit of >= {MIN_UNIT} nuclei):")
        print(res_b.head(15).round(3).to_string())
        foc = res_b[res_b["motif"].str.match(FOCUS)]
        print("  focus motifs (STAT/IRF/NF-kB/AP-1/myeloid), best 10:")
        print(foc.head(10).round(3).to_string())
    else:
        res_b = None
        print("\nB. No donor has an M-6 unit with enough nuclei")

    # C. donor-level macrophage motif accessibility vs RNA signature score
    pools = meta.index[meta["state"] == "MACROPHAGES"]
    zp = z.loc[pools].copy()
    zp.index = meta.loc[pools, "donor"].values
    sig = pd.read_csv(TAB / "multiome_donor_signature.tsv", sep="\t", index_col=0)
    common = [d for d in zp.index if d in sig.index and sig.loc[d, "disease"] == "RA"]
    rs = pd.Series({c: spearmanr(zp.loc[common, c], sig.loc[common, "rna_signature_score"])[0] for c in zp.columns})
    link = pd.DataFrame({"rho_vs_rna_signature": rs, "motif": [c.split("|")[0] for c in rs.index]}).sort_values(
        "rho_vs_rna_signature", ascending=False)
    link["rank"] = np.arange(1, len(link) + 1)
    link.round(4).to_csv(TAB / "chromvar_donor_signature_link.tsv", sep="\t")
    stat_irf = link["motif"].str.match(r"^(STAT|IRF)")
    from scipy.stats import mannwhitneyu
    p_set = mannwhitneyu(link.loc[stat_irf, "rho_vs_rna_signature"], link.loc[~stat_irf, "rho_vs_rna_signature"],
                         alternative="greater").pvalue
    print(f"\nC. Donor level (RA, n = {len(common)}): motifs most correlated with macrophage RNA signature score:")
    print(link.head(12).round(3).to_string())
    print(f"  STAT/IRF motifs ({stat_irf.sum()}): median rho {link.loc[stat_irf, 'rho_vs_rna_signature'].median():.2f} "
          f"vs others {link.loc[~stat_irf, 'rho_vs_rna_signature'].median():.2f}; one-sided Mann-Whitney p = {p_set:.3g}")
    pd.DataFrame([dict(test="STAT/IRF motifs vs all other motifs: correlation with RNA signature (donor level)",
                       n_donors=len(common), n_stat_irf=int(stat_irf.sum()),
                       median_rho_stat_irf=link.loc[stat_irf, "rho_vs_rna_signature"].median(),
                       median_rho_other=link.loc[~stat_irf, "rho_vs_rna_signature"].median(), p_one_sided=p_set)]) \
        .round(4).to_csv(TAB / "chromvar_tests.tsv", sep="\t", index=False)

    # Figure
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    ax = axes[0]
    piv = top.pivot_table(index="motif", columns="cell_type", values="delta_z", aggfunc="max").fillna(0)
    piv = piv.loc[piv.max(axis=1).sort_values(ascending=False).index[:30]]
    im = ax.imshow(piv.values, aspect="auto", cmap="Reds")
    ax.set_xticks(range(piv.shape[1]), piv.columns, rotation=60, fontsize=8)
    ax.set_yticks(range(piv.shape[0]), piv.index, fontsize=7)
    ax.set_title("A. Lineage motifs (delta chromVAR z vs other cell types)", fontsize=10)
    plt.colorbar(im, ax=ax, shrink=0.6)
    ax = axes[1]
    if res_b is not None:
        rb = res_b.head(20).iloc[::-1]
        ax.barh(range(len(rb)), rb["mean_diff_z"], color=["C3" if m else "C7" for m in rb["motif"].str.match(r"^(STAT|IRF)")])
        ax.set_yticks(range(len(rb)), [f"{m} ({p}/{n})" for m, p, n in zip(rb["motif"], rb["n_positive"], rb["n_donors"])], fontsize=7)
        ax.set_xlabel("M-6 minus other macrophages, mean delta z (donors positive / total)")
        ax.set_title("B. Motifs more accessible in STAT1+ CXCL10+ macrophages\n(red = STAT/IRF)", fontsize=10)
    ax = axes[2]
    ax.hist(link["rho_vs_rna_signature"], bins=40, color="C7", alpha=0.7, label="all motifs")
    ax.hist(link.loc[stat_irf, "rho_vs_rna_signature"], bins=20, color="C3", alpha=0.8, label="STAT/IRF motifs")
    ax.set_xlabel("Donor-level Spearman rho with macrophage RNA signature score")
    ax.set_title(f"C. Macrophage motif accessibility vs RNA signature (RA n = {len(common)})\n"
                 f"STAT/IRF higher than other motifs: p = {p_set:.2g}", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "chromvar.png", dpi=150)


if __name__ == "__main__":
    main()
