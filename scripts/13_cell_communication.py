"""Cell-cell communication: who signals to the IFN-activated (STAT1+ CXCL10+) macrophage state?

1. LIANA rank_aggregate (consensus resource; CellPhoneDB, CellChat, NATMI, Connectome, log2FC, SingleCellSignalR,
   geometric mean) on GSE243917 author subtypes with >= 30 nuclei.
2. Focus: interactions received by M-6 STAT1+ CXCL10+ macrophages, the IFNG -> IFNGR axis across senders, and
   M-6 -> lymphocyte chemokine output (CXCL9/10/11 -> CXCR3).
3. Genetics link: across T/NK subtypes, is IFNG expression higher in the states enriched for RA risk loci (step 10)?
4. Donor-level link (testable, n = RA donors): IFNG in lymphocyte pseudobulk vs the IFN-gamma programme in the same
   donor's macrophages (Hallmark IFN-gamma response genes, excluding IFNG itself).
Ligand-receptor co-expression is evidence of possible, not proven, communication.

Input : data/interim/gse243917_annotated.h5ad, results/tables/gwas_enrichment.tsv, data/raw/genesets/h.all*.gmt
Output: results/tables/liana_results.tsv.gz, liana_to_M6.tsv, ifng_senders.tsv, ifng_donor_link.tsv
        results/figures/cell_communication.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import liana as li
from scipy.stats import spearmanr

from config import RAW, INTERIM, FIG, TAB, SEED

M6 = "M-6: STAT1+ CXCL10+"
MIN_GROUP = 30
DC_LIKE = ["M-9", "M-10", "M-11", "M-12", "M-13", "M-14"]


def main():
    adata = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    sub = adata.obs["author_subtype"].astype(str)
    keep_groups = sub.value_counts()[lambda s: s >= MIN_GROUP].index
    adata = adata[sub.isin(keep_groups) & ~sub.str.contains("low quality|MT-high")].copy()
    adata.obs["state"] = adata.obs["author_subtype"].astype(str).astype("category")
    print(f"LIANA input: {adata.n_obs} nuclei, {adata.obs['state'].nunique()} states")

    cache = TAB / "liana_results.tsv.gz"
    if cache.exists():                                    # permutations take minutes; reuse a finished run
        res = pd.read_csv(cache, sep="\t")
    else:
        li.mt.rank_aggregate(adata, groupby="state", resource_name="consensus", expr_prop=0.1, use_raw=True,
                             n_perms=1000, seed=SEED, verbose=False)
        res = adata.uns["liana_res"].copy()
        res.to_csv(cache, sep="\t", index=False)
    print(f"LIANA interactions scored: {len(res):,}")

    # Received by M-6
    to_m6 = res[res["target"] == M6].sort_values("magnitude_rank")
    to_m6.head(200).to_csv(TAB / "liana_to_M6.tsv", sep="\t", index=False)
    print("\nTop interactions received by M-6 (magnitude rank):")
    print(to_m6.head(15)[["source", "ligand_complex", "receptor_complex", "magnitude_rank", "specificity_rank"]]
          .round(4).to_string(index=False))

    # Specificity view: among interactions with reasonable magnitude, which are most specific to M-6 as receiver?
    # many interactions tie at specificity rank 0, so break ties by magnitude
    spec = to_m6[to_m6["magnitude_rank"] < 0.25].sort_values(["specificity_rank", "magnitude_rank"])
    spec.head(100).to_csv(TAB / "liana_to_M6_specific.tsv", sep="\t", index=False)
    print("\nMost M-6-specific received interactions (magnitude rank < 0.25, sorted by specificity):")
    print(spec.head(15)[["source", "ligand_complex", "receptor_complex", "magnitude_rank", "specificity_rank"]]
          .round(4).to_string(index=False))
    cx = res[(res["source"] == M6) & res["ligand_complex"].isin(["CXCL9", "CXCL10", "CXCL11"]) &
             (res["receptor_complex"] == "CXCR3")].sort_values("magnitude_rank")
    print("\nM-6 CXCL9/10/11 -> CXCR3 (best targets):")
    print(cx.head(6)[["target", "ligand_complex", "magnitude_rank", "specificity_rank"]].round(4).to_string(index=False))

    # IFNG -> IFNGR axis: which senders, ranked
    ifng = res[(res["ligand_complex"] == "IFNG") & (res["target"] == M6)].sort_values("magnitude_rank")
    print("\nIFNG -> M-6 by sender:")
    print(ifng.head(10)[["source", "receptor_complex", "lr_means", "magnitude_rank", "specificity_rank"]].round(4).to_string(index=False))

    # M-6 chemokine output to lymphocytes
    chemo = res[(res["source"] == M6) & res["ligand_complex"].isin(["CXCL9", "CXCL10", "CXCL11"])].sort_values("magnitude_rank")
    print("\nM-6 CXCL9/10/11 output (top targets):")
    print(chemo.head(10)[["target", "ligand_complex", "receptor_complex", "magnitude_rank"]].round(4).to_string(index=False))

    # IFNG expression by lymphocyte state vs RA GWAS enrichment (step 10, 25 kb, within-lineage)
    lymph = adata.obs["author_ct"].isin(["Tcell", "NK"])
    expr = sc.get.obs_df(adata[lymph], ["IFNG"], use_raw=True)
    expr["state"] = adata.obs.loc[lymph, "state"].astype(str).values
    senders = expr.groupby("state").agg(mean_IFNG=("IFNG", "mean"), pct_IFNG=("IFNG", lambda s: 100 * (s > 0).mean()),
                                        n=("IFNG", "size"))
    gw = pd.read_csv(TAB / "gwas_enrichment.tsv", sep="\t")
    gw = gw[(gw["trait"] == "RA") & (gw["window_kb"] == 25) & (gw["level"] == "subtype")].set_index("state")
    senders = senders.join(gw[["enrichment", "fdr"]].rename(columns={"enrichment": "RA_gwas_enrichment", "fdr": "RA_gwas_fdr"}))
    senders = senders.sort_values("mean_IFNG", ascending=False)
    senders.round(4).to_csv(TAB / "ifng_senders.tsv", sep="\t")
    ok = senders["RA_gwas_enrichment"].notna()
    r_g, p_g = spearmanr(senders.loc[ok, "mean_IFNG"], senders.loc[ok, "RA_gwas_enrichment"])
    print(f"\nIFNG expression vs RA GWAS enrichment across {ok.sum()} T/NK states: rho = {r_g:.2f}, p = {p_g:.3f}")
    print(senders.head(8).round(3).to_string())

    # Donor-level: lymphocyte IFNG vs macrophage IFN-gamma programme
    hall = {}
    for line in open(next((RAW / "genesets").glob("h.all*.gmt"))):
        f = line.rstrip("\n").split("\t")
        hall[f[0]] = f[2:]
    ifn_genes = [g for g in hall["HALLMARK_INTERFERON_GAMMA_RESPONSE"] if g in adata.raw.var_names and g != "IFNG"]
    full = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    sc_code = full.obs["author_subtype"].astype(str).str.split(":").str[0]
    mac = (full.obs["author_ct"] == "myeloid") & ~sc_code.isin(DC_LIKE)
    lym = full.obs["author_ct"].isin(["Tcell", "NK"])
    rows = []
    for d in sorted(full.obs["donor"].unique()):
        md, ld = mac & (full.obs["donor"] == d), lym & (full.obs["donor"] == d)
        if md.sum() < 50 or ld.sum() < 50:
            continue
        cm = np.asarray(full[md].layers["counts"].sum(axis=0)).ravel()
        cl = np.asarray(full[ld].layers["counts"].sum(axis=0)).ravel()
        cpm_m = pd.Series(np.log2(cm / cm.sum() * 1e6 + 1), index=full.var_names)
        cpm_l = pd.Series(np.log2(cl / cl.sum() * 1e6 + 1), index=full.var_names)
        rows.append(dict(donor=d, disease=full.obs.loc[full.obs.donor == d, "disease"].iloc[0],
                         n_mac=int(md.sum()), n_lymph=int(ld.sum()), lymph_IFNG_logCPM=cpm_l["IFNG"],
                         mac_frac_M6=float((full.obs.loc[md, "author_subtype"] == M6).mean()),
                         _mac_vec=cpm_m[ifn_genes]))
    dd = pd.DataFrame(rows).set_index("donor")
    mat = pd.DataFrame(list(dd.pop("_mac_vec")), index=dd.index)
    dd["mac_IFNg_programme"] = ((mat - mat.mean()) / mat.std().replace(0, np.nan)).mean(axis=1)
    dd.round(4).to_csv(TAB / "ifng_donor_link.tsv", sep="\t")
    ra = dd[dd["disease"] == "RA"]
    r1, p1 = spearmanr(ra["lymph_IFNG_logCPM"], ra["mac_IFNg_programme"])
    r2, p2 = spearmanr(ra["lymph_IFNG_logCPM"], ra["mac_frac_M6"])
    print(f"\nDonor level (RA, >= 50 macrophages and lymphocytes, n = {len(ra)}): lymphocyte IFNG vs macrophage "
          f"IFN-gamma programme rho = {r1:.2f} (p = {p1:.3f}); vs M-6 share of macrophages rho = {r2:.2f} (p = {p2:.3f})")
    pd.DataFrame([dict(test="lymphocyte IFNG vs macrophage IFN-gamma programme", n=len(ra), rho=r1, p=p1),
                  dict(test="lymphocyte IFNG vs M-6 share of macrophages", n=len(ra), rho=r2, p=p2),
                  dict(test="IFNG expression vs RA GWAS enrichment across T/NK states", n=int(ok.sum()), rho=r_g, p=p_g)]) \
        .round(4).to_csv(TAB / "cell_communication_tests.tsv", sep="\t", index=False)

    # Figure
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), gridspec_kw={"width_ratios": [1.3, 1, 0.9]})
    ax = axes[0]
    top = spec.head(20).iloc[::-1]
    ax.barh(range(len(top)), -np.log10(top["magnitude_rank"].clip(lower=1e-4)),
            color=["C4" if v <= 1e-3 else "C7" for v in top["specificity_rank"]])
    ax.set_yticks(range(len(top)), [f"{s.split(':')[0]}  {l} -> {r}" for s, l, r in
                                     zip(top["source"], top["ligand_complex"], top["receptor_complex"])], fontsize=7)
    ax.set_xlabel("-log10 LIANA magnitude rank (purple = specificity rank <= 0.001)")
    ax.set_title("M-6-specific signals received by STAT1+ CXCL10+ macrophages\n(magnitude rank < 0.25, sorted by specificity then magnitude)", fontsize=9)
    ax = axes[1]
    s = senders[ok]
    ax.scatter(s["RA_gwas_enrichment"], s["mean_IFNG"], s=20 + s["n"] / 20, c=["C3" if f < 0.1 else "C7" for f in s["RA_gwas_fdr"]])
    for name, rr in s.sort_values("mean_IFNG", ascending=False).head(6).iterrows():
        ax.annotate(name.split(":")[0] + ":" + name.split(":")[1][:16], (rr["RA_gwas_enrichment"], rr["mean_IFNG"]),
                    fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("RA GWAS locus enrichment (within lineage, 25 kb)"); ax.set_ylabel("Mean IFNG (log-normalised)")
    ax.set_title(f"T/NK states: IFNG vs RA genetics, rho = {r_g:.2f}, p = {p_g:.3f}\n(red = GWAS FDR < 0.1)", fontsize=9)
    ax = axes[2]
    ax.scatter(ra["lymph_IFNG_logCPM"], ra["mac_IFNg_programme"], s=40, color="C3")
    for name, rr in ra.iterrows():
        ax.annotate(name.replace("BRI-", ""), (rr["lymph_IFNG_logCPM"], rr["mac_IFNg_programme"]), fontsize=7,
                    xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("Lymphocyte IFNG (log CPM, donor pseudobulk)"); ax.set_ylabel("Macrophage IFN-gamma programme (z)")
    ax.set_title(f"Donor level, RA n = {len(ra)}: rho = {r1:.2f}, p = {p1:.3f}", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "cell_communication.png", dpi=150)


if __name__ == "__main__":
    main()
