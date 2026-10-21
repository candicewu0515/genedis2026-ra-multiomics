"""TF activity from RNA (CollecTRI regulons) vs TF motif accessibility from chromatin (chromVAR, step 14).

Units: the same subtype x donor pseudobulk units as step 14 (>= 10 nuclei) plus per-donor macrophage pools.
RNA activity: genes z-scored across units (log CPM, genes with mean log CPM > 1); per unit and TF, a univariate
linear model of gene z on the TF's signed CollecTRI targets (t statistic; the ULM of decoupler, Badia-i-Mompel 2022);
TFs with >= 10 measured targets. CollecTRI from OmniPath (cached); sign +1 stimulation, -1 inhibition.
Matching to chromVAR: non-dimer JASPAR motif names equal to the TF symbol (several motif versions averaged).
Questions fixed before analysis:
  T1 sanity: lineage TFs top the RNA activity of their cell type (SPI1/CEBP myeloid, PAX5/EBF1 B, TCF7/LEF1 T ...).
  T2 concordance: within macrophage subtype units (so lineage does not drive it), per-TF Spearman between RNA activity
     and motif deviation; matched TF pairs vs mismatched pairs from different TF families (Mann-Whitney).
  T3 M-6 candidate regulators: TFs more active in M-6 STAT1+ CXCL10+ than in the other macrophage subtypes of the same donor on
     RNA (vs 500 gene-label permutations, which keep regulon sizes and overlaps) and on chromatin (step 14);
     which TFs agree on both layers; where IRF5 (the colocalising RA gene, sections 20 and 29) falls.

Output: data/raw/collectri/collectri_omnipath.tsv, results/tables/tf_activity_units.tsv.gz, tf_crossmodal_concordance.tsv,
        tf_M6_drivers.tsv; results/figures/tf_activity_crossmodal.png
"""
import io
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import mannwhitneyu, spearmanr, rankdata
from statsmodels.stats.multitest import multipletests

from config import RAW, INTERIM, FIG, TAB, SEED

MIN_UNIT = 10
MIN_TARGETS = 10
N_PERM = 500
M6 = "M-6: STAT1+ CXCL10+"
DC_LIKE = ["M-9", "M-10", "M-11", "M-12", "M-13", "M-14"]
LINEAGE = {"myeloid": ["SPI1", "CEBPA", "CEBPB", "IRF8", "MAFB"], "Bcell": ["PAX5", "EBF1", "POU2AF1"],
           "Tcell": ["TCF7", "LEF1", "GATA3", "TBX21", "RUNX3"], "plasma": ["PRDM1", "XBP1", "IRF4"],
           "fibroblast": ["TWIST1", "TWIST2", "PRRX1"], "endothelial": ["ERG", "FLI1", "SOX17"]}
C = RAW / "collectri"
C.mkdir(exist_ok=True)


def collectri():
    p = C / "collectri_omnipath.tsv"
    if not p.exists():
        url = "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&fields=sources,references"
        with urllib.request.urlopen(url, timeout=300) as r:
            p.write_bytes(r.read())
    d = pd.read_csv(p, sep="\t")
    d["weight"] = np.where(d["consensus_inhibition"].astype(bool) & ~d["consensus_stimulation"].astype(bool), -1.0, 1.0)
    d = d.drop(columns=["source", "target"]).rename(columns={"source_genesymbol": "tf", "target_genesymbol": "target"})
    return d[["tf", "target", "weight"]].drop_duplicates(["tf", "target"])


def units_rna(rna):
    obs = rna.obs
    sub = obs["author_subtype"].astype(str)
    code = sub.str.split(":").str[0]
    lab = (sub + "|" + obs["donor"].astype(str)).values
    vc = pd.Series(lab).value_counts()
    keep = list(vc[vc >= MIN_UNIT].index)
    mac = (obs["author_ct"] == "myeloid").values & ~code.isin(DC_LIKE).values
    groups = [(g, np.where(lab == g)[0]) for g in keep]
    groups += [(f"MACROPHAGES|{d}", np.where(mac & (obs["donor"].values == d))[0]) for d in sorted(obs["donor"].unique())]
    groups = [(g, i) for g, i in groups if len(i) >= MIN_UNIT]
    Cn = rna.layers["counts"].tocsr()
    M = sp.csr_matrix((np.ones(sum(len(i) for _, i in groups)),
                       (np.concatenate([[k] * len(i) for k, (_, i) in enumerate(groups)]), np.concatenate([i for _, i in groups]))),
                      shape=(len(groups), Cn.shape[0]))
    X = np.asarray((M @ Cn).todense(), dtype=np.float64)
    ct_of = obs.groupby(sub.values)["author_ct"].first()
    meta = pd.DataFrame({"unit": [g for g, _ in groups], "n": [len(i) for _, i in groups]}).set_index("unit")
    meta["state"] = meta.index.str.split("|").str[0]
    meta["donor"] = meta.index.str.split("|").str[1]
    meta["cell_type"] = [("macrophage_pool" if s == "MACROPHAGES" else ct_of.get(s, "?")) for s in meta["state"]]
    return pd.DataFrame(X, index=meta.index, columns=rna.var_names), meta


def ulm(Y, W):
    """t statistic of the slope of y ~ w per unit (columns of Y) and TF (columns of W); Y genes x units."""
    n = Y.shape[0]
    Wc = W - W.mean(0); Yc = Y - Y.mean(0)
    sxx = (Wc ** 2).sum(0)                                   # TFs
    b = (Wc.T @ Yc) / sxx[:, None]                            # TFs x units
    syy = (Yc ** 2).sum(0)                                    # units
    rss = syy[None, :] - b ** 2 * sxx[:, None]
    se = np.sqrt(np.maximum(rss, 1e-12) / (n - 2) / sxx[:, None])
    return b / se


def main():
    rng = np.random.default_rng(SEED)
    rna = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    X, meta = units_rna(rna)
    del rna
    cpm = X.div(X.sum(axis=1), axis=0) * 1e6
    L = np.log2(cpm + 1)
    L = L.loc[:, L.mean(axis=0) > 1]
    Z = (L - L.mean(0)) / (L.std(0) + 1e-9)
    net = collectri()
    net = net[net.target.isin(Z.columns)]
    nt = net.groupby("tf").size()
    tfs = nt[nt >= MIN_TARGETS].index
    net = net[net.tf.isin(tfs)]
    Wm = net.pivot_table(index="target", columns="tf", values="weight", fill_value=0).reindex(Z.columns, fill_value=0)
    print(f"Units {len(meta)}; genes {Z.shape[1]}; CollecTRI TFs with >= {MIN_TARGETS} targets: {Wm.shape[1]}")
    A = pd.DataFrame(ulm(Z.values.T, Wm.values).T, index=Z.index, columns=Wm.columns)
    A.round(4).to_csv(TAB / "tf_activity_units.tsv.gz", sep="\t")

    # T1 lineage sanity
    print("\nT1 top RNA-inferred TFs per cell type (mean unit activity vs other cell types):")
    sanity = []
    for ct, expect in LINEAGE.items():
        u = meta.index[meta.cell_type == ct]
        if not len(u):
            continue
        d = A.loc[u].mean() - A.loc[meta.index.difference(u)].mean()
        rk = d.rank(ascending=False)
        top = d.sort_values(ascending=False).head(8).index.tolist()
        hits = [t for t in expect if t in rk.index]
        sanity.append(dict(cell_type=ct, top8=", ".join(top), expected=", ".join(f"{t} #{int(rk[t])}" for t in hits),
                           expected_in_top20=sum(rk[t] <= 20 for t in hits), n_expected=len(hits)))
        print(f"  {ct}: {', '.join(top)} | expected: " + ", ".join(f"{t} #{int(rk[t])}" for t in hits))

    # chromVAR deviations on the same units
    cv = pd.read_csv(TAB / "chromvar_unit_deviations.tsv.gz", sep="\t", index_col=0)
    names = pd.Series([c.split("|")[0] for c in cv.columns], index=cv.columns)
    single = names[~names.str.contains("::")].str.upper()
    mot = cv[single.index].T.groupby(single.values).mean().T            # units x TF (motif)
    common_units = A.index.intersection(mot.index)
    shared = sorted(set(A.columns) & set(mot.columns))
    print(f"\nUnits shared with chromVAR: {len(common_units)}; TFs with both RNA activity and a motif: {len(shared)}")

    # T2 concordance within macrophage subtype units
    is_mac = meta.cell_type.eq("myeloid") & ~meta.state.str.split(":").str[0].isin(DC_LIKE)
    mu_ = [u for u in common_units if is_mac.get(u, False)]
    a, m = A.loc[mu_, shared], mot.loc[mu_, shared]
    ra = pd.DataFrame(rankdata(a, axis=0), index=a.index, columns=shared)
    rm = pd.DataFrame(rankdata(m, axis=0), index=m.index, columns=shared)
    zc = lambda df: (df - df.mean()) / (df.std() + 1e-12)
    R = zc(ra).T @ zc(rm) / (len(mu_) - 1)                                # TF (RNA) x TF (motif) Spearman
    matched = pd.Series(np.diag(R.values), index=shared)
    fam = lambda t: t[:3]
    mism = [R.values[i, j] for i in range(len(shared)) for j in range(len(shared)) if fam(shared[i]) != fam(shared[j])]
    u = mannwhitneyu(matched, mism, alternative="greater")
    print(f"T2 within {len(mu_)} macrophage units: matched TF median rho {matched.median():.3f} vs mismatched (other "
          f"families) {np.median(mism):.3f}; Mann-Whitney one-sided p {u.pvalue:.2g}")
    conc = pd.DataFrame({"rho_rna_activity_vs_motif": matched}).sort_values("rho_rna_activity_vs_motif", ascending=False)
    conc.round(4).to_csv(TAB / "tf_crossmodal_concordance.tsv", sep="\t")
    print("  most concordant:", ", ".join(f"{t} {v:.2f}" for t, v in conc.rho_rna_activity_vs_motif.head(12).items()))

    # T3 M-6 drivers (paired within donors that have an M-6 unit)
    def m6_diff(Adf):
        diffs = []
        for d in sorted(meta.donor.unique()):
            m6u = meta.index[(meta.state == M6) & (meta.donor == d)]
            oth = meta.index[is_mac & (meta.state != M6) & (meta.donor == d)]
            if len(m6u) and len(oth):
                diffs.append(Adf.loc[m6u].mean() - Adf.loc[oth].mean())
        return pd.DataFrame(diffs)
    D = m6_diff(A)
    obs = D.mean()
    mac_units = meta.index[is_mac]
    Zm = Z.loc[mac_units].values.T                                          # genes x mac units
    null = np.zeros((N_PERM, A.shape[1]))
    for b in range(N_PERM):
        perm = rng.permutation(Zm.shape[0])
        Ap = pd.DataFrame(ulm(Zm[perm], Wm.values).T, index=mac_units, columns=Wm.columns)
        null[b] = m6_diff(Ap).mean().values
    p_perm = (1 + (null >= obs.values[None, :]).sum(0)) / (1 + N_PERM)
    cvm = pd.read_csv(TAB / "chromvar_M6_vs_macrophages.tsv", sep="\t", index_col=0)
    cvm = cvm[~cvm.motif.str.contains("::")]
    motif_diff = cvm.groupby(cvm.motif.str.upper())["mean_diff_z"].mean()
    drv = pd.DataFrame({"rna_activity_diff": obs, "rna_p_perm": p_perm, "n_donors": len(D),
                        "n_targets": nt.reindex(obs.index).values}).join(motif_diff.rename("motif_diff_z"))
    drv["rna_rank"] = drv.rna_activity_diff.rank(ascending=False)
    drv["motif_rank"] = drv.motif_diff_z.rank(ascending=False)
    # Submission revision: adjust all tested RNA TF permutation P values before candidate selection.
    # This corrects multiplicity within this gene-label null, not donor-level biological uncertainty.
    drv["rna_q_BH"] = multipletests(drv.rna_p_perm, method="fdr_bh")[1]
    motif_candidate=(drv.motif_diff_z > 0) & (drv.motif_rank <= 0.1 * drv.motif_diff_z.notna().sum())
    drv["both_layers_nominal"] = (drv.rna_p_perm < 0.05) & motif_candidate
    drv["both_layers"] = (drv.rna_q_BH < 0.05) & motif_candidate
    drv["complex"] = drv.index.str.contains("_")                        # CollecTRI dimers share their members' regulons
    drv = drv.sort_values("rna_activity_diff", ascending=False)
    drv.to_csv(TAB / "tf_M6_drivers.tsv", sep="\t", float_format="%.9g")
    pd.set_option("display.width", 200)
    print(f"\nT3 M-6 vs other macrophages ({len(D)} donors); top RNA-inferred TFs:")
    print(drv[~drv.complex].head(15).round(3).to_string())
    print("Candidate TFs on both layers (RNA gene-label BH q < 0.05 and motif in the top 10 %):", ", ".join(drv.index[drv.both_layers]))
    for t in ("IRF5", "STAT1", "IRF1", "IRF8", "NFKB1", "RELA", "SPI1"):
        if t in drv.index:
            r = drv.loc[t]
            print(f"  {t}: RNA diff {r.rna_activity_diff:.2f} (rank {int(r.rna_rank)}, perm p {r.rna_p_perm:.3f}); "
                  f"motif diff {r.motif_diff_z if pd.notna(r.motif_diff_z) else float('nan'):.2f}")
    pd.DataFrame(sanity).to_csv(TAB / "tf_lineage_sanity.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.5))
    ax = axes[0]
    ax.hist(mism, bins=60, density=True, alpha=0.5, label="mismatched TF pairs (other families)")
    ax.hist(matched, bins=30, density=True, alpha=0.6, label="matched TF (RNA activity vs own motif)")
    ax.axvline(0, color="k", lw=0.8); ax.legend(fontsize=8)
    ax.set_xlabel("Spearman across macrophage units"); ax.set_title(f"T2 concordance (p {u.pvalue:.2g})", fontsize=10)
    ax = axes[1]
    dd = drv.dropna(subset=["motif_diff_z"])
    ax.scatter(dd.rna_activity_diff, dd.motif_diff_z, s=8, color="0.6")
    lab = dd[dd.both_layers | dd.index.isin(["IRF5", "STAT1", "IRF1", "IRF8", "RELA", "NFKB1"])]
    ax.scatter(lab.rna_activity_diff, lab.motif_diff_z, s=18, color=["C3" if b else "k" for b in lab.both_layers])
    for t, r in lab.iterrows():
        ax.annotate(t, (r.rna_activity_diff, r.motif_diff_z), fontsize=7)
    ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("RNA-inferred activity, M-6 minus other macrophages"); ax.set_ylabel("motif deviation, M-6 minus other (z)")
    ax.set_title(f"T3 M-6 drivers ({len(D)} donors; red: both layers, black: reference TF)", fontsize=10)
    ax = axes[2]
    top = drv[~drv.complex].head(15)
    ax.barh(range(len(top)), top.rna_activity_diff, color=["C3" if b else "C0" for b in top.both_layers])
    ax.set_yticks(range(len(top)), top.index); ax.invert_yaxis()
    ax.set_xlabel("RNA-inferred activity difference (red: also motif top 10 %)"); ax.set_title("Top RNA-inferred M-6 TFs (dimer regulons not shown)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "tf_activity_crossmodal.png", dpi=150)


if __name__ == "__main__":
    main()
