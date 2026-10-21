"""Visium validation (GSE347797; Sisó et al. 2026): IFN-activated macrophage programme across donors and within tissue.

12 FFPE sections: 3 healthy (NH), 6 RA, 3 OA blocks from one donor. Spot scores = mean z (within section) of log-CPM:
  IFN programme: CXCL9, CXCL10, CXCL11, STAT1, GBP1, GBP5, IDO1, IRF1, CD274, APOL6
  lymphoid aggregate: CD3E, MS4A1, CXCL13, CCL19, LTB, CD79A, CR2
Hypotheses fixed before analysis:
  H1  donor level: share of spots above the pooled NH 95th percentile of the IFN score, RA (6) vs NH (3), Mann-Whitney.
  H2  within each RA section: partial Spearman of IFN score with (a) lymphoid aggregate score and (b) IFNG in the spot
      and its 6 neighbours, adjusting for log UMI per spot; specificity vs 200 random expression-matched gene sets.

Output: results/tables/visium_sections.tsv, visium_within_section.tsv, visium_tests.tsv; results/figures/visium.png
"""
import glob
import gzip
import io
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import tarfile
from scipy.stats import mannwhitneyu, spearmanr, rankdata

from config import RAW, FIG, TAB, SEED

D = RAW / "GSE347797"
IFN = ["CXCL9", "CXCL10", "CXCL11", "STAT1", "GBP1", "GBP5", "IDO1", "IRF1", "CD274", "APOL6"]
LYM = ["CD3E", "MS4A1", "CXCL13", "CCL19", "LTB", "CD79A", "CR2"]


def read_member(tf, name):
    return gzip.decompress(tf.extractfile(name).read())


def sections(tf):
    names = tf.getnames()
    stems = sorted({n[:-len("_matrix.mtx.gz")] for n in names if n.endswith("_matrix.mtx.gz")})
    out = []
    for s in stems:
        m = re.search(r"(CT_\d+_[A-Z]-\d+)-\d+_(Normal|OA|RA)_Synovium", s)
        out.append(dict(stem=s, block=m.group(1), condition={"Normal": "NH"}.get(m.group(2), m.group(2))))
    return out


def load(tf, stem):
    X = scipy.io.mmread(io.BytesIO(read_member(tf, stem + "_matrix.mtx.gz"))).T.tocsr().astype(np.float32)
    feats = pd.read_csv(io.BytesIO(read_member(tf, stem + "_features.tsv.gz")), sep="\t", header=None)
    bcs = pd.read_csv(io.BytesIO(read_member(tf, stem + "_barcodes.tsv.gz")), header=None)[0].values
    pos = pd.read_csv(io.BytesIO(read_member(tf, stem + "_tissue_positions_list.csv.gz")), header=None)
    if not str(pos.iloc[0, 1]).isdigit():                      # newer Space Ranger files carry a header row
        pos = pos.iloc[1:]
    pos = pos.iloc[:, :4]
    pos.columns = ["barcode", "in_tissue", "row", "col"]
    pos = pos.set_index("barcode").astype(int)
    genes = feats[1].values
    keep = pos.loc[bcs, "in_tissue"].values == 1
    X, bcs = X[keep], bcs[keep]
    return X, genes, pos.loc[bcs, ["row", "col"]].values


def zscore_cols(M):
    return (M - M.mean(0)) / (M.std(0) + 1e-9)


def partial_spearman(x, y, covar):
    rx, ry, rc = rankdata(x), rankdata(y), rankdata(covar)
    def resid(a):
        A = np.column_stack([np.ones_like(rc), rc])
        return a - A @ np.linalg.lstsq(A, a, rcond=None)[0]
    return np.corrcoef(resid(rx), resid(ry))[0, 1]


def main():
    rng = np.random.default_rng(SEED)
    tf = tarfile.open(D / "GSE347797_RAW.tar")
    secs = sections(tf)
    data = []
    for s in secs:
        X, genes, rc = load(tf, s["stem"])
        lib = np.asarray(X.sum(axis=1)).ravel()
        ok = lib >= 200
        X, rc, lib = X[ok], rc[ok], lib[ok]
        L = np.log1p(X.multiply(1e4 / lib[:, None]).tocsr().toarray())
        gi = {g: i for i, g in enumerate(genes)}
        data.append(dict(**s, L=L, gi=gi, rc=rc, loglib=np.log(lib), genes=genes))
        print(f"{s['block']} {s['condition']}: {L.shape[0]} spots")

    def score(d, gl):
        idx = [d["gi"][g] for g in gl if g in d["gi"]]
        return zscore_cols(d["L"][:, idx]).mean(1) if idx else np.full(d["L"].shape[0], np.nan)

    # H1: pooled-NH threshold on the raw (not within-section z) IFN mean log expression
    def raw_ifn(d):
        idx = [d["gi"][g] for g in IFN if g in d["gi"]]
        return d["L"][:, idx].mean(1)
    thr = np.quantile(np.concatenate([raw_ifn(d) for d in data if d["condition"] == "NH"]), 0.95)
    rows = []
    for d in data:
        r = raw_ifn(d)
        rows.append(dict(block=d["block"], condition=d["condition"], n_spots=len(r), mean_ifn=float(r.mean()),
                         pct_spots_ifn_high=100 * float((r > thr).mean())))
    sec = pd.DataFrame(rows)
    sec.to_csv(TAB / "visium_sections.tsv", sep="\t", index=False, float_format="%.4g")
    ra, nh = sec[sec.condition == "RA"], sec[sec.condition == "NH"]
    u = mannwhitneyu(ra.pct_spots_ifn_high, nh.pct_spots_ifn_high)
    print(sec.round(3).to_string(index=False))
    print(f"H1: IFN-high spot share RA median {ra.pct_spots_ifn_high.median():.1f}% vs NH {nh.pct_spots_ifn_high.median():.1f}% "
          f"(6 vs 3, Mann-Whitney p = {u.pvalue:.3f})")

    # H2 within RA sections
    within = []
    for d in data:
        if d["condition"] != "RA":
            continue
        ifn = score(d, IFN); lym = score(d, LYM)
        g_ifng = d["gi"].get("IFNG")
        if g_ifng is not None:
            ig = d["L"][:, g_ifng]
            key = {tuple(x): i for i, x in enumerate(d["rc"])}
            nbr = np.array([ig[i] + sum(ig[key[(r + dr, c + dc)]] for dr, dc in ((0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1))
                                        if (r + dr, c + dc) in key) for i, (r, c) in enumerate(d["rc"])])
        else:
            nbr = np.full(len(ifn), np.nan)
        r_lym = partial_spearman(ifn, lym, d["loglib"])
        r_ifng = partial_spearman(ifn, nbr, d["loglib"]) if np.isfinite(nbr).all() else np.nan
        # specificity: random gene sets matched on mean expression
        mean_expr = d["L"].mean(0)
        bins = pd.qcut(pd.Series(mean_expr).rank(method="first"), 20, labels=False).values
        want = [bins[d["gi"][g]] for g in IFN if g in d["gi"]]
        null_l, null_i = [], []
        for _ in range(200):
            idx = [rng.choice(np.where(bins == q)[0]) for q in want]
            rs = zscore_cols(d["L"][:, idx]).mean(1)
            null_l.append(partial_spearman(rs, lym, d["loglib"]))
            if np.isfinite(nbr).all():
                null_i.append(partial_spearman(rs, nbr, d["loglib"]))
        within.append(dict(block=d["block"], n_spots=len(ifn), rho_ifn_vs_lymphoid=r_lym,
                           random_median_lymphoid=float(np.median(null_l)), pct_random_ge_lymphoid=float(np.mean(np.array(null_l) >= r_lym)),
                           rho_ifn_vs_IFNG_neighbourhood=r_ifng,
                           random_median_IFNG=float(np.median(null_i)) if null_i else np.nan,
                           pct_random_ge_IFNG=float(np.mean(np.array(null_i) >= r_ifng)) if null_i else np.nan))
    w = pd.DataFrame(within)
    w.to_csv(TAB / "visium_within_section.tsv", sep="\t", index=False, float_format="%.4g")
    print(w.round(3).to_string(index=False))
    pd.DataFrame([dict(test="H1 IFN-high spot share, RA vs NH (donor level)", n=f"{len(ra)} vs {len(nh)}", p=u.pvalue,
                       note=f"RA median {ra.pct_spots_ifn_high.median():.2f}%, NH median {nh.pct_spots_ifn_high.median():.2f}%"),
                  dict(test="H2a IFN vs lymphoid aggregate score (partial on log UMI), RA sections beating 95% of random sets",
                       n=len(w), p=np.nan, note=f"{int((w.pct_random_ge_lymphoid < 0.05).sum())}/{len(w)}"),
                  dict(test="H2b IFN vs IFNG in spot + 6 neighbours, RA sections beating 95% of random sets",
                       n=len(w), p=np.nan, note=f"{int((w.pct_random_ge_IFNG < 0.05).sum())}/{len(w)}")]) \
        .to_csv(TAB / "visium_tests.tsv", sep="\t", index=False, float_format="%.4g")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    ax = axes[0]
    for i, cond in enumerate(["NH", "OA", "RA"]):
        v = sec.loc[sec.condition == cond, "pct_spots_ifn_high"]
        ax.scatter(np.full(len(v), i) + rng.uniform(-.08, .08, len(v)), v, s=45)
    ax.set_xticks(range(3), ["Healthy", "OA (1 donor)", "RA"]); ax.set_ylabel("% spots above healthy 95th pct IFN score")
    ax.set_title(f"H1: RA vs healthy p = {u.pvalue:.3f}", fontsize=10)
    big = max((d for d in data if d["condition"] == "RA"), key=lambda d: d["L"].shape[0])
    for ax, vals, title in ((axes[1], score(big, IFN), "IFN macrophage programme"), (axes[2], score(big, LYM), "Lymphoid aggregate score")):
        sca = ax.scatter(big["rc"][:, 1], -big["rc"][:, 0], c=vals, s=6, cmap="magma", vmin=-1, vmax=2)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{big['block']} (RA): {title}", fontsize=10); plt.colorbar(sca, ax=ax, shrink=0.6)
    fig.tight_layout(); fig.savefig(FIG / "visium.png", dpi=150)


if __name__ == "__main__":
    main()
