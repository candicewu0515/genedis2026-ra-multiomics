"""The CD4 Tph / Treg axis that carries RA heritability (S-LDSC, step 33), tested layer by layer.

States: Tph = T-7 CD4+ Tph and T-3 CD4+ Tfh/Tph; Treg = T-8 CD25-high and T-9 CD25-low Treg (GSE243917 author labels).
Bulk programmes (a priori): Tph = PDCD1, CXCL13, MAF, TOX, ICOS; Treg = FOXP3, IL2RA, CTLA4, IKZF2, TNFRSF18;
T-cell abundance = CD3E, CD3D, CD2, CD247 (covariate).
Hypotheses fixed before analysis:
  P1 composition (GSE243917, 11 RA donors): Tph and Treg shares of T cells rise with lymphocytic inflammation grade;
     Tph share co-varies with the M-6 STAT1+ CXCL10+ share of macrophages (Spearman; prep batch reported).
  P2 bulk stage (GSE89408): Tph and Treg programmes higher in RA than OA, adjusted for age, sex, depth and T-cell
     abundance, beyond 1,000 expression-matched random gene sets; within RA, Tph tracks the IFN macrophage programme
     (partial on T-cell abundance) beyond random sets.
  P3 communication (LIANA, step 13): CXCL13-CXCR5 and IL21-IL21R from Tph states to B-cell states rank in the top 5 %
     of each sender's interactions by specificity; position of Tph among IFNG senders to M-6.
  P4 spatial (Xenium 5K RA sections): B cells within 30 um are more frequent around CXCL13+ T cells than around
     CXCL13- T cells (CXCL13 labels permuted among T cells, 1,000 times); macrophage IFN score rises with CXCL13+ T
     neighbours beyond CXCL13- T neighbours (as step 20).
  P5 anti-TNF (GSE198520, exploratory): Tph programme change week 12 - baseline in good vs non-responders, against
     random gene sets (the cohort's library-size axis makes raw comparisons unreliable, step 11).

Output: results/tables/tph_treg_tests.tsv, tph_treg_composition.tsv, tph_spatial.tsv; results/figures/tph_treg_axis.png
"""
import importlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, norm, spearmanr, wilcoxon

from config import RAW, INTERIM, FIG, TAB, SEED

step20 = importlib.import_module("20_spatial_xenium")
step30 = importlib.import_module("30_gse89408_disease_stage")
TPH_STATES = ["T-7: CD4+ Tph", "T-3: CD4+ Tfh/Tph"]
TREG_STATES = ["T-8: CD4+ CD25-high Treg", "T-9: CD4+ CD25-low Treg"]
M6 = "M-6: STAT1+ CXCL10+"
DC_LIKE = ["M-9", "M-10", "M-11", "M-12", "M-13", "M-14"]
TPH = ["PDCD1", "CXCL13", "MAF", "TOX", "ICOS"]
TREG = ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"]
TCELL = ["CD3E", "CD3D", "CD2", "CD247"]
IFN = step30.IFN
N_RAND, N_PERM, RADIUS = 1000, 1000, 30.0


def zscore_rows(L):
    return L.sub(L.mean(axis=1), axis=0).div(L.std(axis=1) + 1e-9, axis=0)


def p1_composition(rows):
    rna = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad", backed="r")
    obs = rna.obs[["donor", "disease", "author_ct", "author_subtype", "lymphocytic_inflammation_grade", "prep_batch"]].copy()
    rna.file.close()
    obs["code"] = obs.author_subtype.astype(str).str.split(":").str[0]
    out = []
    for d, g in obs.groupby("donor", observed=True):
        t = g[g.author_ct == "Tcell"]; mac = g[(g.author_ct == "myeloid") & ~g.code.isin(DC_LIKE)]
        out.append(dict(donor=d, disease=g.disease.iloc[0], grade=g.lymphocytic_inflammation_grade.iloc[0],
                        prep_batch=g.prep_batch.iloc[0], n_T=len(t),
                        tph_share_of_T=t.author_subtype.isin(TPH_STATES).mean() if len(t) else np.nan,
                        treg_share_of_T=t.author_subtype.isin(TREG_STATES).mean() if len(t) else np.nan,
                        T_share_of_all=len(t) / len(g),
                        m6_share_of_mac=(mac.author_subtype == M6).mean() if len(mac) else np.nan))
    comp = pd.DataFrame(out)
    comp.to_csv(TAB / "tph_treg_composition.tsv", sep="\t", index=False, float_format="%.4g")
    ra = comp[(comp.disease == "RA") & (comp.n_T >= 50)]
    for col in ("tph_share_of_T", "treg_share_of_T"):
        r, p = spearmanr(ra[col], ra.grade.astype(float))
        rows.append(dict(test=f"P1 {col} vs lymphocytic inflammation grade (RA donors, >= 50 T cells)", n=len(ra), stat=r, p=p))
        b = mannwhitneyu(ra.loc[ra.prep_batch == ra.prep_batch.unique()[0], col], ra.loc[ra.prep_batch != ra.prep_batch.unique()[0], col]).pvalue \
            if ra.prep_batch.nunique() == 2 else np.nan
        rows.append(dict(test=f"P1 {col} vs prep batch (technical check)", n=len(ra), stat=np.nan, p=b))
    r, p = spearmanr(ra.tph_share_of_T, ra.m6_share_of_mac)
    rows.append(dict(test="P1 Tph share of T vs M-6 share of macrophages", n=len(ra), stat=r, p=p))
    print(ra.round(3).to_string(index=False))
    return comp


def p2_bulk(rows, rng):
    meta = step30.metadata()
    X = pd.read_csv(step30.D / "GSE89408_GEO_count_matrix_rename.txt.gz", sep="\t", index_col=0)[meta.index].round()
    X = X[~X.index.duplicated()]
    cpm = X / X.sum(axis=0) * 1e6
    keep = (cpm >= 1).mean(axis=1) >= 0.2
    pos = X[(X > 0).all(axis=1)]; lg = np.log(pos)
    sf = np.exp((lg.sub(lg.mean(axis=1), axis=0)).median(axis=0))
    L = np.log2(X[keep] / sf + 1)
    z = zscore_rows(L)
    sc_ = lambda genes: z.loc[[g for g in genes if g in z.index]].mean(axis=0)
    meta = meta.assign(Tph=sc_(TPH), Treg=sc_(TREG), Tab=sc_(TCELL), IFNp=sc_(IFN), log_lib=np.log(X.sum(axis=0)))
    ra = meta.group.isin(["early RA", "established RA"])
    d2 = meta[ra | (meta.group == "OA")].assign(RA=lambda d: d.group.ne("OA").astype(int)).dropna(subset=["age"])
    dra = meta[ra].dropna(subset=["age"])
    bins = pd.qcut(L.mean(axis=1).rank(method="first"), 20, labels=False)
    members = {q: bins.index[bins == q].values for q in range(20)}

    def fit_t(df, formula, term):
        return smf.ols(formula, data=df).fit(cov_type="HC3").tvalues[term]
    for name, genes in (("Tph", TPH), ("Treg", TREG)):
        want = bins.loc[[g for g in genes if g in bins.index]].values
        f = f"{name} ~ RA + age + female + log_lib + Tab"
        t_obs = fit_t(d2, f, "RA")
        fit = smf.ols(f, data=d2).fit(cov_type="HC3")
        null = []
        for _ in range(N_RAND):
            R = z.loc[[rng.choice(members[q]) for q in want]].mean(axis=0)
            null.append(fit_t(d2.assign(R=R[d2.index].values), "R ~ RA + age + female + log_lib + Tab", "RA"))
        null = np.array(null)
        rows.append(dict(test=f"P2 {name} programme RA vs OA (adjusted incl. T-cell abundance)", n=len(d2), stat=t_obs,
                         p=fit.pvalues["RA"], p_vs_random=(1 + (null >= t_obs).sum()) / (1 + N_RAND), random_median=float(np.median(null))))
    # within RA: Tph vs IFN programme, partial on T-cell abundance
    want = bins.loc[[g for g in TPH if g in bins.index]].values
    f = "IFNp ~ Tph + Tab + age + female + log_lib"
    t_obs = fit_t(dra, f, "Tph")
    null = []
    for _ in range(N_RAND):
        R = z.loc[[rng.choice(members[q]) for q in want]].mean(axis=0)
        null.append(fit_t(dra.assign(R=R[dra.index].values), "IFNp ~ R + Tab + age + female + log_lib", "R"))
    null = np.array(null)
    rows.append(dict(test="P2 within RA: IFN macrophage programme ~ Tph programme (partial on T-cell abundance)",
                     n=len(dra), stat=t_obs, p=2 * norm.sf(abs(t_obs)), p_vs_random=(1 + (null >= t_obs).sum()) / (1 + N_RAND),
                     random_median=float(np.median(null))))
    return meta


def p3_communication(rows):
    li = pd.read_csv(TAB / "liana_results.tsv.gz", sep="\t")
    out = []
    for lig, rec in (("CXCL13", "CXCR5"), ("IL21", "IL21R")):
        for s in TPH_STATES:
            sub = li[li.source == s]
            pct = sub.specificity_rank.rank(pct=True)
            hit = sub[(sub.ligand_complex == lig) & (sub.receptor_complex == rec) & sub.target.str.startswith("B-")]
            if len(hit):
                best = hit.loc[pct[hit.index].idxmin()]
                out.append(dict(pair=f"{lig}-{rec}", source=s, best_B_target=best.target,
                                specificity_percentile=float(pct[best.name]), magnitude_rank=best.magnitude_rank))
    com = pd.DataFrame(out)
    print(com.to_string(index=False))
    for pair, g in com.groupby("pair"):
        rows.append(dict(test=f"P3 {pair} Tph -> B: best specificity percentile among each Tph state's interactions",
                         n=len(g), stat=float(g.specificity_percentile.min()), p=np.nan,
                         note="; ".join(f"{r.source} -> {r.best_B_target} {100 * r.specificity_percentile:.1f} %" for r in g.itertuples())))
    ifng = li[(li.ligand_complex == "IFNG") & (li.target == M6)].sort_values("specificity_rank")
    ifng["rank"] = np.arange(1, len(ifng) + 1)
    tph_rank = ifng[ifng.source.isin(TPH_STATES)][["source", "rank"]]
    rows.append(dict(test="P3 IFNG senders to M-6: rank of Tph states (of all senders)", n=len(ifng), stat=float(tph_rank["rank"].min()) if len(tph_rank) else np.nan,
                     p=np.nan, note="; ".join(f"{r.source} #{r.rank}" for r in tph_rank.itertuples()) + f"; top: {ifng.source.iloc[0]}"))
    return com


def p4_spatial(rows, rng):
    out = []
    for sec in step20.sections():
        if sec["panel"] != "5K" or sec["condition"] != "RA":
            continue
        a = step20.annotate(step20.load(sec), "5K")
        if "CXCL13" not in a.var_names:
            continue
        ct = a.obs.cell_type.values
        T = np.where(ct == "T")[0]; B = np.where(ct == "B")[0]; M = np.where(ct == "macrophage")[0]
        cx = np.asarray(a[:, "CXCL13"].layers["counts"].todense()).ravel() > 0
        lab = cx[T]
        xy = a.obs[["x", "y"]].values
        nbB = np.array([len(i) for i in cKDTree(xy[B]).query_ball_point(xy[T], r=RADIUS)]) if len(B) else np.zeros(len(T))
        diff = lambda l: nbB[l].mean() - nbB[~l].mean() if l.any() and (~l).any() else np.nan
        d_obs = diff(lab)
        null = np.array([diff(rng.permutation(lab)) for _ in range(N_PERM)])
        p_b = (1 + np.sum(null >= d_obs)) / (1 + N_PERM)
        # macrophage IFN score ~ CXCL13+ T neighbours + CXCL13- T neighbours
        genes = [g for g in IFN if g in a.var_names]
        Xm = a[M][:, genes].X
        Xm = np.asarray(Xm.todense()) if hasattr(Xm, "todense") else Xm
        score = ((Xm - Xm.mean(0)) / (Xm.std(0) + 1e-9)).mean(1)
        nbT = cKDTree(xy[T]).query_ball_point(xy[M], r=RADIUS)

        def coef(l):
            pos = np.array([l[i].sum() for i in nbT]); neg = np.array([len(i) - l[i].sum() for i in nbT])
            return sm.OLS(score, sm.add_constant(np.column_stack([pos, neg]))).fit().params[1]
        c_obs = coef(lab)
        cnull = np.array([coef(rng.permutation(lab)) for _ in range(N_PERM)])
        p_m = (1 + np.sum(cnull >= c_obs)) / (1 + N_PERM)
        out.append(dict(block=sec["block"], n_T=len(T), pct_T_CXCL13pos=100 * lab.mean(), n_B=len(B), n_macrophages=len(M),
                        B_neighbours_CXCL13pos_T=nbB[lab].mean(), B_neighbours_CXCL13neg_T=nbB[~lab].mean(),
                        z_B=(d_obs - null.mean()) / (null.std() + 1e-12), p_B=p_b,
                        coef_mac_IFN_per_CXCL13pos_T=c_obs, z_mac=(c_obs - cnull.mean()) / (cnull.std() + 1e-12), p_mac=p_m))
        print(f"  {sec['block']}: T {len(T)} ({100 * lab.mean():.1f} % CXCL13+), B {len(B)}; B-neighbour z {out[-1]['z_B']:.2f} "
              f"(p {p_b:.3f}); macrophage IFN z {out[-1]['z_mac']:.2f} (p {p_m:.3f})", flush=True)
    sp = pd.DataFrame(out)
    sp.to_csv(TAB / "tph_spatial.tsv", sep="\t", index=False, float_format="%.4g")
    for key, lab_ in (("B", "P4 B cells within 30 um of CXCL13+ vs CXCL13- T cells"), ("mac", "P4 macrophage IFN score ~ CXCL13+ T neighbours")):
        if key=="B":
            positive=sp.B_neighbours_CXCL13pos_T > sp.B_neighbours_CXCL13neg_T
        else:
            positive=sp.coef_mac_IFN_per_CXCL13pos_T > 0
        rows.append(dict(test=lab_ + " (section-wise empirical P; no combined tail)", n=len(sp),
                         stat=int(positive.sum()), p=np.nan,
                         note=f"{int(positive.sum())}/{len(sp)} positive observed effects; {int((sp[f'p_{key}'] < 0.05).sum())}/{len(sp)} sections empirical p < 0.05; stat is positive-effect count; global label shuffling is exploratory"))
    return sp


def p5_antitnf(rows, rng):
    meta = pd.read_csv(RAW.parent / "processed" / "gse198520_sample_scores.tsv", sep="\t", index_col=0)
    X = pd.read_csv(RAW / "GSE198520" / "GSE198520_Raw_gene_count_matrix.txt.gz", sep="\t", index_col=0)
    X = X[~X.index.str.startswith("MT-")]
    X = X[[c for c in meta.index if c in X.columns]]
    cpm = X / X.sum(axis=0) * 1e6
    L = np.log2(cpm[(cpm >= 1).mean(axis=1) >= 0.2] + 1)
    z = zscore_rows(L)
    bins = pd.qcut(L.mean(axis=1).rank(method="first"), 20, labels=False)
    members = {q: bins.index[bins == q].values for q in range(20)}
    m = meta.loc[z.columns]

    def stat(score):
        w = pd.DataFrame({"s": score, "patient": m.patient, "tp": m.timepoint, "resp": m.response})
        w = w.pivot_table(index=["patient"], columns="tp", values="s").join(m.groupby("patient").response.first())
        w = w.dropna(subset=["Pre", "Post"]); w["d"] = w.Post - w.Pre
        g, n = w[w.response == "good"].d, w[w.response == "non"].d
        return g.median() - n.median(), g, n
    genes = [g for g in TPH if g in z.index]
    s_obs, g, n = stat(z.loc[genes].mean(axis=0))
    want = bins.loc[genes].values
    null = np.array([stat(z.loc[[rng.choice(members[q]) for q in want]].mean(axis=0))[0] for _ in range(N_RAND)])
    rows.append(dict(test="P5 Tph programme change (week 12 - baseline), good minus non-responders (median)", n=len(g) + len(n),
                     stat=s_obs, p=mannwhitneyu(g, n).pvalue, p_vs_random=(1 + (null <= s_obs).sum()) / (1 + N_RAND),
                     random_median=float(np.median(null)), note=f"good median change {g.median():.2f} (Wilcoxon p {wilcoxon(g).pvalue:.3g}); non {n.median():.2f}"))


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    comp = p1_composition(rows)
    bulk = p2_bulk(rows, rng)
    com = p3_communication(rows)
    sp = p4_spatial(rows, rng)
    p5_antitnf(rows, rng)
    res = pd.DataFrame(rows)
    res.to_csv(TAB / "tph_treg_tests.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 110)
    print(res.round(4).to_string(index=False))

    fig, axes = plt.subplots(1, 4, figsize=(22, 5))
    ra = comp[(comp.disease == "RA") & (comp.n_T >= 50)]
    ax = axes[0]
    ax.scatter(ra.grade.astype(float) + np.random.default_rng(SEED).uniform(-.1, .1, len(ra)), 100 * ra.tph_share_of_T, label="Tph (T-3, T-7)")
    ax.scatter(ra.grade.astype(float) + np.random.default_rng(SEED + 1).uniform(-.1, .1, len(ra)), 100 * ra.treg_share_of_T, label="Treg (T-8, T-9)")
    ax.set_xlabel("lymphocytic inflammation grade"); ax.set_ylabel("% of T cells"); ax.legend(fontsize=8)
    ax.set_title("P1 GSE243917 donors", fontsize=10)
    ax = axes[1]
    order = ["healthy", "arthralgia", "UA", "early RA", "established RA", "OA"]
    for i, gname in enumerate(order):
        v = bulk.loc[bulk.group == gname, "Tph"]
        ax.boxplot(v, positions=[i], widths=0.5, showfliers=False)
    ax.set_xticks(range(len(order)), order, rotation=30, fontsize=8); ax.set_ylabel("Tph programme (mean z)")
    t = res.set_index("test")
    ax.set_title(f"P2 GSE89408 Tph: RA vs OA beyond random p {t.loc['P2 Tph programme RA vs OA (adjusted incl. T-cell abundance)', 'p_vs_random']:.2g}", fontsize=10)
    ax = axes[2]
    dra = bulk[bulk.group.isin(["early RA", "established RA"])]
    ax.scatter(dra.Tph, dra.IFNp, s=10)
    ax.set_xlabel("Tph programme"); ax.set_ylabel("IFN macrophage programme"); ax.set_title("P2 within RA biopsies", fontsize=10)
    ax = axes[3]
    if len(sp):
        x = np.arange(len(sp))
        ax.bar(x - 0.2, sp.B_neighbours_CXCL13pos_T, width=0.4, label="around CXCL13+ T")
        ax.bar(x + 0.2, sp.B_neighbours_CXCL13neg_T, width=0.4, label="around CXCL13- T")
        ax.set_xticks(x, sp.block, rotation=20, fontsize=8); ax.set_ylabel("B cells within 30 um"); ax.legend(fontsize=8)
    ax.set_title("P4 Xenium 5K RA sections", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "tph_treg_axis.png", dpi=150)


if __name__ == "__main__":
    main()
