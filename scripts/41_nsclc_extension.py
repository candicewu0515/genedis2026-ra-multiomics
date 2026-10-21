"""Cross-indication extension: do the RA synovial M-6 (IFN-gamma macrophage) and Tph programmes behave alike in NSCLC?

Programmes: M6_programme and Tph_programme, 30 genes each, defined in RA synovium (scripts/state_markers.py) and not
re-fitted here. Data: LuCA core atlas (Salcher et al. 2022 Cancer Cell; CELLxGENE, 298 donors; lung tissue only) and
TCGA LUAD / LUSC bulk RNA-seq (UCSC Xena, as in step 19).

Hypotheses fixed before any NSCLC data were scored (two-sided):
  N1  an M-6-like state exists among primary-tumour macrophages: after clustering non-alveolar tumour macrophages, the
      cluster with the highest M6 score stands out more than the best cluster of 200 expression-matched random gene sets;
  N2  primary: donor-level macrophage pseudobulk M6 score differs between primary tumour and adjacent normal lung
      (paired Wilcoxon over donors with >= 50 macrophages in both; also compared with the same paired difference for
      200 matched random gene sets);
  N3  coupling, as in RA synovium: across tumours, macrophage M6 score correlates with the Tph score of CD4 T cells
      (Spearman on study-adjusted residuals; rho compared with 200 random macrophage gene sets);
  N4  TCGA: M6 and Tph scores of primary tumours vs overall survival (Cox per SD, adjusted for age and stage I-II vs
      III-IV; LUAD and LUSC separately); a second model also adjusts for leukocyte abundance (PTPRC), so that an
      association is not just more immune cells.
Units: non-alveolar macrophages ('Macrophage' in LuCA cell_type_major) and CD4 T cells ('T cell CD4'); a unit is one
donor x origin with >= 50 cells. Scores: per-gene z of log2(CPM + 1) across units, averaged over the programme.

Output: results/tables/nsclc_extension_tests.tsv, nsclc_extension_cox.tsv, luca_m6_cluster_markers.tsv;
        data/processed/luca_unit_scores.tsv; results/figures/nsclc_extension.png
"""
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from anndata.io import sparse_dataset
from scipy.stats import wilcoxon, spearmanr
from statsmodels.duration.hazard_regression import PHReg

from config import RAW, PROCESSED, FIG, TAB, SEED
from state_markers import markers, matched_sets, unit_scores, residualise, logcpm, cluster_best_state

LUCA = RAW / "luca" / "luca_core_atlas.h5ad"
T = RAW / "tcga"
N_RANDOM = 200
MIN_CELLS = 50


# ---------- LuCA ----------
def luca_cells():
    with h5py.File(LUCA, "r") as f:
        obs = {k: f[f"obs/{k}"] for k in ("donor_id", "origin", "tissue", "cell_type_major", "study", "disease")}
        def col(name):
            g = obs[name]
            if isinstance(g, h5py.Group):                       # categorical
                cats = g["categories"][:].astype(str); codes = g["codes"][:]
                return np.where(codes >= 0, cats[np.maximum(codes, 0)], "nan")
            return g[:].astype(str)
        o = pd.DataFrame({k: col(k) for k in obs})
        sym = f["raw/var/feature_name"]
        sym = (sym["categories"][:].astype(str)[sym["codes"][:]] if isinstance(sym, h5py.Group) else sym[:].astype(str))
    keep = (o.tissue == "lung") & o.origin.isin(["tumor_primary", "normal_adjacent", "normal"]) & \
           o.cell_type_major.isin(["Macrophage", "T cell CD4"])
    return o, keep.to_numpy(), pd.Index(sym)


def load_rows(idx):
    with h5py.File(LUCA, "r") as f:
        X = sparse_dataset(f["raw/X"])
        return X[np.sort(idx)].tocsr()


def pseudobulk(X, obs, genes):
    key = obs["donor_id"] + "|" + obs["origin"]
    units = key.value_counts()
    units = units[units >= MIN_CELLS].index
    rows = []
    for u in units:
        m = (key == u).to_numpy()
        rows.append(np.asarray(X[m].sum(axis=0)).ravel())
    pb = pd.DataFrame(np.vstack(rows), index=units, columns=genes)
    pb = pb.T.groupby(level=0).sum().T                          # duplicated symbols summed
    meta = obs.assign(unit=key).drop_duplicates("unit").set_index("unit").loc[units, ["donor_id", "origin", "study", "disease"]]
    meta["n_cells"] = key.value_counts().loc[units].values
    return pb, meta


def luca(m6, tph, rng):
    o, keep, sym = luca_cells()
    idx = np.where(keep)[0]
    print(f"LuCA: loading {len(idx):,} lung macrophages and CD4 T cells")
    X = load_rows(idx)
    obs = o.iloc[np.sort(idx)].reset_index(drop=True)
    rows = []

    mac = (obs.cell_type_major == "Macrophage").to_numpy()
    mobs = obs[mac].reset_index(drop=True)
    sel = np.where(mobs.origin.to_numpy() == "tumor_primary")[0]
    if len(sel) > 40000:
        sel = np.sort(rng.choice(sel, 40000, replace=False))
    row1, mk, cellres = cluster_best_state(X[mac][sel], mobs.iloc[sel].reset_index(drop=True), sym, m6, "study",
                                           "donor_id", rng, N_RANDOM)
    row1["test"] = "N1 M-6-like cluster among primary-tumour macrophages"
    rows.append(row1); mk.to_csv(TAB / "luca_m6_cluster_markers.tsv", sep="\t", index=False, float_format="%.4g")
    print(row1)

    pbm, metam = pseudobulk(X[mac], obs[mac].reset_index(drop=True), sym)
    pbt, metat = pseudobulk(X[~mac], obs[~mac].reset_index(drop=True), sym)
    Lm, Lt = logcpm(pbm), logcpm(pbt)
    expressed = Lm.columns[(Lm > 1).mean() >= 0.5]
    Lm = Lm[expressed]
    Lt = Lt[Lt.columns[(Lt > 1).mean() >= 0.5]]
    mean_m = Lm.mean()
    rand_m = matched_sets([g for g in m6 if g in Lm.columns], mean_m, N_RANDOM, rng)
    sets = {"M6": m6, **{f"r{i}": r for i, r in enumerate(rand_m)}}
    Sm = unit_scores(Lm, sets)
    St = unit_scores(Lt, {"Tph": tph})
    units = metam.join(Sm[["M6"]]).rename(columns={"M6": "macrophage_M6"})
    units = units.join(metat[[]].join(St).rename(columns={"Tph": "cd4_Tph"}), how="left")
    units.round(4).to_csv(PROCESSED / "luca_unit_scores.tsv", sep="\t")

    # N2 tumour vs adjacent normal, paired
    piv = Sm.join(metam[["donor_id", "origin"]])
    t = piv[piv.origin == "tumor_primary"].set_index("donor_id")
    nn = piv[piv.origin == "normal_adjacent"].set_index("donor_id")
    both = t.index.intersection(nn.index)
    d_obs = (t.loc[both, "M6"] - nn.loc[both, "M6"])
    d_null = np.array([np.median(t.loc[both, f"r{i}"] - nn.loc[both, f"r{i}"]) for i in range(N_RANDOM)])
    rows.append(dict(test="N2 macrophage M6 score, tumour vs adjacent normal (paired)", n=int(len(both)),
                     statistic=float(np.median(d_obs)), null_median=float(np.median(d_null)),
                     p=float(wilcoxon(d_obs).pvalue),
                     p_vs_random=(1 + (np.abs(d_null - np.median(d_null)) >= abs(np.median(d_obs) - np.median(d_null))).sum()) / (1 + N_RANDOM),
                     detail=f"tumour higher in {(d_obs > 0).sum()}/{len(both)} donors"))

    # N3 coupling within tumours
    tt = units[units.origin == "tumor_primary"].dropna(subset=["cd4_Tph"])
    rx = residualise(tt["macrophage_M6"], tt["study"]); ry = residualise(tt["cd4_Tph"], tt["study"])
    rho = spearmanr(rx, ry)
    tm = Sm.loc[tt.index]
    null_rho = np.array([spearmanr(residualise(tm[f"r{i}"], tt["study"]), ry)[0] for i in range(N_RANDOM)])
    rows.append(dict(test="N3 tumour macrophage M6 vs CD4 T Tph score (study-adjusted)", n=int(len(tt)),
                     statistic=float(rho[0]), null_median=float(np.median(null_rho)), p=float(rho[1]),
                     p_vs_random=(1 + (null_rho >= rho[0]).sum()) / (1 + N_RANDOM),
                     detail=f"{tt.study.nunique()} studies"))
    return rows, units, (d_obs, d_null, tt, rx, ry, cellres)


# ---------- TCGA ----------
def stage_group(s):
    s = str(s).upper()
    return 1 if ("IV" in s or "III" in s) else (0 if "I" in s else np.nan)


def tcga(m6, tph):
    surv = pd.read_csv(T / "Survival_SupplementalTable_S1_20171025_xena_sp.tsv", sep="\t", index_col=0)
    rows, tn = [], []
    for cancer in ("LUAD", "LUSC"):
        x = pd.read_csv(T / f"{cancer}_HiSeqV2.gz", sep="\t", index_col=0)
        clin = pd.read_csv(T / f"{cancer}_clinicalMatrix.tsv", sep="\t", index_col=0, low_memory=False)
        stype = x.columns.str[13:15]
        tum = x.columns[stype == "01"]
        z = x.sub(x[tum].mean(axis=1), axis=0).div(x[tum].std(axis=1).replace(0, np.nan), axis=0)
        s = pd.DataFrame({"M6": z.loc[[g for g in m6 if g in z.index]].mean(),
                          "Tph": z.loc[[g for g in tph if g in z.index]].mean(), "PTPRC": z.loc["PTPRC"]})
        norm = x.columns[stype == "11"]
        pairs = [(a, b) for b in norm for a in tum if a[:12] == b[:12]]
        for k in ("M6", "Tph"):
            d = np.array([s.loc[a, k] - s.loc[b, k] for a, b in pairs])
            tn.append(dict(test=f"N4 TCGA {cancer} {k}: tumour vs adjacent normal (paired)", n=len(pairs),
                           statistic=float(np.median(d)), p=float(wilcoxon(d).pvalue), detail=f"tumour higher in {(d > 0).sum()}/{len(d)}"))
        t = s.loc[tum].copy(); t["patient"] = t.index.str[:15]
        t = t.join(surv[["OS", "OS.time"]], on="patient")
        t = t.join(clin[["age_at_initial_pathologic_diagnosis", "pathologic_stage"]], how="left")
        t["age"] = pd.to_numeric(t["age_at_initial_pathologic_diagnosis"], errors="coerce")
        t["stage_III_IV"] = t["pathologic_stage"].map(stage_group)
        t = t.dropna(subset=["OS", "OS.time", "age", "stage_III_IV"])
        t = t[(t["OS.time"] > 0) & ~t.index.str[:12].duplicated()]
        for k in ("M6", "Tph"):
            for adj in (False, True):
                X = pd.DataFrame({"score_per_SD": (t[k] - t[k].mean()) / t[k].std(), "age_per_10y": t["age"] / 10,
                                  "stage_III_IV": t["stage_III_IV"]})
                if adj:
                    X["PTPRC_per_SD"] = (t["PTPRC"] - t["PTPRC"].mean()) / t["PTPRC"].std()
                r = PHReg(t["OS.time"].values, X.values, status=t["OS"].values.astype(int)).fit()
                hr = np.exp(r.params); ci = np.exp(r.conf_int())
                rows.append(dict(cancer=cancer, score=k, adjusted_for_PTPRC=adj, n=len(t), events=int(t["OS"].sum()),
                                 HR_per_SD=hr[0], CI_low=ci[0, 0], CI_high=ci[0, 1], p=r.pvalues[0]))
    return pd.DataFrame(rows), tn


def main():
    rng = np.random.default_rng(SEED)
    mk = markers()
    m6 = mk.loc[mk.programme == "M6_programme", "gene"].tolist()
    tph = mk.loc[mk.programme == "Tph_programme", "gene"].tolist()
    rows, units, aux = luca(m6, tph, rng)
    cox, tn = tcga(m6, tph)
    tests = pd.DataFrame(rows + tn)
    tests.to_csv(TAB / "nsclc_extension_tests.tsv", sep="\t", index=False, float_format="%.4g")
    cox.to_csv(TAB / "nsclc_extension_cox.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 90)
    print(tests.round(4).to_string(index=False)); print(cox.round(3).to_string(index=False))

    d_obs, d_null, tt, rx, ry, cellres = aux
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    ax[0].hist(d_null, bins=30, color="0.75"); ax[0].axvline(np.median(d_obs), color="C3")
    ax[0].set_title("N2 tumour - adjacent normal (median), M6 vs random sets", fontsize=9)
    ax[1].scatter(rx, ry, s=12); ax[1].set_xlabel("macrophage M6 (study-adjusted)"); ax[1].set_ylabel("CD4 T Tph (study-adjusted)")
    ax[1].set_title("N3 coupling across tumours", fontsize=9)
    c = cox[~cox.adjusted_for_PTPRC]
    ax[2].errorbar(c.HR_per_SD, range(len(c)), xerr=[c.HR_per_SD - c.CI_low, c.CI_high - c.HR_per_SD], fmt="o")
    ax[2].set_yticks(range(len(c)), c.cancer + " " + c.score); ax[2].axvline(1, ls="--", color="0.5")
    ax[2].set_title("N4 TCGA overall survival (HR per SD)", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "nsclc_extension.png", dpi=150)


if __name__ == "__main__":
    main()
