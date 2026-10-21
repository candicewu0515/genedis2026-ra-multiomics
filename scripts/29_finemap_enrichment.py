"""Fine-mapped RA risk variants in cell-state-specific chromatin (upgrade of the locus-window test in step 10).

Credible sets (Open Targets Platform 26.09 GraphQL, cached):
  primary      GCST90132223  Ishigaki 2022 RA, European ancestry (SuSiE-inf, harmonised GWAS Catalog sumstats)
  replication  FINNGEN_R12_M13_RHEUMA  FinnGen R12 RA (SuSiE). FinnGen is not among the Ishigaki 2022 cohorts
               (no mention in the paper's full text), so the two are treated as independent.
  ancestry     GCST90132224  Ishigaki 2022 RA, East Asian ancestry (SuSiE-inf; separate cohorts from the European
               analysis); tests whether the same states and loci emerge under a different LD structure.
  control      GCST90018961  LDL cholesterol (Sakaue 2021; 83 % European, 17 % East Asian), same pipeline
Variants in 95 % credible sets; MHC (chr6:28.51-33.48 Mb) excluded.
Peak sets: the same top-4,000 state-specific peaks as step 10 (cell types vs all; subtypes vs same lineage).
Statistic: posterior inclusion probability (PIP) mass of credible-set variants within +/-250 bp of a state-specific
peak midpoint (each variant counted once). Null: 1,000 random peak sets matched on accessibility decile x promoter
status (as step 10). Enrichment = observed / null mean; empirical p; BH within level x study.
Hypotheses fixed before analysis:
  F1 T-cell peaks are enriched for Ishigaki PIP mass (step 10: 1.59x, FDR 0.008)
  F2 the T-cell enrichment replicates in FinnGen (one-sided p < 0.05), and state-level enrichments correlate between
     the two RA studies (Spearman across states)
  F3 LDL shows no enrichment in synovial immune states
  F4 (added in a later round) the T-cell enrichment is also present in East Asian ancestry (one-sided p < 0.05), and
     state-level enrichments correlate with the European analysis
Candidate causal variants: PIP >= 0.1 inside a state-specific peak, with nearest and linked genes.

Output: data/raw/opentargets/finemap_<study>.json, results/tables/finemap_variants.tsv, finemap_enrichment.tsv,
        finemap_candidates.tsv; results/figures/finemap_enrichment.png
"""
import importlib
import json
import time
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

from config import RAW, INTERIM, FIG, TAB, SEED

step10 = importlib.import_module("10_gwas_peak_enrichment")
OT = "https://api.platform.opentargets.org/api/v4/graphql"
STUDIES = {"GCST90132223": "RA Ishigaki 2022 (EUR)", "FINNGEN_R12_M13_RHEUMA": "RA FinnGen R12",
           "GCST90132224": "RA Ishigaki 2022 (EAS)",
           "GCST90018961": "LDL Sakaue 2021 (control)"}
HALF = 250
N_PERM = 1000
MHC = step10.MHC
C = RAW / "opentargets"
C.mkdir(exist_ok=True)

Q = """query($id: String!, $i: Int!) { study(studyId: $id) { credibleSets(page: {index: $i, size: 10}) { count rows {
  studyLocusId finemappingMethod purityMinR2 variant { id rsIds }
  locus(page: {index: 0, size: 1000}) { count rows { variant { id rsIds } posteriorProbability is95CredibleSet } } } } } }"""


def gql(query, variables):
    req = urllib.request.Request(OT, data=json.dumps({"query": query, "variables": variables}).encode(),
                                 headers={"Content-Type": "application/json"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                d = json.load(r)
            if "errors" in d:
                raise RuntimeError(d["errors"])
            return d["data"]
        except Exception:
            if attempt == 4:
                raise
            time.sleep(10 * (attempt + 1))


def credible_sets(study):
    p = C / f"finemap_{study}.json"
    if p.exists():
        return json.loads(p.read_text())
    rows, i = [], 0
    while True:
        cs = gql(Q, {"id": study, "i": i})["study"]["credibleSets"]
        rows += cs["rows"]
        if len(rows) >= cs["count"] or not cs["rows"]:
            break
        i += 1
    for r in rows:
        if r["locus"]["count"] > len(r["locus"]["rows"]):
            raise RuntimeError(f"{study} {r['studyLocusId']}: locus has more than 1,000 variants; paginate")
    p.write_text(json.dumps(rows))
    return rows


def variants_table(study, rows):
    out = []
    for r in rows:
        lead = r["variant"]["id"]
        for v in r["locus"]["rows"]:
            if not v["is95CredibleSet"]:
                continue
            chrom, pos, ref, alt = v["variant"]["id"].split("_")
            out.append(dict(study=study, studyLocusId=r["studyLocusId"], method=r["finemappingMethod"], lead=lead,
                            variant=v["variant"]["id"], rsid=";".join(v["variant"]["rsIds"] or []), chrom="chr" + chrom,
                            pos=int(pos), pip=float(v["posteriorProbability"])))
    d = pd.DataFrame(out)
    mhc = (d.chrom == MHC[0]) & d.pos.between(MHC[1], MHC[2])
    return d[~mhc]


def main():
    rng = np.random.default_rng(SEED)
    atac = sc.read_h5ad(INTERIM / "gse243917_atac.h5ad")
    peaks = atac.var[["chrom", "start", "end", "frac_cells"]].copy()
    ok_peak = ~((peaks["chrom"] == MHC[0]) & (peaks["end"] >= MHC[1]) & (peaks["start"] <= MHC[2])).values
    peaks = peaks[ok_peak]
    mid = ((peaks["start"] + peaks["end"]) // 2).values
    tss = step10.tss_table()
    near_tss = np.zeros(len(peaks), bool)
    nearest = np.empty(len(peaks), dtype=object)
    for chrom, g in tss.groupby("chrom"):
        idx = np.where(peaks["chrom"].values == chrom)[0]
        if not len(idx):
            continue
        g = g.sort_values("tss"); t, names = g["tss"].values, g["gene"].values
        j = np.clip(np.searchsorted(t, mid[idx]), 1, len(t) - 1)
        jj = np.where(np.abs(mid[idx] - t[j - 1]) <= np.abs(mid[idx] - t[j]), j - 1, j)
        near_tss[idx] = np.abs(mid[idx] - t[jj]) <= 2000
        nearest[idx] = names[jj]
    strata = pd.qcut(peaks["frac_cells"].rank(method="first"), 10, labels=False).values * 2 + near_tss.astype(int)
    pools = {s: np.where(strata == s)[0] for s in np.unique(strata)}
    sets = step10.specific_sets(atac, ok_peak)
    print(f"Peak sets: {len(sets)}; universe {len(peaks):,} non-MHC peaks")

    var = pd.concat([variants_table(s, credible_sets(s)) for s in STUDIES], ignore_index=True)
    var.to_csv(TAB / "finemap_variants.tsv", sep="\t", index=False)
    for s, g in var.groupby("study"):
        print(f"{STUDIES[s]}: {g.studyLocusId.nunique()} non-MHC credible sets, {len(g)} variants, "
              f"{(g.pip >= 0.5).sum()} with PIP >= 0.5")

    links = pd.read_csv(TAB / "peak_gene_links.tsv.gz", sep="\t")
    links = links[(links["r"] > 0.45) & (links["fdr"] < 0.05)]
    peak_names = peaks.index.values
    rows, cands = [], []
    for study, v in var.groupby("study"):
        v = v.reset_index(drop=True)
        # sparse peak x variant overlap (variant within +/-250 bp of the peak midpoint)
        pi, vi = [], []
        for chrom, g in v.groupby("chrom"):
            idx = np.where(peaks["chrom"].values == chrom)[0]
            if not len(idx):
                continue
            o = np.argsort(mid[idx]); m = mid[idx][o]
            lo = np.searchsorted(m, g.pos.values - HALF, side="left")
            hi = np.searchsorted(m, g.pos.values + HALF, side="right")
            for k, a, b in zip(g.index, lo, hi):
                for j in range(a, b):
                    pi.append(idx[o[j]]); vi.append(k)
        M = sp.csr_matrix((np.ones(len(pi)), (pi, vi)), shape=(len(peaks), len(v)))
        pip = v.pip.values
        total = pip.sum()

        def mass(sel_idx):
            sel = np.zeros(len(peaks)); sel[sel_idx] = 1
            return float(pip @ ((M.T @ sel) > 0))
        for (level, state), idx in sets.items():
            obs = mass(idx)
            sel = np.zeros(len(peaks)); sel[idx] = 1
            cov = np.where((M.T @ sel) > 0)[0]
            per_set = pd.Series(pip[cov]).groupby(v.studyLocusId.values[cov]).sum()
            need = pd.Series(strata[idx]).value_counts()
            null = np.array([mass(np.concatenate([rng.choice(pools[s], n, replace=False) for s, n in need.items()]))
                             for _ in range(N_PERM)])
            rows.append(dict(study=study, label=STUDIES[study], level=level, state=state, pip_mass=obs,
                             pct_of_total_pip=100 * obs / total, null_mean=null.mean(),
                             enrichment=obs / max(null.mean(), 1e-9), p_empirical=(1 + np.sum(null >= obs)) / (1 + N_PERM),
                             n_sets_contributing=int((per_set > 0.01).sum()),
                             top_set_share=float(per_set.max() / obs) if obs > 0 else np.nan,
                             top_set_lead=v.lead[v.studyLocusId == per_set.idxmax()].iloc[0] if obs > 0 else ""))
            if study != "GCST90018961":
                hit = np.unique(M[idx].nonzero()[1])
                for k in hit[pip[hit] >= 0.1]:
                    pk = [peak_names[j] for j in idx if M[j, k]]
                    for p_ in pk:
                        j = peaks.index.get_loc(p_)
                        cands.append(dict(study=STUDIES[study], level=level, state=state, variant=v.variant[k],
                                          rsid=v.rsid[k], pip=pip[k], credible_set=v.studyLocusId[k], peak=p_,
                                          nearest_protein_coding_gene=nearest[j],
                                          linked_genes=";".join(links.loc[links.peak == p_, "gene"])))
        print(f"  {STUDIES[study]} done", flush=True)
    res = pd.DataFrame(rows)
    res["fdr"] = np.nan
    for _, g in res.groupby(["study", "level"]):
        res.loc[g.index, "fdr"] = multipletests(g.p_empirical, method="fdr_bh")[1]
    res = res.sort_values(["study", "level", "p_empirical"])
    res.round(4).to_csv(TAB / "finemap_enrichment.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    for s in STUDIES:
        print(f"\n== {STUDIES[s]}")
        print(res[res.study == s].groupby("level").head(6)[["level", "state", "pip_mass", "null_mean", "enrichment", "p_empirical",
                                                             "fdr", "n_sets_contributing", "top_set_share", "top_set_lead"]].round(3).to_string(index=False))
    ish = res[res.study == "GCST90132223"].set_index(["level", "state"])
    fin = res[res.study == "FINNGEN_R12_M13_RHEUMA"].set_index(["level", "state"])
    j = ish.join(fin, lsuffix="_ish", rsuffix="_fin")
    rho, prho = spearmanr(np.log2(j.enrichment_ish + 0.05), np.log2(j.enrichment_fin + 0.05))
    t = j.loc[("cell type", "Tcell")]
    print(f"\nF1 T cells, Ishigaki: {t.enrichment_ish:.2f}x, p {t.p_empirical_ish:.3f}, FDR {t.fdr_ish:.3f}; "
          f"{t.n_sets_contributing_ish} credible sets contribute, largest share {t.top_set_share_ish:.2f} ({t.top_set_lead_ish})")
    print(f"F2 T cells, FinnGen: {t.enrichment_fin:.2f}x, p {t.p_empirical_fin:.3f}, FDR {t.fdr_fin:.3f}; "
          f"{t.n_sets_contributing_fin} credible sets, largest share {t.top_set_share_fin:.2f} ({t.top_set_lead_fin}); "
          f"state-level Spearman Ishigaki vs FinnGen rho {rho:.2f} (p {prho:.2g}, {len(j)} states)")
    eas = res[res.study == "GCST90132224"].set_index(["level", "state"])
    je = ish.join(eas, lsuffix="_ish", rsuffix="_eas")
    rho_e, p_e = spearmanr(np.log2(je.enrichment_ish + 0.05), np.log2(je.enrichment_eas + 0.05))
    te = eas.loc[("cell type", "Tcell")]
    print(f"F4 T cells, East Asian: {te.enrichment:.2f}x, p {te.p_empirical:.3f}, FDR {te.fdr:.3f}; {te.n_sets_contributing} "
          f"credible sets, largest share {te.top_set_share:.2f} ({te.top_set_lead}); state-level Spearman EUR vs EAS "
          f"rho {rho_e:.2f} (p {p_e:.2g})")
    ldl = res[res.study == "GCST90018961"]
    print(f"F3 LDL: states with FDR < 0.05: {int((ldl.fdr < 0.05).sum())} of {len(ldl)}; max enrichment {ldl.enrichment.max():.2f}")
    pd.DataFrame([dict(test="F1 T-cell PIP enrichment, Ishigaki", value=t.enrichment_ish, p=t.p_empirical_ish, fdr=t.fdr_ish),
                  dict(test="F2 T-cell PIP enrichment, FinnGen replication", value=t.enrichment_fin, p=t.p_empirical_fin, fdr=t.fdr_fin),
                  dict(test="F2 state-level concordance Ishigaki vs FinnGen (Spearman rho)", value=rho, p=prho, fdr=np.nan),
                  dict(test="F3 LDL states with FDR < 0.05", value=int((ldl.fdr < 0.05).sum()), p=np.nan, fdr=np.nan),
                  dict(test="F4 T-cell PIP enrichment, East Asian ancestry", value=te.enrichment, p=te.p_empirical, fdr=te.fdr),
                  dict(test="F4 state-level concordance EUR vs EAS (Spearman rho)", value=rho_e, p=p_e, fdr=np.nan)]) \
        .to_csv(TAB / "finemap_tests.tsv", sep="\t", index=False, float_format="%.4g")

    cand = pd.DataFrame(cands)
    if len(cand):
        cand = cand.sort_values(["pip"], ascending=False)
        old = set(pd.read_csv(TAB / "gwas_candidate_variants.tsv", sep="\t").rs_id)
        cand["in_step10_candidates"] = cand.rsid.apply(lambda s: any(x in old for x in str(s).split(";")))
        cand.to_csv(TAB / "finemap_candidates.tsv", sep="\t", index=False)
        top = cand.drop_duplicates(["variant"]).head(20)
        print(f"\nFine-mapped variants (PIP >= 0.1) in state-specific peaks: {cand.variant.nunique()} variants "
              f"({cand[cand.in_step10_candidates].variant.nunique()} also step-10 lead SNPs)")
        print(cand.groupby("variant").agg(rsid=("rsid", "first"), pip=("pip", "max"), study=("study", lambda s: ",".join(sorted(set(s)))),
                                           states=("state", lambda s: ",".join(sorted(set(s)))[:70]),
                                           gene=("nearest_protein_coding_gene", "first"), linked=("linked_genes", "first"))
              .sort_values("pip", ascending=False).head(25).to_string())

    fig, axes = plt.subplots(1, 2, figsize=(16, 9), gridspec_kw={"width_ratios": [1, 2.4]})
    for ax, level in zip(axes, ("cell type", "subtype")):
        p = res[res.level == level].pivot_table(index="state", columns="label", values="enrichment")
        f = res[res.level == level].pivot_table(index="state", columns="label", values="fdr")
        p = p.sort_values(STUDIES["GCST90132223"])
        y = np.arange(len(p))
        hgt = 0.8 / len(STUDIES)
        for k, (lab, col) in enumerate(zip(STUDIES.values(), ("C3", "C1", "C2", "C7"))):
            off = hgt * (len(STUDIES) / 2 - 0.5 - k)
            ax.barh(y + off, p[lab], height=hgt, color=col, label=lab)
            for i, st in enumerate(p.index):
                if f.loc[st, lab] < 0.05:
                    ax.text(p.loc[st, lab] + 0.03, i + off, "*", va="center", fontsize=9)
        ax.axvline(1, color="k", lw=0.8, ls="--")
        ax.set_yticks(y, p.index, fontsize=7 if level == "subtype" else 9)
        ax.set_xlabel("PIP mass in specific peaks / matched-null mean")
        ax.set_title(f"{level}: fine-mapped variants (* FDR < 0.05)", fontsize=10)
        ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "finemap_enrichment.png", dpi=150)


if __name__ == "__main__":
    main()
