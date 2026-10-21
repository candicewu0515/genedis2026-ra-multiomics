"""Exploratory expression-associated Wald ratios at colocalising RA loci: odds of RA per SD change in a gene's expression or
protein level in immune cells.

Instruments: for each RA signal x gene pair that colocalises (H4 >= 0.8, cis) with an immune-cell or blood level QTL
(eQTL, single-cell eQTL, cis-pQTL) in step 32, the lead variant of the QTL credible set. QTL effect (beta, SE, per SD
of normalised expression / protein) from Open Targets 26.09; GWAS effect for the same variant and allele from the full
summary statistics (Ishigaki 2022 EUR or FinnGen R12, matching the RA signal). Wald ratio b_GWAS / b_QTL with
first-order SE; reported as OR per SD increase.
Datasets measuring the same variant are not independent instruments, so per signal x gene the primary estimate is the
dataset with the highest CLPP (immune cells preferred over blood/plasma), and the share of datasets with the same sign
is reported alongside. Colocalisation supports a shared variant; a single-instrument ratio cannot exclude pleiotropy.
Hypotheses fixed before analysis:
  M1 CTLA4: OR per SD < 1. Abatacept binds CD80/CD86; this is not a genetic proxy for the drug.
  M2 IL6ST (ANKRD55-shared instrument): OR per SD > 1. IL6ST is not the IL6R target of tocilizumab.
  M3 IRF5, CD40, TYK2, CCR6, FCGR2A: report estimates and cross-dataset consistency (no direction fixed)

Output: data/raw/opentargets/qtl_lead_<studyLocusId>.json (cached), results/tables/mr_wald_all.tsv, mr_primary.tsv;
        results/figures/drug_target_mr.png
"""
import json
import time
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import RAW, FIG, TAB

OT = "https://api.platform.opentargets.org/api/v4/graphql"
C = RAW / "opentargets"
LEVEL = {"eqtl", "sceqtl", "pqtl", "scpqtl"}
CAT_ORDER = {"myeloid": 0, "T/NK": 0, "B/LCL": 0, "granulocyte": 0, "blood/plasma": 1}
SUMSTATS = {"Ishigaki 2022 (EUR)": RAW / "gwas_sumstats" / "GCST90132223.h.tsv.gz",
            "FinnGen R12": RAW / "gwas_sumstats" / "finngen_R12_M13_RHEUMA.gz"}
FOCUS = ["CTLA4", "IL6ST", "IRF5", "CD40", "TYK2", "CCR6", "FCGR2A", "ANKRD55", "TRAF1", "IL2RA", "PTPN22"]


def gql(query):
    req = urllib.request.Request(OT, data=json.dumps({"query": query}).encode(), headers={"Content-Type": "application/json"})
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


def qtl_leads(ids):
    """Lead variant, beta and SE of QTL credible sets (cached per id; 10 per request)."""
    out, todo = {}, []
    for i in ids:
        p = C / f"qtl_lead_{i}.json"
        if p.exists():
            out[i] = json.loads(p.read_text())
        else:
            todo.append(i)
    for k in range(0, len(todo), 10):
        chunk = todo[k:k + 10]
        q = "{ " + " ".join(f'a{j}: credibleSet(studyLocusId: "{sl}") {{ studyLocusId variant {{ id }} beta standardError }}'
                            for j, sl in enumerate(chunk)) + " }"
        d = gql(q)
        for j, sl in enumerate(chunk):
            rec = d.get(f"a{j}") or {}
            (C / f"qtl_lead_{sl}.json").write_text(json.dumps(rec))
            out[sl] = rec
    return out


def gwas_effect(path, fmt_finngen, variants):
    """Beta and SE of the alt allele of 'chr_pos_ref_alt' variants from a bgzipped/tabix-less file: scan by position."""
    want = {}
    for v in variants:
        c, p, ref, alt = v.split("_")
        want.setdefault(c, {})[int(p)] = (ref, alt, v)
    res = {}
    for d in pd.read_csv(path, sep="\t", chunksize=2_000_000, dtype=str):
        if fmt_finngen:
            d = d.rename(columns={"#chrom": "c", "pos": "p", "ref": "a2", "alt": "a1", "beta": "b", "sebeta": "se"})
        else:
            d = d.rename(columns={"chromosome": "c", "base_pair_location": "p", "other_allele": "a2", "effect_allele": "a1",
                                  "beta": "b", "standard_error": "se"})
        d = d[d.c.isin(want.keys())]
        if d.empty:
            continue
        d = d.assign(p=d.p.astype(int))
        for c, g in d.groupby("c"):
            sub = g[g.p.isin(want[c].keys())]
            for r in sub.itertuples():
                ref, alt, v = want[c][r.p]
                b, se = float(r.b), float(r.se)
                if (r.a1, r.a2) == (alt, ref):
                    res[v] = (b, se)
                elif (r.a1, r.a2) == (ref, alt):
                    res[v] = (-b, se)
    return res


def main():
    allc = pd.read_csv(TAB / "coloc_multisignal_all.tsv", sep="\t")
    hit = allc[(allc.h4 >= 0.8) & allc.cis & allc.qtl_type.isin(LEVEL) & allc.category.isin(CAT_ORDER.keys())].copy()
    # recover the QTL credible-set ids from the cached step-32 responses
    rows = []
    for sl in hit.ra_set.unique():
        for r in json.loads((C / f"coloc_{sl}.json").read_text()):
            o = r["otherStudyLocus"]
            sig4 = lambda x: float(f"{x:.4g}") if x is not None else np.nan     # step 32 wrote %.4g
            rows.append(dict(ra_set=sl, qtl_study=o["study"]["id"], gene=(o["study"].get("target") or {}).get("approvedSymbol"),
                             qtl_locus=o["studyLocusId"], h4=sig4(r["h4"]), clpp=sig4(r["clpp"])))
    ids = pd.DataFrame(rows).drop_duplicates(["ra_set", "qtl_study", "gene", "h4", "clpp"])
    hit = hit.merge(ids, on=["ra_set", "qtl_study", "gene", "h4", "clpp"], how="left").dropna(subset=["qtl_locus"])
    leads = qtl_leads(sorted(hit.qtl_locus.unique()))
    hit["qtl_variant"] = hit.qtl_locus.map(lambda i: ((leads.get(i) or {}).get("variant") or {}).get("id"))
    hit["b_qtl"] = hit.qtl_locus.map(lambda i: (leads.get(i) or {}).get("beta"))
    hit["se_qtl"] = hit.qtl_locus.map(lambda i: (leads.get(i) or {}).get("standardError"))
    hit = hit.dropna(subset=["qtl_variant", "b_qtl"])
    print(f"Level-QTL colocalisations with a QTL lead variant: {len(hit)} ({hit.gene.nunique()} genes)")

    parts = []
    for gw, path in SUMSTATS.items():
        sub = hit[hit.gwas == gw].copy()
        eff = gwas_effect(path, "FinnGen" in gw, sub.qtl_variant.unique())
        sub["b_gwas"] = sub.qtl_variant.map(lambda v: eff.get(v, (np.nan, np.nan))[0])
        sub["se_gwas"] = sub.qtl_variant.map(lambda v: eff.get(v, (np.nan, np.nan))[1])
        parts.append(sub)
        print(f"  {gw}: GWAS effects found for {sub.b_gwas.notna().sum()} of {len(sub)} instruments")
    mr = pd.concat(parts).dropna(subset=["b_gwas"])
    mr = mr[mr.b_qtl.abs() > 0]
    mr["wald"] = mr.b_gwas / mr.b_qtl
    mr["wald_se"] = mr.se_gwas / mr.b_qtl.abs()
    mr["OR_per_SD"] = np.exp(mr.wald)
    mr["OR_lo"] = np.exp(mr.wald - 1.96 * mr.wald_se); mr["OR_hi"] = np.exp(mr.wald + 1.96 * mr.wald_se)
    mr["z"] = mr.wald / mr.wald_se
    mr.to_csv(TAB / "mr_wald_all.tsv", sep="\t", index=False, float_format="%.4g")

    prim = []
    for (gw, sl, gene), g in mr.groupby(["gwas", "ra_set", "gene"]):
        g = g.assign(o=g.category.map(CAT_ORDER)).sort_values(["o", "clpp"], ascending=[True, False])
        r = g.iloc[0]
        prim.append(dict(gene=gene, gwas=gw, ra_lead=r.ra_lead_rsid, instrument=r.qtl_variant, dataset=f"{r.qtl_project} {r.biosample} "
                         f"{r.condition if isinstance(r.condition, str) else ''}".strip(), qtl_type=r.qtl_type, category=r.category,
                         clpp=r.clpp, OR_per_SD=r.OR_per_SD, OR_lo=r.OR_lo, OR_hi=r.OR_hi, z=r.z,
                         n_datasets=len(g), share_same_sign=float((np.sign(g.wald) == np.sign(r.wald)).mean())))
    pr = pd.DataFrame(prim).sort_values(["gene", "gwas"])
    pr.to_csv(TAB / "mr_primary.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
    print("\nFocus genes (primary instrument per signal):")
    print(pr[pr.gene.isin(FOCUS)].round(3).to_string(index=False))
    for gene, hyp in (("CTLA4", "M1 OR < 1"), ("IL6ST", "M2 OR > 1")):
        x = pr[pr.gene == gene]
        print(f"{hyp} {gene}: " + ("; ".join(f"{r.gwas} OR {r.OR_per_SD:.2f} ({r.OR_lo:.2f}-{r.OR_hi:.2f}), "
                                              f"{r.share_same_sign:.0%} of {r.n_datasets} datasets same sign" for r in x.itertuples()) or "no instrument"))

    f = pr[pr.gene.isin(FOCUS)].copy()
    f["label"] = (f.gene + "  (" + f.gwas.str.split().str[0] + "; " + f.category + "; "
                  + (f.share_same_sign * f.n_datasets).round().astype(int).astype(str) + "/" + f.n_datasets.astype(str) + " same sign)")
    f = f.sort_values("OR_per_SD")
    fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * len(f))))
    y = np.arange(len(f))
    ax.errorbar(f.OR_per_SD, y, xerr=[f.OR_per_SD - f.OR_lo, f.OR_hi - f.OR_per_SD], fmt="o", color="k", ecolor="0.5")
    ax.axvline(1, color="C3", lw=0.8, ls="--")
    ax.set_xscale("log"); ax.set_yticks(y, f.label, fontsize=8)
    ticks = [0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75]
    ax.set_xticks(ticks, [str(t) for t in ticks]); ax.minorticks_off()
    ax.set_xlabel("OR of RA per SD higher expression / protein (Wald ratio, 95 % CI)")
    ax.set_title("Drug-target MR at colocalising loci (primary instrument per signal)", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "drug_target_mr.png", dpi=150)


if __name__ == "__main__":
    main()
