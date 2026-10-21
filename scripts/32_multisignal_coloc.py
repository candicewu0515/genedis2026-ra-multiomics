"""Multi-signal colocalisation of fine-mapped RA credible sets with molecular QTLs, and direction of effect.

Step 25 used coloc.abf, which assumes one causal variant per locus; step 29 showed that IRF5 carries two independent
RA signals. Here each fine-mapped RA credible set (step 29: Ishigaki 2022 EUR, FinnGen R12; MHC excluded) is matched
to fine-mapped QTL credible sets by Open Targets Platform 26.09 (coloc on SuSiE credible sets plus eCAVIAR CLPP;
eQTL Catalogue incl. single-cell eQTL, GTEx v10, UKB-PPP pQTL), so each signal is tested separately.
Colocalisation: H4 >= 0.8, cis only (QTL gene's TSS within 1 Mb of the RA lead variant; UKB-PPP also reports
trans-pQTLs, which would otherwise attach dozens of distant proteins to every immune locus).
Direction: sign of beta(GWAS) / beta(QTL) for the same allele (Open Targets betaRatioSignAverage): +1 = the
risk-raising allele raises the gene's expression or protein level, -1 = it lowers it. Taken from eQTL / single-cell
eQTL / cis-pQTL only (splicing and transcript-usage QTLs do not measure level), one call per QTL dataset (the record
with the highest CLPP), then the consensus across datasets.
Biosamples are grouped by name: myeloid (monocyte, macrophage, dendritic cell), granulocyte, T/NK, B/LCL,
blood/plasma, other tissues.
Questions fixed before analysis:
  C1 does each of the two IRF5 signals colocalise with myeloid IRF5 expression, and with which direction?
  C2 how many RA signals colocalise with any QTL, and how many with an immune-cell QTL; which genes?
  C3 for prioritised genes (step 22), is the risk allele's effect on expression in immune cells consistent in sign,
     and does it agree with inhibiting the gene (late-stage inhibitors from step 15)?
  Replication: genes colocalising in both Ishigaki and FinnGen.

Output: data/raw/opentargets/coloc_<studyLocusId>.json (cached), results/tables/coloc_multisignal_all.tsv,
        coloc_multisignal_genes.tsv, coloc_direction.tsv; results/figures/coloc_multisignal.png
"""
import json
import re
import time
import urllib.request

import gzip

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import RAW, FIG, TAB

OT = "https://api.platform.opentargets.org/api/v4/graphql"
C = RAW / "opentargets"
STUDIES = {"GCST90132223": "Ishigaki 2022 (EUR)", "FINNGEN_R12_M13_RHEUMA": "FinnGen R12"}
H4 = 0.8
Q = """query($id: String!, $i: Int!) { credibleSet(studyLocusId: $id) { colocalisation(
  studyTypes: [eqtl, pqtl, sqtl, tuqtl, sceqtl, scpqtl, scsqtl, sctuqtl], page: {index: $i, size: 100}) {
  count rows { colocalisationMethod h4 clpp betaRatioSignAverage numberColocalisingVariants
  otherStudyLocus { studyLocusId studyType variant { id } study { id projectId condition
  target { approvedSymbol } biosample { biosampleName } } } } } } }"""
CIS = 1_000_000
LEVEL_QTL = {"eqtl", "sceqtl", "pqtl", "scpqtl"}
CATS = [("myeloid", r"monocyte|macrophage|dendritic|CD14|CD16"), ("granulocyte", r"neutrophil|granulocyte|eosinophil"),
        ("T/NK", r"T cell|T-cell|CD4|CD8|regulatory T|natural killer|NK cell|T-helper|thymocyte|lymphocyte"),
        ("B/LCL", r"B cell|B-cell|lymphoblastoid|plasma cell"), ("blood/plasma", r"blood|PBMC|mononuclear|platelet|plasma")]
IMMUNE = {"myeloid", "granulocyte", "T/NK", "B/LCL"}


def gql(variables):
    req = urllib.request.Request(OT, data=json.dumps({"query": Q, "variables": variables}).encode(),
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


def colocs(sl):
    p = C / f"coloc_{sl}.json"
    if p.exists():
        return json.loads(p.read_text())
    rows, i = [], 0
    while True:
        c = gql({"id": sl, "i": i})["credibleSet"]["colocalisation"]
        rows += c["rows"]
        if len(rows) >= c["count"] or not c["rows"]:
            break
        i += 1
    p.write_text(json.dumps(rows))
    return rows


def gene_tss():
    rows = []
    with gzip.open(RAW / "annotation" / "Homo_sapiens.GRCh38.98.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene":
                continue
            nm = re.search(r'gene_name "([^"]+)"', f[8]); gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            tss = int(f[3]) if f[6] == "+" else int(f[4])
            rows.append((gid, f[0], tss))
            if nm:
                rows.append((nm.group(1), f[0], tss))
    return pd.DataFrame(rows, columns=["gene", "chrom", "tss"]).drop_duplicates("gene").set_index("gene")


def per_dataset(df):
    """One direction call per QTL dataset (study x biosample x condition): the level-QTL record with the highest CLPP."""
    d = df[df.qtl_type.isin(LEVEL_QTL) & df.direction.notna() & (df.direction != 0)]
    if d.empty:
        return d
    return d.sort_values("clpp", ascending=False).drop_duplicates(["ra_set", "gene", "qtl_study"])


def category(name):
    for cat, pat in CATS:
        if name and re.search(pat, name, re.I):
            return cat
    return "other tissue"


def main():
    var = pd.read_csv(TAB / "finemap_variants.tsv", sep="\t")
    sets = var[var.study.isin(STUDIES)].sort_values("pip", ascending=False).drop_duplicates("studyLocusId")
    rows = []
    for _, s in sets.iterrows():
        for r in colocs(s.studyLocusId):
            o = r["otherStudyLocus"]; st = o["study"]
            bio = (st.get("biosample") or {}).get("biosampleName") or ""
            rows.append(dict(gwas=STUDIES[s.study], ra_set=s.studyLocusId, ra_lead=s.lead, ra_lead_rsid=s.rsid,
                             qtl_type=o["studyType"], qtl_project=st["projectId"], qtl_study=st["id"],
                             gene=(st.get("target") or {}).get("approvedSymbol"), biosample=bio,
                             condition=st.get("condition"), category=category(bio), h4=r["h4"], clpp=r["clpp"],
                             direction=r["betaRatioSignAverage"], n_coloc_variants=r["numberColocalisingVariants"]))
    allc = pd.DataFrame(rows)
    tss = gene_tss()
    lead_chr = allc.ra_lead.str.split("_").str[0]; lead_pos = allc.ra_lead.str.split("_").str[1].astype(int)
    gchr = allc.gene.map(tss.chrom); gpos = allc.gene.map(tss.tss)
    allc["cis"] = (gchr == lead_chr) & ((gpos - lead_pos).abs() <= CIS)
    allc.to_csv(TAB / "coloc_multisignal_all.tsv", sep="\t", index=False, float_format="%.4g")
    hit = allc[(allc.h4 >= H4) & allc.gene.notna() & allc.cis].copy()
    print(f"RA credible sets queried: {len(sets)}; QTL colocalisation records: {len(allc)}; with H4 >= {H4}: "
          f"{int((allc.h4 >= H4).sum())}; of these cis: {len(hit)} (trans or unplaced dropped: {int(((allc.h4 >= H4) & ~allc.cis).sum())})")

    # gene-level summary per RA signal
    g = hit.groupby(["gwas", "ra_set", "ra_lead_rsid", "gene"]).agg(
        n_records=("h4", "size"), best_h4=("h4", "max"),
        categories=("category", lambda x: ",".join(sorted(set(x)))),
        immune_records=("category", lambda x: int(x.isin(IMMUNE).sum())),
        myeloid_biosamples=("biosample", lambda x: ""),
    ).reset_index()
    for k, r in g.iterrows():
        sub = hit[(hit.ra_set == r.ra_set) & (hit.gene == r.gene)]
        my = sub[sub.category == "myeloid"]
        g.loc[k, "myeloid_biosamples"] = "; ".join(sorted(set(my.qtl_project + " " + my.biosample + " " + my.condition.fillna(""))))[:200]
        imm = per_dataset(sub[sub.category.isin(IMMUNE)])
        g.loc[k, "immune_direction_mean"] = imm.direction.mean() if len(imm) else np.nan
        g.loc[k, "immune_direction_datasets"] = len(imm)
    g = g.sort_values(["immune_records", "best_h4"], ascending=False)
    g.to_csv(TAB / "coloc_multisignal_genes.tsv", sep="\t", index=False, float_format="%.3g")

    sets_with_any = hit.ra_set.nunique()
    sets_with_imm = hit[hit.category.isin(IMMUNE)].ra_set.nunique()
    print(f"C2: RA signals with any QTL colocalisation {sets_with_any}/{len(sets)}; with an immune-cell QTL {sets_with_imm}/{len(sets)}")
    for gw in STUDIES.values():
        x = g[(g.gwas == gw) & (g.immune_records > 0)]
        print(f"  {gw}: genes with immune-cell colocalisation ({x.gene.nunique()}): {', '.join(sorted(x.gene.unique()))}")
    both = set(g[g.gwas == STUDIES['GCST90132223']].gene) & set(g[g.gwas == STUDIES['FINNGEN_R12_M13_RHEUMA']].gene)
    print(f"Replication: genes colocalising in both studies ({len(both)}): {', '.join(sorted(both))}")

    # C1 IRF5
    print("\nC1 IRF5 signals:")
    irf5 = hit[(hit.gene == "IRF5")]
    for (gw, rs), sub in irf5.groupby(["gwas", "ra_lead_rsid"]):
        my = sub[sub.category == "myeloid"]
        print(f"  {gw} {rs}: {len(sub)} IRF5 QTL records (H4 >= {H4}); myeloid {len(my)}: "
              + "; ".join(f"{r.qtl_project} {r.biosample} {r.condition if isinstance(r.condition, str) else ''} H4 {r.h4:.2f} dir {r.direction:+.0f}"
                          for r in my.sort_values("h4", ascending=False).itertuples()))
    for (gw, rs), sub in allc[allc.ra_lead_rsid.isin(irf5.ra_lead_rsid.unique()) & (allc.gene == "IRF5") &
                               (allc.category == "myeloid") & (allc.h4 < H4)].groupby(["gwas", "ra_lead_rsid"]):
        print(f"  {gw} {rs}: myeloid IRF5 records below threshold: {len(sub)} (max H4 {sub.h4.max():.2f})")

    # C3 direction for prioritised genes
    pr = pd.read_csv(TAB / "gene_prioritisation.tsv", sep="\t", index_col=0)
    mine = pd.read_csv(TAB / "target_drug_mining.tsv", sep="\t").set_index("gene")
    rows = []
    for gene, sub in per_dataset(hit[hit.category.isin(IMMUNE)]).groupby("gene"):
        d = sub.direction
        frac_up = float((d > 0).mean())
        rows.append(dict(gene=gene, immune_datasets=len(sub), signals=sub.ra_set.nunique(), frac_risk_raises_expression=frac_up,
                         consistent=bool(frac_up >= 0.8 or frac_up <= 0.2),
                         risk_allele_effect=("raises expression" if frac_up >= 0.8 else "lowers expression" if frac_up <= 0.2
                                             else "mixed") if len(sub) >= 2 else "single dataset",
                         in_prioritisation=gene in pr.index, evidence_score=pr.evidence_score.get(gene, np.nan),
                         late_stage_inhibitors=mine.late_stage_inhibitors.get(gene, "") if gene in mine.index else ""))
    dr = pd.DataFrame(rows).sort_values(["in_prioritisation", "immune_datasets"], ascending=False)
    dr["inhibition_fits"] = np.where(dr.risk_allele_effect == "raises expression", "yes",
                                     np.where(dr.risk_allele_effect == "lowers expression", "no (risk lowers expression)", "unclear"))
    dr.to_csv(TAB / "coloc_direction.tsv", sep="\t", index=False, float_format="%.3g")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
    print("\nC3 direction in immune cells (genes with immune-cell colocalisation):")
    print(dr.head(30).to_string(index=False))

    # figure: genes x category (best H4), direction annotated
    imm_genes = dr.sort_values(["in_prioritisation", "immune_datasets"], ascending=False).gene.tolist()
    top = [x for x in imm_genes if not x.startswith("ENSG")][:30]
    cats = [c for c, _ in CATS] + ["other tissue"]
    M = pd.DataFrame(0.0, index=top, columns=cats)
    for gene in top:
        sub = hit[hit.gene == gene]
        for c in cats:
            v = sub[sub.category == c].h4
            M.loc[gene, c] = v.max() if len(v) else 0
    fig, axes = plt.subplots(1, 2, figsize=(14, max(5, 0.3 * len(top))), gridspec_kw={"width_ratios": [2, 1]})
    ax = axes[0]
    ax.imshow(M.values, aspect="auto", cmap="Reds", vmin=0, vmax=1)
    ax.set_xticks(range(len(cats)), cats, rotation=30, ha="right"); ax.set_yticks(range(len(top)), top, fontsize=7)
    ax.set_title(f"RA signals x cis-QTLs: best H4 by cell group (prioritised genes first)", fontsize=10)
    ax = axes[1]
    d = dr.set_index("gene").reindex(top)
    ax.barh(range(len(top)), d.frac_risk_raises_expression.fillna(0.5) - 0.5,
            color=["C3" if x == "raises expression" else "C0" if x == "lowers expression" else "C7" for x in d.risk_allele_effect])
    ax.axvline(0, color="k", lw=0.8); ax.set_xlim(-0.55, 0.55)
    ax.set_yticks(range(len(top)), top, fontsize=7)
    ax.set_xlabel("risk allele raises level: share of immune QTL datasets - 0.5", fontsize=8)
    ax.set_title("Direction (red: raises, blue: lowers, grey: mixed)", fontsize=10)
    axes[0].set_ylim(len(top) - 0.5, -0.5); ax.set_ylim(len(top) - 0.5, -0.5)
    fig.tight_layout(); fig.savefig(FIG / "coloc_multisignal.png", dpi=150)


if __name__ == "__main__":
    main()
