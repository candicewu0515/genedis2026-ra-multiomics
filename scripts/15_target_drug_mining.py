"""Knowledge mining: what is already known and druggable about the genes this workflow surfaced?

Candidates are taken programmatically from earlier outputs (no hand picking):
  signature      FDR genes up in active-RA macrophages (step 05)
  gwas           nearest protein-coding gene of RA lead SNPs inside state-specific peaks (step 10)
  communication  IFN-gamma and CXCR3 axes + ligands/receptors of the top-20 M-6-specific interactions (step 13)
  pathway        leading-edge genes of HALLMARK_INTERFERON_GAMMA_RESPONSE shared by both GSEA contrasts (step 12)
  control        targets of approved RA drugs (IL6R, JAK1, TNF, MS4A1); must come back as "approved in RA"
Open Targets Platform GraphQL API (data release recorded in the output): RA association score
(MONDO_0008383; overall and genetic-association datatype), tractability, drugs and clinical candidates.

Output: data/raw/opentargets/*.json (cached responses), results/tables/target_drug_mining.tsv,
        results/figures/target_drug_mining.png
"""
import gzip
import json
import re
import time
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import RAW, FIG, TAB

API = "https://api.platform.opentargets.org/api/v4/graphql"
RA_ID = "MONDO_0008383"
OT = RAW / "opentargets"
OT.mkdir(exist_ok=True)
CONTROLS = ["IL6R", "JAK1", "TNF", "MS4A1"]
STAGE_ORDER = ["APPROVAL", "PHASE_4", "PHASE_3", "PHASE_2_3", "PHASE_2", "PHASE_1_2", "PHASE_1", "EARLY_PHASE_1",
               "IND", "PRECLINICAL", "UNKNOWN"]
RA_DISEASE = re.compile(r"rheumatoid arthritis", re.I)


def gql(query, cache_name):
    path = OT / f"{cache_name}.json"
    if path.exists():
        return json.loads(path.read_text())
    req = urllib.request.Request(API, data=json.dumps({"query": query}).encode(),
                                 headers={"Content-Type": "application/json"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                d = json.load(r)
            break
        except Exception:
            if attempt == 4:
                raise
            time.sleep(5 * (attempt + 1))
    if "errors" in d:
        raise RuntimeError(d["errors"])
    path.write_text(json.dumps(d))
    time.sleep(0.2)
    return d


def stage_rank(s):
    s = (s or "UNKNOWN").upper().replace("/", "_")
    return STAGE_ORDER.index(s) if s in STAGE_ORDER else len(STAGE_ORDER)


def candidates():
    rows = []
    sg = pd.read_csv(TAB / "signature_genes.tsv", sep="\t")
    for g in sg.loc[(sg.signature == "primary_FDR") & (sg.direction == "up_in_active"), "gene"]:
        rows.append((g, "signature", "up in active-RA macrophages"))
    cv = pd.read_csv(TAB / "gwas_candidate_variants.tsv", sep="\t")
    for g, s in cv.groupby("nearest_protein_coding_gene"):
        rows.append((g, "gwas", ";".join(sorted(set(s.rs_id))) + " in " + "; ".join(sorted(set(s.state)))[:120]))
    for g in ["IFNG", "IFNGR1", "IFNGR2", "CXCL9", "CXCL10", "CXCL11", "CXCR3"]:
        rows.append((g, "communication", "IFN-gamma / CXCR3 axis around M-6"))
    sp = pd.read_csv(TAB / "liana_to_M6_specific.tsv", sep="\t").head(20)
    for _, r in sp.iterrows():
        for g in str(r.ligand_complex).split("_"):
            rows.append((g, "communication", f"ligand {r.ligand_complex} from {r.source.split(':')[0]} -> M-6"))
        for g in str(r.receptor_complex).split("_"):
            rows.append((g, "communication", f"receptor {r.receptor_complex} on M-6"))
    gs = pd.read_csv(TAB / "gsea_results.tsv", sep="\t")
    ifn = gs[(gs.pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE") & (gs.library == "hallmark")]
    le = [set(str(x).split(";")) for x in ifn["Lead_genes"]]
    for g in sorted(set.intersection(*le)) if len(le) == 2 else []:
        rows.append((g, "pathway", "IFN-gamma response leading edge in both GSEA contrasts"))
    for g in CONTROLS:
        rows.append((g, "control", "target of an approved RA drug"))
    df = pd.DataFrame(rows, columns=["gene", "source", "evidence"])
    return df.groupby("gene").agg(sources=("source", lambda s: ",".join(sorted(set(s)))),
                                  evidence=("evidence", lambda s: " | ".join(dict.fromkeys(s))[:300])).reset_index()


def symbol_to_ensembl():
    m = {}
    with gzip.open(RAW / "annotation" / "Homo_sapiens.GRCh38.98.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene":
                continue
            name = re.search(r'gene_name "([^"]+)"', f[8]); gid = re.search(r'gene_id "([^"]+)"', f[8])
            if name and name.group(1) not in m:
                m[name.group(1)] = gid.group(1)
    return m


def ra_associations():
    rows, page = [], 0
    while True:
        d = gql(f'query {{ disease(efoId: "{RA_ID}") {{ associatedTargets(page: {{index: {page}, size: 1000}}) {{ count '
                f'rows {{ target {{ id }} score datatypeScores {{ id score }} }} }} }} }}', f"ra_assoc_page{page}")
        at = d["data"]["disease"]["associatedTargets"]
        for r in at["rows"]:
            dt = {x["id"]: x["score"] for x in r["datatypeScores"]}
            rows.append(dict(ensembl=r["target"]["id"], ot_ra_score=r["score"],
                             ot_ra_genetic=dt.get("genetic_association", 0.0), ot_ra_clinical=dt.get("clinical", 0.0)))
        if (page + 1) * 1000 >= at["count"]:
            break
        page += 1
    a = pd.DataFrame(rows).set_index("ensembl")
    a["ot_ra_rank"] = a["ot_ra_score"].rank(ascending=False).astype(int)
    return a, at["count"]


INHIBITORY = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR", "NEGATIVE MODULATOR",
              "ANTISENSE INHIBITOR", "RNAI INHIBITOR", "DEGRADER", "INVERSE AGONIST"}
ACTIVATING = {"AGONIST", "ACTIVATOR", "PARTIAL AGONIST", "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE MODULATOR",
              "OPENER", "STABILISER"}


def drug_info(chembl, ens):
    """RA-specific stage (from the drug's indications) and the drug's action type on this target."""
    d = gql(f'query {{ drug(chemblId: "{chembl}") {{ indications {{ rows {{ maxClinicalStage disease {{ name }} }} }} '
            f'mechanismsOfAction {{ rows {{ actionType targets {{ id }} }} }} }} }}', f"drug2_{chembl}")
    drug = d["data"]["drug"] or {}
    rows = (drug.get("indications") or {}).get("rows") or []
    ra = [r["maxClinicalStage"] for r in rows if r["disease"] and RA_DISEASE.fullmatch(r["disease"]["name"].strip())]
    acts = {m["actionType"] for m in ((drug.get("mechanismsOfAction") or {}).get("rows") or [])
            if m["actionType"] and any(t["id"] == ens for t in m["targets"])}
    action = "INHIBITORY" if acts & INHIBITORY else ("ACTIVATING" if acts & ACTIVATING else ("OTHER" if acts else "UNKNOWN"))
    return (min(ra, key=stage_rank) if ra else None), action


def target_info(ens):
    d = gql(f'query {{ target(ensemblId: "{ens}") {{ approvedSymbol tractability {{ label modality value }} '
            f'drugAndClinicalCandidates {{ count rows {{ maxClinicalStage drug {{ id name drugType }} }} }} }} }}',
            f"target2_{ens}")
    t = d["data"]["target"]
    if t is None:
        return {}
    tract = [x for x in t["tractability"] if x["value"]]
    rows = t["drugAndClinicalCandidates"]["rows"] if t["drugAndClinicalCandidates"] else []
    best = min(rows, key=lambda r: stage_rank(r["maxClinicalStage"])) if rows else None
    ra, late_inh, late_act = [], set(), set()
    for r in rows:
        st, action = drug_info(r["drug"]["id"], ens)
        if st:
            ra.append((st, r["drug"]["name"]))
        if stage_rank(r["maxClinicalStage"]) <= 2:          # approval, phase 4, phase 3
            (late_inh if action == "INHIBITORY" else late_act if action == "ACTIVATING" else set()).add(r["drug"]["name"])
    best_ra = min(ra, key=lambda x: stage_rank(x[0])) if ra else None
    return dict(ot_symbol=t["approvedSymbol"],
                tractable_small_molecule=any(x["modality"] == "SM" for x in tract),
                tractable_antibody=any(x["modality"] == "AB" for x in tract),
                n_drugs=len({r["drug"]["name"] for r in rows}),
                max_stage_any=best["maxClinicalStage"] if best else "",
                late_stage_inhibitors="; ".join(sorted(late_inh))[:200],
                late_stage_activators="; ".join(sorted(late_act))[:200],
                max_stage_RA=best_ra[0] if best_ra else "",
                RA_drugs="; ".join(f"{n} ({st})" for st, n in sorted(ra, key=lambda x: stage_rank(x[0])))[:250])


def main():
    cand = candidates()
    ens = symbol_to_ensembl()
    cand["ensembl"] = cand["gene"].map(ens)
    print(f"Candidates: {len(cand)} genes ({cand.ensembl.isna().sum()} without Ensembl ID)")
    meta = gql("query { meta { apiVersion { x y z } dataVersion { year month } } }", "meta")["data"]["meta"]
    release = f"{meta['dataVersion']['year']}.{meta['dataVersion']['month']}"
    assoc, n_assoc = ra_associations()
    info = pd.DataFrame([dict(ensembl=e, **target_info(e)) for e in cand["ensembl"].dropna()]).set_index("ensembl")
    out = cand.set_index("ensembl", drop=False).join(assoc).join(info)
    out["ot_ra_score"] = out["ot_ra_score"].fillna(0.0)
    for c in ("max_stage_RA", "max_stage_any", "RA_drugs"):
        out[c] = out[c].fillna("")
    for c in ("late_stage_inhibitors", "late_stage_activators"):
        out[c] = out[c].fillna("")
    disease_up = out["sources"].str.contains("signature|pathway|communication")
    has_inh = out["late_stage_inhibitors"] != ""
    late_any = out["max_stage_any"].isin(["APPROVAL", "PHASE_4", "PHASE_3"])
    tract = out["tractable_small_molecule"].fillna(False).astype(bool) | out["tractable_antibody"].fillna(False).astype(bool)
    out["category"] = np.select(
        [out["max_stage_RA"] == "APPROVAL",
         out["max_stage_RA"] != "",
         has_inh & disease_up,
         has_inh,
         late_any,
         tract],
        ["approved in RA", "tested in RA, not approved",
         "repurposing lead: late-stage inhibitor of a disease-up gene",
         "late-stage inhibitor, direction to verify (GWAS gene)",
         "late-stage drug, activating/unclear direction",
         "tractable, no late-stage drug"], default="no tractability evidence")
    out = out.sort_values(["category", "ot_ra_score"], ascending=[True, False])
    cols = ["gene", "sources", "category", "ot_ra_score", "ot_ra_rank", "ot_ra_genetic", "ot_ra_clinical",
            "max_stage_RA", "RA_drugs", "max_stage_any", "n_drugs", "late_stage_inhibitors", "late_stage_activators",
            "tractable_small_molecule",
            "tractable_antibody", "evidence", "ensembl"]
    out[cols].round(4).to_csv(TAB / "target_drug_mining.tsv", sep="\t", index=False)
    (OT / "README.txt").write_text(f"Open Targets Platform GraphQL API {API}; data release {release}; "
                                   f"RA = {RA_ID} ({n_assoc} associated targets). Responses cached as JSON here.\n")
    print(f"Open Targets data release {release}; RA associated targets: {n_assoc}")
    print("\nControls:")
    print(out[out.sources.str.contains("control")][["gene", "category", "max_stage_RA", "RA_drugs"]].to_string(index=False))
    print("\nCategory counts (non-control):")
    print(out[~out.sources.str.contains("control")]["category"].value_counts().to_string())
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
    print("\nNon-control genes with RA clinical evidence or late-stage drugs:")
    show = out[~out.sources.str.contains("control") & out.category.str.match(r"approved|tested|late|repurposing")]
    print(show[["gene", "sources", "category", "ot_ra_score", "RA_drugs", "late_stage_inhibitors", "late_stage_activators"]]
          .round(3).to_string(index=False))

    # Figure: RA association score vs best clinical stage, coloured by source
    fig, ax = plt.subplots(figsize=(11, 7))
    colors = {"signature": "C3", "gwas": "C0", "communication": "C2", "pathway": "C4", "control": "k"}
    for _, r in out.iterrows():
        src = r["sources"].split(",")[0]
        yy = (len(STAGE_ORDER) - stage_rank(r["max_stage_any"])) if r["max_stage_any"] else 0
        ax.scatter(r["ot_ra_score"], yy + np.random.default_rng(hash(r["gene"]) % 2**32).uniform(-.25, .25),
                   s=45, color=colors.get(src, "grey"), marker="*" if r["category"] == "approved in RA" else "o")
        if r["ot_ra_score"] > 0.3 or yy >= len(STAGE_ORDER) - 2 or src == "control":
            ax.annotate(r["gene"], (r["ot_ra_score"], yy), fontsize=6, xytext=(3, 2), textcoords="offset points")
    ax.set_yticks(range(0, len(STAGE_ORDER) + 1), ["no drug"] + STAGE_ORDER[::-1], fontsize=7)
    ax.set_xlabel(f"Open Targets RA association score (release {release})")
    ax.set_ylabel("Most advanced drug for the target (any indication)")
    for k, c in colors.items():
        ax.scatter([], [], color=c, label=k)
    ax.scatter([], [], color="grey", marker="*", label="approved in RA")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("Workflow-derived candidate genes: existing evidence and druggability", fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / "target_drug_mining.png", dpi=150)


if __name__ == "__main__":
    main()
