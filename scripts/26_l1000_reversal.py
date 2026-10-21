"""Signature-reversal drug search (L1000CDS2, LINCS L1000) for the active-RA macrophage programme.

Queries (reverse mode, top-50 signatures each):
  Q1 primary: up = 15 FDR genes up in active-RA macrophages + next-ranked up genes to 50 (E-MTAB-8322 Wald statistic);
     down = top 50 genes up in remission.
  Q2 robustness: up = Hallmark IFN-gamma leading edge shared by both GSEA contrasts; down = a priori homeostatic
     macrophage programme.
Compounds are aggregated across signatures (count, best score) and annotated with ChEMBL mechanisms of action.
Positive-control classes (expected if the programme is IFN/JAK-STAT and NF-kB driven): JAK, glucocorticoid receptor,
HSP90, proteasome, NF-kB/IKK, BET, HDAC, mTOR/PI3K. L1000 profiles come from cell lines, not synovial macrophages.

Output: data/raw/l1000cds2/*.json (cached), results/tables/l1000_reversal_signatures.tsv, l1000_reversal_compounds.tsv
"""
import json
import re
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

from config import RAW, TAB, MAC_HOMEOSTATIC

API = "https://maayanlab.cloud/L1000CDS2/query"
CH = "https://www.ebi.ac.uk/chembl/api/data"
C = RAW / "l1000cds2"
C.mkdir(exist_ok=True)
# Tool compounds often lack ChEMBL mechanism records, so classes are also assigned from well-known compound names.
NAME_CLASSES = {"JAK": r"ruxolitinib|tofacitinib|baricitinib|upadacitinib|filgotinib|momelotinib|fedratinib|AZD-?1480|TG-?101348|CP-?690550",
                "glucocorticoid": r"dexamethasone|prednisolone|prednisone|methylprednisolone|hydrocortisone|budesonide|betamethasone|triamcinolone|fluocinolone|clobetasol|halcinonide|beclomethasone",
                "HSP90": r"geldanamycin|radicicol|tanespimycin|alvespimycin|17-?AAG|17-?DMAG|NVP-?AUY922|luminespib|ganetespib|PU-?H71|SNX-?2112|BIIB021",
                "proteasome": r"bortezomib|carfilzomib|MG-?132|ixazomib|epoxomicin|MG-?262",
                "NF-kB/IKK": r"BAY ?11-?7(082|085|821)|withaferin|parthenolide|IKK|BMS-?345541|TPCA-?1|SC-?514",
                "BET": r"I-?BET|JQ-?1|PFI-?1\b|OTX-?015|birabresib|molibresib",
                "HDAC": r"vorinostat|trichostatin|panobinostat|belinostat|entinostat|mocetinostat|romidepsin|givinostat|HDAC",
                "mTOR/PI3K": r"sirolimus|rapamycin|temsirolimus|everolimus|torin|AZD-?8055|PP-?242|NVP-?BEZ235|BGT226|GDC-?0941|wortmannin|LY-?294002|PI-?103"}
CLASSES = {"JAK": r"\bJAK|Janus", "glucocorticoid": r"[Gg]lucocorticoid", "HSP90": r"HSP ?90|Heat shock protein 90",
           "proteasome": r"[Pp]roteasome", "NF-kB/IKK": r"NF-?kappa|IKK|Inhibitor of nuclear factor kappa",
           "BET": r"BET|Bromodomain", "HDAC": r"[Hh]istone deacetylase|HDAC", "mTOR/PI3K": r"mTOR|PI3|[Pp]hosphatidylinositol"}


def cached_json(name, fn):
    p = C / f"{name}.json"
    if p.exists():
        return json.loads(p.read_text())
    d = fn()
    p.write_text(json.dumps(d))
    return d


def l1000(up, dn, name):
    payload = {"data": {"upGenes": up, "dnGenes": dn},
               "config": {"aggravate": False, "searchMethod": "geneSet", "share": False, "combination": False,
                          "db-version": "latest"}}

    def call():
        req = urllib.request.Request(API, data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=180) as r:
            return json.load(r)
    return cached_json(f"query_{name}", call)


def chembl_moa(name):
    def call():
        out = {"query": name, "chembl_id": None, "mechanisms": []}
        q = urllib.parse.quote(name)
        try:
            with urllib.request.urlopen(f"{CH}/molecule/search.json?q={q}&limit=1", timeout=60) as r:
                mols = json.load(r).get("molecules", [])
            if not mols:
                return out
            cid = mols[0]["molecule_chembl_id"]
            out["chembl_id"] = cid
            out["pref_name"] = mols[0].get("pref_name")
            with urllib.request.urlopen(f"{CH}/mechanism.json?molecule_chembl_id={cid}&limit=20", timeout=60) as r:
                out["mechanisms"] = [m.get("mechanism_of_action") for m in json.load(r).get("mechanisms", [])]
        except Exception as e:
            out["error"] = str(e)
        time.sleep(0.3)
        return out
    safe = re.sub(r"[^A-Za-z0-9]+", "_", name)[:60]
    return cached_json(f"chembl_{safe}", call)


def main():
    de = pd.read_csv(TAB / "emtab8322_DE_active_vs_remission.tsv", sep="\t", index_col=0).dropna(subset=["stat"])
    de = de[~de.index.str.match(r"^(ENSG|MT-)")]
    fdr_up = de[(de.padj < 0.05) & (de.log2FoldChange > 0)].index.tolist()
    up1 = (fdr_up + [g for g in de.sort_values("stat", ascending=False).index if g not in fdr_up])[:50]
    dn1 = de.sort_values("stat").index[:50].tolist()
    gs = pd.read_csv(TAB / "gsea_results.tsv", sep="\t")
    le = [set(str(x).split(";")) for x in gs[(gs.pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE") & (gs.library == "hallmark")]["Lead_genes"]]
    up2 = sorted(set.intersection(*le))
    dn2 = list(MAC_HOMEOSTATIC) + ["SELENOP", "LYVE1", "DAB2", "SLC40A1", "MAF"]

    rows = []
    for name, up, dn in (("Q1_active_vs_remission", up1, dn1), ("Q2_IFNg_vs_homeostatic", up2, dn2)):
        r = l1000(up, dn, name)
        meta = r.get("topMeta", [])
        print(f"{name}: up {len(up)} genes, down {len(dn)} genes; signatures returned {len(meta)}")
        for i, m in enumerate(meta):
            rows.append(dict(query=name, rank=i + 1, compound=m.get("pert_desc"), pert_id=m.get("pert_id"),
                             cell=m.get("cell_id"), dose=m.get("pert_dose"), dose_unit=m.get("pert_dose_unit"),
                             time=m.get("pert_time"), score=m.get("score"), sig_id=m.get("sig_id")))
    sig = pd.DataFrame(rows)
    sig.to_csv(TAB / "l1000_reversal_signatures.tsv", sep="\t", index=False)

    agg = sig[sig.compound.notna() & (sig.compound != "-666")].groupby(["query", "compound"]).agg(
        n_signatures=("sig_id", "size"), best_rank=("rank", "min"), best_score=("score", "max"),
        cells=("cell", lambda s: ",".join(sorted(set(s))))).reset_index()
    moa = {c: chembl_moa(c) for c in sorted(agg.compound.unique())}
    agg["chembl_id"] = agg.compound.map(lambda c: moa[c].get("chembl_id"))
    agg["mechanism_of_action"] = agg.compound.map(lambda c: "; ".join(m for m in moa[c].get("mechanisms", []) if m))
    for k, pat in CLASSES.items():
        agg[f"class_{k}"] = agg.mechanism_of_action.str.contains(pat, regex=True) | \
            agg.compound.str.contains(NAME_CLASSES[k], case=False, regex=True)
    agg["classes"] = agg.apply(lambda r: ",".join(k for k in CLASSES if r[f"class_{k}"]), axis=1)
    agg["in_both_queries"] = agg.compound.isin(set(agg[agg["query"].str.startswith("Q1")].compound) &
                                               set(agg[agg["query"].str.startswith("Q2")].compound))
    agg = agg.sort_values(["query", "best_rank"])
    base = random_baseline(de, 50, 50)
    obs = {q: float((agg[agg["query"] == q].classes != "").mean()) for q in agg["query"].unique()}
    print(f"Random gene lists (n = {len(base)}): share of named top reversers in positive-control classes "
          f"median {np.median(base):.2f} (range {base.min():.2f}-{base.max():.2f}); observed "
          + ", ".join(f"{q} {v:.2f} (empirical p = {(1 + (base >= v).sum()) / (1 + len(base)):.3f})" for q, v in obs.items()))
    pd.DataFrame([dict(query=q, observed_share=v, random_median=float(np.median(base)), random_max=float(base.max()),
                       empirical_p=(1 + (base >= v).sum()) / (1 + len(base))) for q, v in obs.items()]) \
        .to_csv(TAB / "l1000_reversal_class_enrichment.tsv", sep="\t", index=False, float_format="%.4g")
    agg.to_csv(TAB / "l1000_reversal_compounds.tsv", sep="\t", index=False)
    pd.set_option("display.width", 240); pd.set_option("display.max_colwidth", 70)
    for q in agg["query"].unique():
        a = agg[agg["query"] == q]
        print(f"\n{q}: {len(a)} named compounds; annotated with a ChEMBL mechanism: {(a.mechanism_of_action != '').sum()}")
        print("  positive-control classes hit:", {k: int(a[f'class_{k}'].sum()) for k in CLASSES})
        print(a.head(15)[["compound", "n_signatures", "best_rank", "cells", "classes", "in_both_queries"]].to_string(index=False))
        cls = a[a.classes != ""]
        print(f"  compounds in positive-control classes: {len(cls)}/{len(a)} "
              f"({cls.n_signatures.sum()} of {a.n_signatures.sum()} named signatures); best ranks: "
              + ", ".join(f"{c} ({k}, #{r})" for c, k, r in zip(cls.compound, cls.classes, cls.best_rank)))


def random_baseline(de, n_up, n_dn, n_rep=20):
    """Base rate of positive-control classes among top reversers for random gene lists of the same sizes."""
    import numpy as np
    rng = np.random.default_rng(0)
    pool = de[de.baseMean > 20].index.tolist()
    fr = []
    for i in range(n_rep):
        genes = list(rng.choice(pool, n_up + n_dn, replace=False))
        r = l1000(genes[:n_up], genes[n_up:], f"random_{i:02d}")
        names = [m.get("pert_desc") for m in r.get("topMeta", []) if m.get("pert_desc") not in (None, "-666")]
        names = sorted(set(names))
        hit = [n for n in names if any(re.search(NAME_CLASSES[k], n, re.I) for k in NAME_CLASSES)]
        fr.append(len(hit) / max(len(names), 1))
    return np.array(fr)


if __name__ == "__main__":
    main()
