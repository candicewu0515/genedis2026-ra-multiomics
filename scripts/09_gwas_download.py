"""Fetch genome-wide significant associations from the NHGRI-EBI GWAS Catalog (REST API v2).

Traits: rheumatoid arthritis (test) and LDL cholesterol (negative control: hepatic/lipid biology,
not expected to enrich in synovial chromatin). Coordinates are GRCh38, as served by the Catalog.

Output: data/raw/gwas/<trait>_associations.tsv (one row per association x SNP), data/raw/gwas/README.txt
"""
import json
import time
import urllib.parse
import urllib.request
from datetime import date

import pandas as pd

from config import RAW

API = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
TRAITS = {"RA": "rheumatoid arthritis", "LDL": "low density lipoprotein cholesterol measurement"}
OUT = RAW / "gwas"
OUT.mkdir(parents=True, exist_ok=True)


def fetch(key, trait, size=50):
    """Page through the API; every page is cached so an interrupted download resumes where it stopped."""
    cache = OUT / f"cache_{key}"
    cache.mkdir(exist_ok=True)
    rows, page = [], 0
    while True:
        cf = cache / f"page_{page:05d}.json"
        if cf.exists():
            d = json.loads(cf.read_text())
        else:
            url = f"{API}?{urllib.parse.urlencode({'efo_trait': trait, 'size': size, 'page': page})}"
            for attempt in range(6):
                try:
                    with urllib.request.urlopen(url, timeout=300) as r:
                        d = json.load(r)
                    break
                except Exception:                        # the API is slow / flaky: back off and retry
                    if attempt == 5:
                        raise
                    time.sleep(10 * (attempt + 1))
            cf.write_text(json.dumps(d))
            time.sleep(0.5)
        for a in d.get("_embedded", {}).get("associations", []):
            snps = a.get("snp_allele") or []
            locs = a.get("locations") or []
            for i, s in enumerate(snps):
                loc = locs[i] if i < len(locs) else (locs[0] if len(locs) == 1 else "")
                rows.append(dict(association_id=a["association_id"], rs_id=s.get("rs_id"), effect_allele=s.get("effect_allele"),
                                 location=loc,
                                 p_value=a.get("p_value"), efo_traits=";".join(t["efo_trait"] for t in a.get("efo_traits", [])),
                                 efo_ids=";".join(t["efo_id"] for t in a.get("efo_traits", [])),
                                 reported_trait=";".join(a.get("reported_trait") or []),
                                 mapped_genes=";".join(a.get("mapped_genes") or []),
                                 accession_id=a.get("accession_id"), pubmed_id=a.get("pubmed_id"),
                                 first_author=a.get("first_author")))
        pg = d["page"]
        print(f"  {trait}: page {pg['number'] + 1}/{pg['totalPages']}", flush=True)
        if pg["number"] + 1 >= pg["totalPages"]:
            break
        page += 1
    return pd.DataFrame(rows)


def main():
    notes = [f"Downloaded {date.today().isoformat()} from {API} (NHGRI-EBI GWAS Catalog, REST API v2).",
             "Query parameter efo_trait=<trait label>; one row per association x SNP; locations are GRCh38 chr:pos."]
    for key, trait in TRAITS.items():
        out = OUT / f"{key}_associations.tsv"
        if out.exists() and "effect_allele" in pd.read_csv(out, sep="\t", nrows=1).columns:
            df = pd.read_csv(out, sep="\t")
        else:
            df = fetch(key, trait)
            df.to_csv(out, sep="\t", index=False)
        notes.append(f"{key}: efo_trait='{trait}' -> {df['association_id'].nunique()} associations, "
                     f"{df['rs_id'].nunique()} unique rsIDs, {df['accession_id'].nunique()} studies")
        print(notes[-1])
    # Cached reference provenance records its original retrieval date.
    if not (OUT / "README.txt").exists():
        (OUT / "README.txt").write_text("\n".join(notes) + "\n")


if __name__ == "__main__":
    main()
