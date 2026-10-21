"""Variant-to-function: do RA lead SNPs inside cell-state-specific peaks change TF motif matches?

For each RA lead SNP within 250 bp of a state-specific peak (step 10):
  - alleles from Ensembl REST (GRCh38 allele_string; first allele = reference), checked against UCSC hg38;
  - risk (effect) allele from the GWAS Catalog association (step 09);
  - every JASPAR 2026 CORE vertebrate motif is scored (log-odds, both strands) on all windows overlapping the SNP
    for the reference and each alternative allele; a motif "hit" is score >= MOODS threshold at p < 1e-4;
  - a motif effect is reported when the hit exists for one allele only and |delta score| >= 2;
  - the TF must be expressed in the peak's cell state (detected in >= 10 % of nuclei of that state, GSE243917).
These are sequence-based predictions, not measured allelic effects.

Input : results/tables/gwas_candidate_variants.tsv, data/raw/gwas/RA_associations.tsv, data/interim/hg38.fa,
        data/interim/gse243917_annotated.h5ad
Output: data/raw/ensembl/variants.json (cached), results/tables/variant_motif_effects.tsv, variant_alleles.tsv
"""
import json
import time
import urllib.request

import numpy as np
import pandas as pd
import pysam
import scanpy as sc
import MOODS.tools
from pyjaspar import jaspardb

from config import RAW, INTERIM, TAB

FLANK = 30
P_HIT = 1e-4
MIN_DELTA = 2.0
MIN_EXPR_FRAC = 0.10
ENS = RAW / "ensembl"
ENS.mkdir(exist_ok=True)
BASES = "ACGT"


def ensembl_alleles(rsids):
    path = ENS / "variants.json"
    if path.exists():
        return json.loads(path.read_text())
    out = {}
    for i in range(0, len(rsids), 100):
        req = urllib.request.Request("https://rest.ensembl.org/variation/homo_sapiens",
                                     data=json.dumps({"ids": rsids[i:i + 100]}).encode(),
                                     headers={"Content-Type": "application/json", "Accept": "application/json"})
        for attempt in range(5):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    out.update(json.load(r))
                break
            except Exception:
                if attempt == 4:
                    raise
                time.sleep(5 * (attempt + 1))
        time.sleep(1)
    path.write_text(json.dumps(out))
    return out


def motif_matrices():
    mats = []
    bg = MOODS.tools.flat_bg(4)
    for m in jaspardb(release="JASPAR2026").fetch_motifs(collection="CORE", tax_group=["vertebrates"]):
        counts = [[float(x) for x in m.counts[b]] for b in BASES]
        lo = np.array(MOODS.tools.log_odds(counts, bg, 0.8))           # 4 x L
        thr = MOODS.tools.threshold_from_p(lo.tolist(), bg, P_HIT)
        mats.append((m.name, m.matrix_id, lo, thr))
    return mats


def onehot(seq):
    idx = {b: i for i, b in enumerate(BASES)}
    x = np.zeros((4, len(seq)))
    for j, ch in enumerate(seq):
        if ch in idx:
            x[idx[ch], j] = 1
    return x


def best_score(lo, seq, centre):
    """Max log-odds over windows (both strands) that overlap position `centre` of `seq`."""
    L = lo.shape[1]
    x = onehot(seq)
    rc = lo[::-1, ::-1]                                               # reverse complement of the matrix
    best = -np.inf
    for start in range(max(0, centre - L + 1), min(centre, len(seq) - L) + 1):
        w = x[:, start:start + L]
        best = max(best, float((lo * w).sum()), float((rc * w).sum()))
    return best


def main():
    cand = pd.read_csv(TAB / "gwas_candidate_variants.tsv", sep="\t")
    cand = cand[cand["rs_id"].str.startswith("rs")]
    snps = cand.drop_duplicates("rs_id")[["rs_id", "location"]].copy()
    gw = pd.read_csv(RAW / "gwas" / "RA_associations.tsv", sep="\t")
    if "effect_allele" in gw.columns:
        risk = gw.dropna(subset=["effect_allele"]).sort_values("p_value").drop_duplicates("rs_id").set_index("rs_id")["effect_allele"]
    else:                                   # older download without risk alleles: orientation left unresolved
        print("WARNING: RA_associations.tsv has no effect_allele column; re-run step 09 to orient effects")
        risk = pd.Series(dtype=object)
    ens = ensembl_alleles(snps["rs_id"].tolist())
    fa = pysam.FastaFile(str(INTERIM / "hg38.fa"))

    rows = []
    for _, s in snps.iterrows():
        e = ens.get(s["rs_id"])
        maps = [m for m in (e or {}).get("mappings", []) if m.get("assembly_name") == "GRCh38" and "_" not in m["seq_region_name"]]
        if not maps or maps[0]["start"] != maps[0]["end"]:
            rows.append(dict(rs_id=s["rs_id"], status="not a mappable SNV")); continue
        m = maps[0]
        chrom, pos = "chr" + m["seq_region_name"], int(m["start"])
        alleles = m["allele_string"].split("/")
        ref_fa = fa.fetch(chrom, pos - 1, pos).upper()
        if ref_fa != alleles[0]:
            rows.append(dict(rs_id=s["rs_id"], status=f"reference mismatch ({alleles[0]} vs hg38 {ref_fa})")); continue
        ra = risk.get(s["rs_id"], None)
        rows.append(dict(rs_id=s["rs_id"], chrom=chrom, pos=pos, ref=alleles[0], alts=",".join(alleles[1:]),
                         risk_allele=ra, risk_is=("ref" if ra == alleles[0] else "alt" if ra in alleles[1:] else "unknown"),
                         status="ok"))
    al = pd.DataFrame(rows)
    al.to_csv(TAB / "variant_alleles.tsv", sep="\t", index=False)
    ok = al[al["status"] == "ok"]
    print(f"Candidate SNPs: {len(snps)}; usable SNVs: {len(ok)}; risk allele resolved for {(ok.risk_is != 'unknown').sum()}")

    # TF expression per state (fraction of nuclei detecting the TF gene)
    rna = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    mats = motif_matrices()
    tf_genes = sorted({g.upper() for n, _, _, _ in mats for g in n.split("::")} & set(rna.raw.var_names))
    det = sc.get.obs_df(rna, tf_genes, use_raw=True) > 0
    det["author_subtype"] = rna.obs["author_subtype"].astype(str).values
    det["author_ct"] = rna.obs["author_ct"].astype(str).values
    frac_sub = det.groupby("author_subtype")[tf_genes].mean()
    frac_ct = det.groupby("author_ct")[tf_genes].mean()

    def expressed(tf, state):
        genes = [g.upper() for g in tf.split("::")]
        tab = frac_sub if state in frac_sub.index else frac_ct if state in frac_ct.index else None
        if tab is None or not all(g in tab.columns for g in genes):
            return np.nan
        return float(min(tab.loc[state, g] for g in genes))

    states = cand.groupby("rs_id")["state"].apply(list)
    out = []
    for _, v in ok.iterrows():
        seq_ref = fa.fetch(v["chrom"], v["pos"] - 1 - FLANK, v["pos"] + FLANK).upper()
        for alt in v["alts"].split(","):
            if len(alt) != 1 or alt not in BASES:
                continue
            seq_alt = seq_ref[:FLANK] + alt + seq_ref[FLANK + 1:]
            for name, mid, lo, thr in mats:
                sr, sa = best_score(lo, seq_ref, FLANK), best_score(lo, seq_alt, FLANK)
                hit_r, hit_a = sr >= thr, sa >= thr
                if hit_r == hit_a or abs(sa - sr) < MIN_DELTA:
                    continue
                for st in states[v["rs_id"]]:
                    out.append(dict(rs_id=v["rs_id"], location=f"{v['chrom']}:{int(v['pos'])}", ref=v["ref"], alt=alt,
                                    risk_allele=v["risk_allele"], risk_is=v["risk_is"], state=st, tf=name, jaspar_id=mid,
                                    score_ref=round(sr, 2), score_alt=round(sa, 2), delta_alt_minus_ref=round(sa - sr, 2),
                                    effect="alt gains motif" if hit_a else "alt loses motif",
                                    tf_detected_frac_in_state=expressed(name, st)))
    res = pd.DataFrame(out)
    res["tf_expressed_in_state"] = res["tf_detected_frac_in_state"] >= MIN_EXPR_FRAC
    res["risk_allele_effect"] = np.where(res["risk_is"] == "alt", res["effect"].str.replace("alt ", "risk "),
                                np.where(res["risk_is"] == "ref", res["effect"].map({"alt gains motif": "risk loses motif",
                                                                                     "alt loses motif": "risk gains motif"}),
                                         "risk allele unresolved"))
    res = res.sort_values(["tf_expressed_in_state", "rs_id", "delta_alt_minus_ref"], ascending=[False, True, True])
    res.to_csv(TAB / "variant_motif_effects.tsv", sep="\t", index=False)
    keep = res[res["tf_expressed_in_state"]].drop_duplicates(["rs_id", "tf", "alt"])
    print(f"Motif-changing allele pairs: {res[['rs_id', 'tf', 'alt']].drop_duplicates().shape[0]} "
          f"across {res.rs_id.nunique()} SNPs; with the TF expressed in the peak's state: {len(keep)} across {keep.rs_id.nunique()} SNPs")
    pd.set_option("display.width", 250)
    print(keep.groupby("rs_id").apply(lambda g: pd.Series({
        "location": g.location.iloc[0], "states": "; ".join(sorted(set(res.loc[res.rs_id == g.name, "state"])))[:70],
        "risk_allele_effects": "; ".join(f"{t}: {e.replace('risk ', '')}" for t, e in zip(g.tf, g.risk_allele_effect))[:160]}),
        include_groups=False).to_string())


if __name__ == "__main__":
    main()
