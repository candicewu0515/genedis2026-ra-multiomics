"""eQTL colocalisation: do RA risk loci and myeloid eQTLs share a causal variant at the prioritised genes?

GWAS: Ishigaki et al. 2022 RA, European ancestry (GWAS Catalog GCST90132223, harmonised GRCh38; read from the local
copy fetched by 00_download.sh when present, else by remote tabix).
eQTL Catalogue (tabix FTP; the REST API is no longer available), gene expression ('ge') datasets:
  QTD000001 Alasoo_2018 macrophage naive (n = 84); QTD000006 Alasoo_2018 macrophage IFN-gamma 18 h (n = 84);
  QTD000384 Nedelec_2016 macrophage naive (n = 163); QTD000409 Quach_2016 monocyte naive (n = 200).
Genes: prioritisation genes with genetic evidence (step 22, L5 > 0); region = TSS +/- 500 kb; MHC excluded.
Method: coloc.abf (Giambartolomei 2014) with Wakefield approximate Bayes factors (prior SD 0.2 GWAS log-OR, 0.15 eQTL),
priors p1 = p2 = 1e-4, p12 = 1e-5; variants matched on GRCh38 position and allele pair; >= 50 shared variants.
PP.H4 > 0.8 = shared causal variant; PP.H3 > 0.8 = distinct variants. Single-causal-variant assumption.

Output: data/raw/eqtl/tabix_ftp_paths.tsv (+ cached .tbi indexes and regions/: every remote region read, so reruns
        work offline), results/tables/eqtl_coloc.tsv, results/figures/eqtl_coloc.png
"""
import gzip
import os
import re
import time
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from scipy.special import logsumexp

from config import RAW, FIG, TAB

GWAS = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132223/"
        "harmonised/GCST90132223.h.tsv.gz")
DATASETS = ["QTD000001", "QTD000006", "QTD000384", "QTD000409"]
WINDOW = 500_000
MHC = ("6", 28_510_120, 33_480_577)
P1, P2, P12 = 1e-4, 1e-4, 1e-5
E = RAW / "eqtl"
E.mkdir(exist_ok=True)


def header(url, name):
    """Column names of a remote bgzipped table (cached, retried)."""
    cache = E / "regions" / f"{name}.header.txt"
    if cache.exists():
        return cache.read_text().rstrip("\n").split("\t")
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                with gzip.GzipFile(fileobj=r) as g:
                    h = g.readline().decode().rstrip("\n")
            cache.parent.mkdir(exist_ok=True)
            cache.write_text(h + "\n")
            return h.split("\t")
        except OSError:
            if attempt == 4:
                raise
            time.sleep(15 * (attempt + 1))


class Tabix:
    """Remote/local tabix reader that re-opens and retries on transient read errors (HTTP range requests can fail).
    A contig absent from the file returns None (skipped, as before); any other persistent error stops the step."""
    def __init__(self, path, index=None, name=None):
        self.path, self.index, self.name, self.t = path, index, name, None

    def _open(self):
        here = os.getcwd()
        os.chdir(E)                    # pysam caches remote .tbi indexes in the working directory
        try:
            self.t = pysam.TabixFile(self.path, index=self.index) if self.index else pysam.TabixFile(self.path)
        finally:
            os.chdir(here)

    def fetch(self, chrom, lo, hi):
        cache = E / "regions" / f"{self.name}_{chrom}_{lo}_{hi}.tsv.gz" if self.name else None
        if cache is not None and cache.exists():           # remote regions are cached so later runs work offline
            with gzip.open(cache, "rt") as fh:
                first = fh.readline().rstrip("\n")
                if first == "#absent_contig":
                    return None
                return ([first] if first else []) + [l.rstrip("\n") for l in fh]
        for attempt in range(8):                            # EBI refuses connections for minutes at a time
            try:
                if self.t is None:
                    self._open()
                contig = f"chr{chrom}" if self.t.contigs and self.t.contigs[0].startswith("chr") else chrom
                rows = None if contig not in self.t.contigs else list(self.t.fetch(contig, lo, hi))
                time.sleep(1)                               # be gentle with the server
                break
            except (ValueError, OSError):
                if attempt == 7:
                    raise
                time.sleep(min(30 * (attempt + 1), 180))
                self.t = None
        if cache is not None:
            cache.parent.mkdir(exist_ok=True)
            with gzip.open(cache, "wt") as fh:
                fh.write("#absent_contig\n" if rows is None else "".join(r + "\n" for r in rows))
        return rows


def labf(beta, se, w):
    v = se ** 2
    r = w ** 2 / (w ** 2 + v)
    z = beta / se
    return 0.5 * np.log(1 - r) + 0.5 * r * z ** 2


def coloc_abf(l1, l2):
    s1, s2, s12 = logsumexp(l1), logsumexp(l2), logsumexp(l1 + l2)
    h = np.array([0.0, np.log(P1) + s1, np.log(P2) + s2,
                  np.log(P1) + np.log(P2) + (s1 + s2 + np.log1p(-np.exp(min(s12 - (s1 + s2), -1e-12)))),
                  np.log(P12) + s12])
    return np.exp(h - logsumexp(h))


def genes_tss():
    rows = []
    with gzip.open(RAW / "annotation" / "Homo_sapiens.GRCh38.98.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "gene":
                continue
            nm = re.search(r'gene_name "([^"]+)"', f[8]); gid = re.search(r'gene_id "([^"]+)"', f[8])
            if nm:
                rows.append((nm.group(1), gid.group(1), f[0], int(f[3]) if f[6] == "+" else int(f[4])))
    return pd.DataFrame(rows, columns=["gene", "ensembl", "chrom", "tss"]).drop_duplicates("gene").set_index("gene")


def main():
    paths = E / "tabix_ftp_paths.tsv"
    if not paths.exists():
        urllib.request.urlretrieve("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/"
                                   "tabix/tabix_ftp_paths.tsv", paths)
    tp = pd.read_csv(paths, sep="\t").set_index("dataset_id")
    pr = pd.read_csv(TAB / "gene_prioritisation.tsv", sep="\t", index_col=0)
    genes = pr.index[pr["L5_genetics"] > 0].tolist()
    tss = genes_tss()
    genes = [g for g in genes if g in tss.index and not (tss.loc[g, "chrom"] == MHC[0] and MHC[1] <= tss.loc[g, "tss"] <= MHC[2])]
    print(f"Genes with genetic evidence for colocalisation: {len(genes)}")

    local = RAW / "gwas_sumstats" / "GCST90132223.h.tsv.gz"       # same file, downloaded for step 33
    if local.exists() and (E / "GCST90132223.h.tsv.gz.tbi").exists():
        with gzip.open(local, "rt") as fh:
            gh = fh.readline().rstrip("\n").split("\t")
        gw_t = Tabix(str(local), index=str(E / "GCST90132223.h.tsv.gz.tbi"))
    else:
        gh = header(GWAS, "GCST90132223")
        gw_t = Tabix(GWAS, name="GCST90132223")
    eq = {}
    for ds in DATASETS:
        url = tp.loc[ds, "ftp_path"].replace("ftp://", "https://")
        eq[ds] = (Tabix(url, name=ds), header(url, ds), tp.loc[ds, "study_label"] + " " + tp.loc[ds, "sample_group"])
    print("GWAS columns:", gh[:9])

    rows = []
    for g in genes:
        chrom, pos, ens = str(tss.loc[g, "chrom"]), int(tss.loc[g, "tss"]), tss.loc[g, "ensembl"]
        lo, hi = max(1, pos - WINDOW), pos + WINDOW
        gw = pd.DataFrame([r.split("\t") for r in (gw_t.fetch(chrom, lo, hi) or [])], columns=gh)
        if gw.empty:
            continue
        gw = gw.rename(columns={"base_pair_location": "pos", "standard_error": "se"})
        gw = gw[(gw.beta != "NA") & (gw.se != "NA")]
        gw["pos"] = gw.pos.astype(int); gw["beta"] = gw.beta.astype(float); gw["se"] = gw.se.astype(float)
        gw["p_value"] = pd.to_numeric(gw.p_value, errors="coerce")
        gw["key"] = gw.pos.astype(str) + ":" + [":".join(sorted(x)) for x in zip(gw.effect_allele, gw.other_allele)]
        gw = gw.drop_duplicates("key")
        for ds, (t, h, label) in eq.items():
            rows_e = t.fetch(chrom, lo, hi)
            if rows_e is None:
                continue
            e = pd.DataFrame([r.split("\t") for r in rows_e], columns=h)
            e = e[e.gene_id == ens]
            if e.empty:
                rows.append(dict(gene=g, dataset=ds, label=label, n_shared=0)); continue
            e = e[(e.beta != "NA") & (e.se != "NA")]
            e["pos"] = e.position.astype(int); e["beta"] = e.beta.astype(float); e["se"] = e.se.astype(float)
            e["pvalue"] = pd.to_numeric(e.pvalue, errors="coerce")
            e["key"] = e.pos.astype(str) + ":" + [":".join(sorted(x)) for x in zip(e.ref, e.alt)]
            e = e.drop_duplicates("key")
            if "rsid" not in e.columns:                      # older dataset files lack rsid; keep the variant id
                e["rsid"] = e["variant"] if "variant" in e.columns else ""
            m = gw[["key", "pos", "beta", "se", "p_value"]].merge(
                e[["key", "beta", "se", "pvalue", "rsid"]].rename(columns={"rsid": "variant_label"}), on="key", suffixes=("_gwas", "_eqtl"))
            m = m[(m.se_gwas > 0) & (m.se_eqtl > 0)]
            if len(m) < 50:
                rows.append(dict(gene=g, dataset=ds, label=label, n_shared=len(m))); continue
            l1 = labf(m.beta_gwas.values, m.se_gwas.values, 0.2)
            l2 = labf(m.beta_eqtl.values, m.se_eqtl.values, 0.15)
            pp = coloc_abf(l1, l2)
            lead = m.iloc[int(np.argmax(l1 + l2))]
            rows.append(dict(gene=g, dataset=ds, label=label, n_shared=len(m), gwas_min_p=float(m.p_value.min()),
                             eqtl_min_p=float(m.pvalue.min()), PP_H0=pp[0], PP_H1=pp[1], PP_H2=pp[2], PP_H3=pp[3], PP_H4=pp[4],
                             top_shared_variant=lead.variant_label, top_variant_pos=f"{chrom}:{int(lead.pos)}"))
        best = max((r.get("PP_H4", 0) for r in rows if r["gene"] == g), default=0)
        print(f"  {g}: done (best PP.H4 {best:.2f})", flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(TAB / "eqtl_coloc.tsv", sep="\t", index=False, float_format="%.4g")
    tested = res.dropna(subset=["PP_H4"])
    print(f"\nTests run: {len(tested)} (gene x dataset with >= 50 shared variants)")
    hits = tested[tested.PP_H4 > 0.8].sort_values("PP_H4", ascending=False)
    print("PP.H4 > 0.8 (shared causal variant):")
    print(hits[["gene", "label", "n_shared", "gwas_min_p", "eqtl_min_p", "PP_H4", "top_shared_variant"]].round(4).to_string(index=False)
          if len(hits) else "  none")
    print("PP.H3 > 0.8 (distinct variants):", ", ".join(f"{r.gene} ({r.label})" for r in tested[tested.PP_H3 > 0.8].itertuples()) or "none")
    print("Suggestive 0.5 < PP.H4 <= 0.8:", ", ".join(f"{r.gene} ({r.label}, {r.PP_H4:.2f})" for r in tested[(tested.PP_H4 > 0.5) & (tested.PP_H4 <= 0.8)].itertuples()) or "none")

    if len(tested):
        piv = tested.pivot_table(index="gene", columns="label", values="PP_H4")
        piv = piv.loc[piv.max(axis=1).sort_values(ascending=False).index]
        fig, ax = plt.subplots(figsize=(8, max(4, 0.3 * len(piv))))
        im = ax.imshow(piv.values, aspect="auto", cmap="Reds", vmin=0, vmax=1)
        ax.set_xticks(range(piv.shape[1]), piv.columns, rotation=30, ha="right", fontsize=8)
        ax.set_yticks(range(piv.shape[0]), piv.index, fontsize=7)
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                v = piv.values[i, j]
                if np.isfinite(v) and v > 0.5:
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6)
        ax.set_title("RA GWAS (EUR) x myeloid eQTL colocalisation, PP.H4", fontsize=10)
        plt.colorbar(im, ax=ax, shrink=0.5)
        fig.tight_layout(); fig.savefig(FIG / "eqtl_coloc.png", dpi=150)


if __name__ == "__main__":
    main()
