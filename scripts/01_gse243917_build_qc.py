"""GSE243917: attach sample labels to raw counts, compute QC metrics, flag doublets.

Input : data/raw/GSE243917/*.h5ad.gz, series_matrix.txt.gz
Output: data/interim/gse243917_qc.h5ad
        data/processed/gse243917_sample_labels.tsv
        results/tables/gse243917_qc_per_sample.tsv
        results/figures/gse243917_qc.png
"""
import gzip
import re
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from config import RAW, INTERIM, PROCESSED, FIG, TAB, SEED

RAWDIR = RAW / "GSE243917"


def gunzip_once(name):
    src = RAWDIR / f"GSE243917_{name}.h5ad.gz"
    dst = INTERIM / f"GSE243917_{name}.h5ad"
    if not dst.exists():
        with gzip.open(src, "rb") as fi, open(dst, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    return dst


def grade(text):
    """Pull the numeric grade out of strings like 'Moderate (...) (2)'; OA/'Moderate' fall back by keyword."""
    m = re.search(r"\((\d)\)\s*$", text)
    if m:
        return int(m.group(1))
    if text.startswith("Mild"):
        return 1
    if text.startswith("Moderate"):
        return 2
    return np.nan


def sample_labels():
    rows = {}
    with gzip.open(RAWDIR / "GSE243917_series_matrix.txt.gz", "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title") or line.startswith("!Sample_characteristics_ch1"):
                key, *vals = line.rstrip("\n").split("\t")
                vals = [v.strip('"') for v in vals]
                if key == "!Sample_title":
                    rows["title"] = vals
                else:
                    field = vals[0].split(":")[0]
                    rows[field] = [v.split(": ", 1)[1] for v in vals]
    df = pd.DataFrame(rows)
    df = df[df["title"].str.endswith("_RNA")].copy()
    df["donor"] = "BRI-" + df["title"].str.extract(r"S(\d+)_")[0]
    df = df.rename(columns={"synovial lymphocytic_inflammation": "lymphocytic_inflammation_text",
                            "synovial lining_hyperplasia": "lining_hyperplasia_text"})
    df["lymphocytic_inflammation_grade"] = df["lymphocytic_inflammation_text"].map(grade)
    df["lining_hyperplasia_grade"] = df["lining_hyperplasia_text"].map(grade)
    # OA sample is graded "Mild (OA)" for both; keep as grade 1 but flag it
    return df[["donor", "disease", "lymphocytic_inflammation_grade", "lining_hyperplasia_grade",
               "lymphocytic_inflammation_text", "lining_hyperplasia_text"]].set_index("donor")


def main():
    labels = sample_labels()
    labels.to_csv(PROCESSED / "gse243917_sample_labels.tsv", sep="\t")
    print(labels[["disease", "lymphocytic_inflammation_grade", "lining_hyperplasia_grade"]])

    adata = sc.read_h5ad(gunzip_once("genesXcells_ALL_raw"))
    for n in ["B_RNA_processed", "T_RNA_processed", "endothelial_RNA_processed",
              "fibroblast_RNA_processed", "myeloid_RNA_processed"]:
        gunzip_once(n)
    adata.var_names_make_unique()
    adata.obs = adata.obs.rename(columns={"ct": "author_ct", "ct_subtype": "author_subtype"})
    adata.obs["donor"] = adata.obs_names.str.split("_").str[0]
    adata.obs = adata.obs.join(labels[["disease", "lymphocytic_inflammation_grade",
                                       "lining_hyperplasia_grade"]], on="donor")
    assert adata.obs["disease"].notna().all(), "unmatched donor IDs"

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.match(r"^RP[SL]\d")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

    # Authors already removed low-quality nuclei; we apply the same floor checks and report what would drop
    adata.obs["qc_pass"] = (adata.obs["n_genes_by_counts"] >= 200) & (adata.obs["pct_counts_mt"] < 20)

    # Doublet scores per sample (report-only; authors' doublet removal is trusted as primary)
    sc.pp.scrublet(adata, batch_key="donor", random_state=SEED)

    per_sample = adata.obs.groupby("donor").agg(
        n_cells=("donor", "size"),
        median_genes=("n_genes_by_counts", "median"),
        median_umis=("total_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
        pct_qc_fail=("qc_pass", lambda s: 100 * (1 - s.mean())),
        pct_scrublet_doublet=("predicted_doublet", lambda s: 100 * s.mean()),
    ).join(labels[["disease", "lymphocytic_inflammation_grade", "lining_hyperplasia_grade"]])
    per_sample.round(2).to_csv(TAB / "gse243917_qc_per_sample.tsv", sep="\t")
    print(per_sample.round(2))

    fig, axes = plt.subplots(1, 4, figsize=(18, 4))
    order = per_sample.index.tolist()
    for ax, col, title in zip(axes,
                              ["n_genes_by_counts", "total_counts", "pct_counts_mt", "doublet_score"],
                              ["Genes per nucleus", "UMIs per nucleus", "% mitochondrial", "Scrublet score"]):
        data = [adata.obs.loc[adata.obs["donor"] == d, col].values for d in order]
        ax.boxplot(data, showfliers=False)
        ax.set_xticks(range(1, len(order) + 1), [d.replace("BRI-", "") for d in order], rotation=90)
        ax.set_title(title)
        if col == "total_counts":
            ax.set_yscale("log")
    fig.suptitle("GSE243917 per-sample QC (sample 1281 = OA)")
    fig.tight_layout()
    fig.savefig(FIG / "gse243917_qc.png", dpi=150)

    adata.write_h5ad(INTERIM / "gse243917_qc.h5ad", compression="gzip")
    print(f"Saved {adata.n_obs} nuclei x {adata.n_vars} genes; "
          f"{(~adata.obs['qc_pass']).sum()} fail floor QC, {adata.obs['predicted_doublet'].sum()} scrublet-flagged")


if __name__ == "__main__":
    main()
