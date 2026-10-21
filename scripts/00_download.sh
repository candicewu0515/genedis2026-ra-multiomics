#!/usr/bin/env bash
# Download open RA synovium single-cell datasets and record checksums.
# Usage: bash scripts/00_download.sh   (run from the project root)
set -euo pipefail

cd "$(dirname "$0")/.."
PY="${PY:-python}"
"$PY" scripts/restore_references.py
RAW=data/raw
mkdir -p "$RAW/GSE243917" "$RAW/E-MTAB-8322"

# GSE243917: Gupta, Weinand et al. 2023 Nat Genet, single-nucleus multiome (RNA + ATAC) of RA/OA synovium
G=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243917
for f in genesXcells_ALL_raw peaksXcells_ALL_raw B_RNA_processed T_RNA_processed endothelial_RNA_processed \
         fibroblast_RNA_processed myeloid_RNA_processed; do
  [ -s "$RAW/GSE243917/GSE243917_${f}.h5ad.gz" ] || \
    curl -sS -L --retry 3 -o "$RAW/GSE243917/GSE243917_${f}.h5ad.gz" "$G/suppl/GSE243917_${f}.h5ad.gz"
done
[ -s "$RAW/GSE243917/GSE243917_series_matrix.txt.gz" ] || \
  curl -sS -L --retry 3 -o "$RAW/GSE243917/GSE243917_series_matrix.txt.gz" "$G/matrix/GSE243917_series_matrix.txt.gz"

# E-MTAB-8322: Alivernini et al. 2020, synovial tissue macrophages (EBI Single Cell Expression Atlas)
B=https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/E-MTAB-8322
for f in aggregated_counts.mtx.gz aggregated_counts.mtx_rows.gz aggregated_counts.mtx_cols.gz \
         aggregated_counts.decorated.mtx_rows cell_metadata.tsv condensed-sdrf.tsv idf.txt clusters.tsv; do
  [ -s "$RAW/E-MTAB-8322/E-MTAB-8322.$f" ] || \
    curl -sS -L --retry 3 -o "$RAW/E-MTAB-8322/E-MTAB-8322.$f" "$B/E-MTAB-8322.$f"
done

# GSE89408: Guo et al. 2017 J Immunol / Walsh et al. 2017 PLoS One: 218 synovial biopsies across disease stages
mkdir -p "$RAW/GSE89408"
G3=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89408
[ -s "$RAW/GSE89408/GSE89408_GEO_count_matrix_rename.txt.gz" ] || \
  curl -sS -L --retry 3 -o "$RAW/GSE89408/GSE89408_GEO_count_matrix_rename.txt.gz" "$G3/suppl/GSE89408_GEO_count_matrix_rename.txt.gz"
[ -s "$RAW/GSE89408/GSE89408_series_matrix.txt.gz" ] || \
  curl -sS -L --retry 3 -o "$RAW/GSE89408/GSE89408_series_matrix.txt.gz" "$G3/matrix/GSE89408_series_matrix.txt.gz"

# GSE198520: Wang et al. 2022 Arthritis Rheumatol, synovial bulk RNA-seq before/after anti-TNF (response labels)
mkdir -p "$RAW/GSE198520"
G2=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198520
[ -s "$RAW/GSE198520/GSE198520_Raw_gene_count_matrix.txt.gz" ] || \
  curl -sS -L --retry 3 -o "$RAW/GSE198520/GSE198520_Raw_gene_count_matrix.txt.gz" "$G2/suppl/GSE198520_Raw_gene_count_matrix.txt.gz"
[ -s "$RAW/GSE198520/GSE198520_series_matrix.txt.gz" ] || \
  curl -sS -L --retry 3 -o "$RAW/GSE198520/GSE198520_series_matrix.txt.gz" "$G2/matrix/GSE198520_series_matrix.txt.gz"

# Ensembl 98 gene annotation (same as 10x GRCh38-2020-A) for ATAC gene activity
mkdir -p "$RAW/annotation"
[ -s "$RAW/annotation/Homo_sapiens.GRCh38.98.gtf.gz" ] || \
  curl -sS -L --retry 3 -o "$RAW/annotation/Homo_sapiens.GRCh38.98.gtf.gz" \
    "https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz"

# Blood (CELLxGENE, step 18) and TCGA lung (UCSC Xena, step 19)
mkdir -p "$RAW/cellxgene_blood" "$RAW/tcga"
[ -s "$RAW/cellxgene_blood/RA_disease_activity_blood.h5ad" ] || \
  curl -sS -L --retry 3 -o "$RAW/cellxgene_blood/RA_disease_activity_blood.h5ad" \
    "https://datasets.cellxgene.cziscience.com/3b751975-34bb-409a-a9b7-98380f0450ea.h5ad"
X=https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download
for c in LUAD LUSC; do
  [ -s "$RAW/tcga/${c}_HiSeqV2.gz" ] || curl -sS -L --retry 3 -o "$RAW/tcga/${c}_HiSeqV2.gz" "$X/TCGA.${c}.sampleMap%2FHiSeqV2.gz"
  [ -s "$RAW/tcga/${c}_clinicalMatrix.tsv" ] || \
    curl -sS -L --retry 3 -o "$RAW/tcga/${c}_clinicalMatrix.tsv" "$X/TCGA.${c}.sampleMap%2F${c}_clinicalMatrix"
done
[ -s "$RAW/tcga/Survival_SupplementalTable_S1_20171025_xena_sp.tsv" ] || \
  curl -sS -L --retry 3 -o "$RAW/tcga/Survival_SupplementalTable_S1_20171025_xena_sp.tsv" \
    "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp"

# Xenium per-sample matrices and cell tables (GSE347958; step 20); file list from GEO filelist.txt
mkdir -p "$RAW/GSE347958" "$RAW/PXD062696"
if [ ! -s "$RAW/GSE347958/selected_files.tsv" ]; then
  curl -sS -L "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE347nnn/GSE347958/suppl/filelist.txt" | \
    awk -F'\t' 'NR>1 && $1=="File" && ($2 ~ /(matrix\.mtx|features\.tsv|barcodes\.tsv|cells\.parquet)\.gz$/){print $2"\t"$4}' \
    > "$RAW/GSE347958/selected_files.tsv"
fi
while IFS=$'\t' read -r f size; do
  gsm=${f%%_*}
  [ -s "$RAW/GSE347958/$f" ] || curl -sS -L --retry 3 -o "$RAW/GSE347958/$f" \
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/${gsm:0:${#gsm}-3}nnn/$gsm/suppl/$f"
done < "$RAW/GSE347958/selected_files.tsv"
# Visium raw count matrices and positions (step 24)
mkdir -p "$RAW/GSE347797"
[ -s "$RAW/GSE347797/GSE347797_RAW.tar" ] || curl -f -sS -L --retry 3 -o "$RAW/GSE347797/GSE347797_RAW.tar" "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE347nnn/GSE347797/suppl/GSE347797_RAW.tar"
# LCM proteomics, DIA-NN processed matrices (PRIDE PXD062696; step 21)
for f in experiment_design.xlsx report-first-pass.gg_matrix.tsv report-first-pass.pg_matrix.tsv report-first-pass.stats.tsv; do
  [ -s "$RAW/PXD062696/$f" ] || curl -sS -L --retry 3 -o "$RAW/PXD062696/$f" \
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/06/PXD062696/$f"
done

# UCSC hg38 genome sequence (chr-prefixed, matches 10x GRCh38 peak names) for TF motif analysis (step 14)
mkdir -p "$RAW/genome"
[ -s "$RAW/genome/hg38.fa.gz" ] || \
  curl -sS -L --retry 3 -o "$RAW/genome/hg38.fa.gz" "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
[ -s data/interim/hg38.fa ] || { mkdir -p data/interim && gunzip -c "$RAW/genome/hg38.fa.gz" > data/interim/hg38.fa; }

# S-LDSC reference (Broad, Zenodo 10515792: GRCh38 1000G EUR plink files, baselineLD v2.2, HapMap3 weights; ~1.4 GB,
# 3.9 GB extracted) and full GWAS summary statistics for step 33
mkdir -p "$RAW/ldsc" "$RAW/gwas_sumstats"
Z=https://zenodo.org/records/10515792/files
[ -s "$RAW/ldsc/GRCh38.tgz" ] || curl -sS -L --retry 3 -o "$RAW/ldsc/GRCh38.tgz" "$Z/GRCh38.tgz?download=1"
[ -s "$RAW/ldsc/hm3_no_MHC.list.txt" ] || curl -sS -L --retry 3 -o "$RAW/ldsc/hm3_no_MHC.list.txt" "$Z/hm3_no_MHC.list.txt?download=1"
if [ ! -d "$RAW/ldsc/GRCh38/baselineLD_v2.2" ]; then
  (cd "$RAW/ldsc" && tar xzf GRCh38.tgz && cd GRCh38 && for f in baselineLD_v2.2 plink_files weights; do tar xzf $f.tgz && rm $f.tgz; done && rm -f baseline_v1.2.tgz)
fi
[ -s "$RAW/gwas_sumstats/GCST90132223.h.tsv.gz" ] || curl -sS -L --retry 3 -o "$RAW/gwas_sumstats/GCST90132223.h.tsv.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132223/harmonised/GCST90132223.h.tsv.gz"
for e in M13_RHEUMA T2D; do
  [ -s "$RAW/gwas_sumstats/finngen_R12_$e.gz" ] || curl -sS -L --retry 3 -o "$RAW/gwas_sumstats/finngen_R12_$e.gz" \
    "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_$e.gz"
done
# Python 3 LDSC (polyfun), pinned commit; run in the environment of environment_ldsc.yml
TOOLS="$(cd "$(dirname "$0")/.." && pwd)/tools"
if [ ! -d "$TOOLS/polyfun" ]; then
  mkdir -p "$TOOLS" && git clone -q https://github.com/omerwe/polyfun.git "$TOOLS/polyfun" && \
    git -C "$TOOLS/polyfun" checkout -q 8b17b4ec87fde11f8d216c049c08d1420cd25de6
fi

# MAGMA v1.10 (macOS binary; use magma_v1.10_static.zip on Linux) and NCBI GRCh38 gene locations for step 36
if [ ! -x "$TOOLS/magma/magma" ]; then
  mkdir -p "$TOOLS/magma" && curl -sS -L -o "$TOOLS/magma/magma.zip" "https://vu.data.surfsara.nl/index.php/s/1M1d9vHtVidEwvU/download" && \
    (cd "$TOOLS/magma" && unzip -o -q magma.zip && rm magma.zip)
fi
mkdir -p "$RAW/magma"
[ -s "$RAW/magma/NCBI38.gene.loc" ] || (curl -sS -L -o "$RAW/magma/NCBI38.zip" "https://vu.data.surfsara.nl/index.php/s/yj952iHqy5anYhH/download" && \
  cd "$RAW/magma" && unzip -o -q NCBI38.zip && rm NCBI38.zip)

# Second RA single-cell cohort for step 40: GSE296117 synovial fluid before/after adalimumab or tofacitinib (He C et al.
# 2025 Sci Data). Deposited only as a Seurat object (2.5 GB); converted once with R + Seurat (scripts/convert_gse296117.R).
mkdir -p "$RAW/gse296117"
[ -s "$RAW/gse296117/GSE296117_RA_geo.rds" ] || curl -sS -L --retry 5 -C - -o "$RAW/gse296117/GSE296117_RA_geo.rds" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296117/suppl/GSE296117_RA_geo.rds"
[ -s "$RAW/gse296117/converted/meta.tsv.gz" ] || \
  Rscript scripts/convert_gse296117.R "$RAW/gse296117/GSE296117_RA_geo.rds" "$RAW/gse296117/converted"

# NSCLC single-cell for step 41: LuCA core atlas (Salcher et al. 2022 Cancer Cell; CELLxGENE, 12.9 GB h5ad)
mkdir -p "$RAW/luca"
[ -s "$RAW/luca/luca_core_atlas.h5ad" ] || curl -sS -L --retry 5 -C - -o "$RAW/luca/luca_core_atlas.h5ad" \
  "https://datasets.cellxgene.cziscience.com/46e0287b-9a33-4e83-99f3-8c044131bfdc.h5ad"

"$PY" - <<'CHECK'
from pathlib import Path
import hashlib,json,subprocess
spec=json.loads(Path('env/external_tools.json').read_text())
commit=subprocess.check_output(['git','-C','tools/polyfun','rev-parse','HEAD'],text=True).strip()
assert commit==spec['polyfun_commit'],'Polyfun commit differs from the pinned version'
assert hashlib.sha256(Path('tools/magma/magma').read_bytes()).hexdigest()==spec['magma_sha256'],'MAGMA binary differs from the tested version/platform'
print('PASS: pinned LDSC and MAGMA tools')
CHECK
"$PY" scripts/verify_raw_inputs.py

echo "Downloads complete; checksums in $RAW/SHA256SUMS.txt"
echo "Fetched by their own steps: GWAS Catalog (09), MSigDB (12), JASPAR (14), Open Targets (15), Ensembl variants (16)."
