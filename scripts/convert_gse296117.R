# Convert the GSE296117 Seurat object (RA synovial fluid, pre/post adalimumab or tofacitinib) to plain files for Python.
# Called by scripts/00_download.sh; needs R with the Seurat and Matrix packages.
# Writes data/raw/gse296117/converted/{counts.mtx.gz, genes.tsv.gz, barcodes.tsv.gz, meta.tsv.gz}
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
src <- args[1]; out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
obj <- readRDS(src)
obj <- UpdateSeuratObject(obj)
assay <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)
cnt <- tryCatch(LayerData(obj, assay = assay, layer = "counts"),
                error = function(e) GetAssayData(obj, assay = assay, slot = "counts"))
if (is.list(cnt)) cnt <- do.call(cbind, cnt)          # split layers (Seurat v5)
cnt <- as(cnt, "CsparseMatrix")
meta <- obj@meta.data[colnames(cnt), , drop = FALSE]
writeMM(cnt, file.path(out, "counts.mtx")); system2("gzip", c("-f", file.path(out, "counts.mtx")))
g <- gzfile(file.path(out, "genes.tsv.gz"), "w"); writeLines(rownames(cnt), g); close(g)
b <- gzfile(file.path(out, "barcodes.tsv.gz"), "w"); writeLines(colnames(cnt), b); close(b)
m <- gzfile(file.path(out, "meta.tsv.gz"), "w")
write.table(cbind(barcode = rownames(meta), meta), m, sep = "\t", quote = FALSE, row.names = FALSE); close(m)
cat(sprintf("converted: %d genes x %d cells; metadata columns: %s\n", nrow(cnt), ncol(cnt), paste(colnames(meta), collapse = ", ")))
