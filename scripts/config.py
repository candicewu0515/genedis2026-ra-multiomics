"""Shared paths and constants for the RA synovium pipeline."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROCESSED = ROOT / "data" / "processed"
# Per-step diagnostic figures; GENEDIS_FIG_DIR redirects them (verify_quick.sh writes them to a temporary folder)
FIG = Path(os.environ["GENEDIS_FIG_DIR"]) if os.environ.get("GENEDIS_FIG_DIR") else ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"

for p in (INTERIM, PROCESSED, FIG, TAB):
    p.mkdir(parents=True, exist_ok=True)

SEED = 0

# Broad synovial cell-type markers (nuclear RNA friendly; avoid ribosomal/mito genes)
BROAD_MARKERS = {
    "fibroblast": ["PDGFRA", "COL1A1", "COL3A1", "PRG4", "FAP", "DCN"],
    "mural": ["RGS5", "MYH11", "ACTA2", "NOTCH3", "PDGFRB", "KCNJ8"],
    # PECAM1 left out on purpose: monocytes express it and were mis-called endothelial
    "endothelial": ["VWF", "CDH5", "EMCN", "PTPRB", "FLT1", "PLVAP"],
    # covers both tissue macrophages (C1QA, CSF1R) and infiltrating monocytes (CD14, LYZ)
    "myeloid": ["CD14", "LYZ", "CSF1R", "C1QA", "TYROBP", "FCER1G", "AIF1"],
    "Tcell": ["CD3E", "CD3D", "CD2", "CD247", "IL7R", "THEMIS"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1", "NCAM1"],
    "Bcell": ["MS4A1", "BANK1", "CD79A", "PAX5", "BLK", "FCRL1"],
    "plasma": ["MZB1", "JCHAIN", "XBP1", "PRDM1", "TNFRSF17", "IGHG1"],
}

# Macrophage programs used across both datasets (Alivernini 2020 Nat Med; Zhang 2023 Nature)
MAC_HOMEOSTATIC = ["MERTK", "LYVE1", "FOLR2", "SELENOP", "TREM2", "F13A1", "MRC1", "CD163", "STAB1", "TIMD4"]
MAC_INFLAMMATORY = ["S100A8", "S100A9", "S100A12", "SPP1", "IL1B", "FCN1", "CD48", "ICAM1", "CXCL10", "CCL20"]
