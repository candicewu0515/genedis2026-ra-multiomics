"""Integrated, transparent gene prioritisation across all layers of the workflow.

Universe: the 115 workflow-derived genes of step 15 (signature, GWAS, communication, pathway sources) plus controls.
Each layer contributes 1 (strong), 0.5 (supporting) or 0, using only outputs of earlier steps:
  L1 scRNA disease     E-MTAB-8322 active vs remission macrophages: FDR < 0.05 & up = 1; p < 0.05 & up = 0.5
  L2 tissue state      GSE243917: log2 enrichment in M-6 STAT1+ CXCL10+ vs other macrophages > 1 and detected in >= 10 % of M-6 nuclei
  L3 treatment         GSE198520 good responders week 12 vs baseline: FDR < 0.05 & down = 1; p < 0.05 & down = 0.5
  L4 chromatin         significant peak-gene link = 1; RNA-gene-activity r > 0.5 = 0.5
  L5 genetics          nearest gene of an RA lead SNP inside a state-specific peak = 1; Open Targets genetic score >= 0.5 = 0.5
  L6 variant motif     that SNP changes a motif of a TF expressed in the peak's state = 1
  L7 pathway           IFN-gamma response leading edge in both GSEA contrasts = 1
  L8 communication     ligand/receptor of an M-6 interaction (IFN-gamma, CXCR3 axes or top-20 M-6-specific) = 1
  L9 protein           LCM proteomics RA > NH (lining or sublining donor means, Mann-Whitney p < 0.05) = 1
  L10 colocalisation   RA GWAS x myeloid eQTL (step 25): PP.H4 > 0.8 in any dataset = 1; 0.5-0.8 = 0.5; or a fine-mapped
                       RA signal colocalising (H4 >= 0.8, cis) with an immune-cell QTL of the gene (step 32) = 1
  L11 fine-mapping     nearest gene of a fine-mapped variant (PIP >= 0.5 = 1; 0.1-0.5 = 0.5) inside a state-specific
                       peak, in any RA fine-mapping of step 29 (Ishigaki 2022 EUR or EAS, FinnGen R12)
evidence_score = sum of L1-L11 (max 11). Druggability and novelty (Open Targets) are reported alongside, not summed.

Output: results/tables/gene_prioritisation.tsv; results/figures/gene_prioritisation.png
"""
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu

from config import RAW, INTERIM, FIG, TAB

M6 = "M-6: STAT1+ CXCL10+"
DC_LIKE = ["M-9", "M-10", "M-11", "M-12", "M-13", "M-14"]
LAYERS = ["L1_scRNA_disease", "L2_tissue_M6", "L3_treatment", "L4_chromatin", "L5_genetics", "L6_variant_motif",
          "L7_pathway", "L8_communication", "L9_protein", "L10_coloc", "L11_finemap"]


def main():
    mine = pd.read_csv(TAB / "target_drug_mining.tsv", sep="\t").set_index("gene")
    genes = mine.index.tolist()
    ev = pd.DataFrame(0.0, index=genes, columns=LAYERS)

    de = pd.read_csv(TAB / "emtab8322_DE_active_vs_remission.tsv", sep="\t", index_col=0)
    for g in genes:
        if g in de.index:
            r = de.loc[g]
            ev.loc[g, "L1_scRNA_disease"] = 1.0 if (r.padj < 0.05 and r.log2FoldChange > 0) else 0.5 if (r.pvalue < 0.05 and r.log2FoldChange > 0) else 0

    rna = sc.read_h5ad(INTERIM / "gse243917_annotated.h5ad")
    code = rna.obs["author_subtype"].astype(str).str.split(":").str[0]
    mac = rna[(rna.obs["author_ct"] == "myeloid") & ~code.isin(DC_LIKE)]
    present = [g for g in genes if g in mac.raw.var_names]
    df = sc.get.obs_df(mac, present, use_raw=True)
    is6 = (mac.obs["author_subtype"] == M6).values
    m6, oth = df[is6], df[~is6]
    lfc = np.log2((np.expm1(m6).mean() + 1e-3) / (np.expm1(oth).mean() + 1e-3))
    det = (m6 > 0).mean()
    for g in present:
        ev.loc[g, "L2_tissue_M6"] = 1.0 if (lfc[g] > 1 and det[g] >= 0.10) else 0

    tr = pd.read_csv(TAB / "gse198520_DE_post_vs_pre_good.tsv", sep="\t", index_col=0)
    for g in genes:
        if g in tr.index:
            r = tr.loc[g]
            ev.loc[g, "L3_treatment"] = 1.0 if (r.padj < 0.05 and r.log2FoldChange < 0) else 0.5 if (r.pvalue < 0.05 and r.log2FoldChange < 0) else 0

    links = pd.read_csv(TAB / "peak_gene_links.tsv.gz", sep="\t")
    linked = set(links[(links.r > 0.45) & (links.fdr < 0.05)].gene)
    conc = pd.read_csv(TAB / "rna_atac_concordance.tsv", sep="\t").set_index("gene")["r_rna_vs_gene_activity"]
    for g in genes:
        ev.loc[g, "L4_chromatin"] = 1.0 if g in linked else 0.5 if (g in conc.index and conc[g] > 0.5) else 0

    cv = pd.read_csv(TAB / "gwas_candidate_variants.tsv", sep="\t")
    gwas_genes = set(cv.nearest_protein_coding_gene)
    for g in genes:
        ev.loc[g, "L5_genetics"] = 1.0 if g in gwas_genes else 0.5 if mine.loc[g, "ot_ra_genetic"] >= 0.5 else 0

    vm = pd.read_csv(TAB / "variant_motif_effects.tsv", sep="\t")
    vm = vm[vm.tf_expressed_in_state & ~vm.tf.str.upper().str.match(r"^(ZNF|ZSCAN|ZKSCAN|ZBTB|ZBED|ZFP)")]
    snp_gene = cv.drop_duplicates("rs_id").set_index("rs_id")["nearest_protein_coding_gene"]
    vm_genes = set(snp_gene.reindex(vm.rs_id.unique()).dropna())
    for g in genes:
        ev.loc[g, "L6_variant_motif"] = 1.0 if g in vm_genes else 0

    gs = pd.read_csv(TAB / "gsea_results.tsv", sep="\t")
    le = [set(str(x).split(";")) for x in gs[(gs.pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE") & (gs.library == "hallmark")]["Lead_genes"]]
    le_both = set.intersection(*le) if len(le) == 2 else set()
    for g in genes:
        ev.loc[g, "L7_pathway"] = 1.0 if g in le_both else 0

    comm = set(mine[mine.sources.str.contains("communication")].index)
    for g in genes:
        ev.loc[g, "L8_communication"] = 1.0 if g in comm else 0

    # protein: per-gene RA vs NH, donor means, lining or sublining
    P = RAW / "PXD062696"
    design = pd.read_excel(P / "experiment_design.xlsx").set_index("Proteomics Sample ID")
    gg = pd.read_csv(P / "report-first-pass.gg_matrix.tsv", sep="\t")
    gg = gg[gg["Genes"].notna() & (gg["Genes"] != "-")]
    gg["gene"] = gg["Genes"].str.split(";").str[0]
    gg = gg.drop(columns="Genes").groupby("gene").max()
    gg.columns = [re.search(r"(P\d+-[LSI])_", c).group(1) for c in gg.columns]
    meta = design.loc[gg.columns]
    use = meta.index[meta["Disease"].isin(["NH", "RA"])]
    X = np.log2(gg[use]); X = X - X.median(axis=0)
    meta = meta.loc[use].assign(lining=lambda m: m["Tissue description"].str.contains("Lining"),
                                donor=lambda m: m["CT#"].astype(str) + "_" + m["LabLion Block ID"].astype(str))
    prot_p = {}
    for g in [g for g in genes if g in X.index]:
        best = 1.0
        for lin in (True, False):
            sel = meta[meta.lining == lin]
            v = X.loc[g, sel.index].groupby(sel["donor"]).mean()
            ra = v[sel.groupby("donor")["Disease"].first() == "RA"].dropna()
            nh = v[sel.groupby("donor")["Disease"].first() == "NH"].dropna()
            if len(ra) >= 3 and len(nh) >= 2 and ra.median() > nh.median():
                best = min(best, mannwhitneyu(ra, nh).pvalue)
        prot_p[g] = best
        ev.loc[g, "L9_protein"] = 1.0 if best < 0.05 else 0

    cpath = TAB / "eqtl_coloc.tsv"
    if cpath.exists():
        co = pd.read_csv(cpath, sep="\t").dropna(subset=["PP_H4"]).groupby("gene")["PP_H4"].max()
        for g in genes:
            if g in co.index:
                ev.loc[g, "L10_coloc"] = 1.0 if co[g] > 0.8 else 0.5 if co[g] > 0.5 else 0
    mpath = TAB / "coloc_multisignal_genes.tsv"
    if mpath.exists():
        ms = pd.read_csv(mpath, sep="\t")
        imm = set(ms.loc[ms.immune_records > 0, "gene"])
        for g in genes:
            if g in imm:
                ev.loc[g, "L10_coloc"] = 1.0
    fpath = TAB / "finemap_candidates.tsv"
    if fpath.exists():
        fm = pd.read_csv(fpath, sep="\t").groupby("nearest_protein_coding_gene")["pip"].max()
        for g in genes:
            if g in fm.index:
                ev.loc[g, "L11_finemap"] = 1.0 if fm[g] >= 0.5 else 0.5
    out = ev.copy()
    out["evidence_score"] = ev.sum(axis=1)
    out["n_layers_with_evidence"] = (ev > 0).sum(axis=1)
    out = out.join(mine[["sources", "category", "ot_ra_score", "RA_drugs", "late_stage_inhibitors"]])
    out["novelty"] = np.where(out["ot_ra_score"] < 0.1, "new to RA (OT < 0.1)",
                     np.where(out["ot_ra_score"] < 0.3, "weakly known", "established"))
    out["protein_p_best"] = pd.Series(prot_p)
    out = out.sort_values(["evidence_score", "n_layers_with_evidence"], ascending=False)
    out.index.name = "gene"
    out.round(3).to_csv(TAB / "gene_prioritisation.tsv", sep="\t")
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 50)
    print(out.head(25)[LAYERS + ["evidence_score", "category", "novelty"]].round(2).to_string())

    top = out.head(30)
    fig, ax = plt.subplots(figsize=(11, 10))
    im = ax.imshow(top[LAYERS].values, aspect="auto", cmap="Reds", vmin=0, vmax=1)
    ax.set_xticks(range(len(LAYERS)), [l.split("_", 1)[1].replace("_", " ") for l in LAYERS], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(top)), [f"{g}  ({s:.1f}; {c.split(':')[0][:28]}; {n.split(' (')[0]})"
                                    for g, s, c, n in zip(top.index, top.evidence_score, top.category, top.novelty)], fontsize=7)
    ax.set_title("Top 30 genes: evidence by layer (score; druggability; novelty in RA)", fontsize=10)
    plt.colorbar(im, ax=ax, shrink=0.5)
    fig.tight_layout(); fig.savefig(FIG / "gene_prioritisation.png", dpi=150)


if __name__ == "__main__":
    main()
