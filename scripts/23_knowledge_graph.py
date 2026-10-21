"""Knowledge graph of the workflow's findings (nodes and typed edges with provenance), plus an interactive viewer.

Node types: Gene, CellState, Variant, TF, Pathway, Score, Drug, Disease. Every edge carries the step that produced it, the
statistic and its value. Only results that passed each step's own threshold are included.

Output: results/knowledge_graph/nodes.tsv, edges.tsv, ra_multiomics_kg.graphml, kg_viewer.html
"""
import json
import shutil

import networkx as nx
import numpy as np
import pandas as pd

from config import ROOT, TAB

OUT = ROOT / "results" / "knowledge_graph"
OUT.mkdir(parents=True, exist_ok=True)
M6 = "M-6: STAT1+ CXCL10+"

nodes, edges = {}, []


def node(nid, ntype, label=None, **attrs):
    if nid not in nodes:
        nodes[nid] = dict(id=nid, type=ntype, label=label or nid, **attrs)
    return nid


def edge(src, dst, rel, step, statistic="", value=None, **attrs):
    edges.append(dict(source=src, target=dst, relation=rel, step=step, statistic=statistic,
                      value=(round(float(value), 4) if value is not None and pd.notna(value) else None), **attrs))


def main():
    ra = node("RA", "Disease", "Rheumatoid arthritis")
    lc = node("NSCLC", "Disease", "Non-small cell lung cancer")
    s2 = node("score:S2_IFN10", "Score", "S2 predefined 10-gene IFN score")
    for gene in ["CXCL9","CXCL10","CXCL11","STAT1","GBP1","GBP5","IDO1","IRF1","CD274","APOL6"]:
        edge(node(f"gene:{gene}", "Gene", gene), s2, "member_of_predefined_score", "24/30")

    pr = pd.read_csv(TAB / "gene_prioritisation.tsv", sep="\t", index_col=0)
    for g, r in pr.iterrows():
        node(f"gene:{g}", "Gene", g, evidence_score=float(r.evidence_score), druggability=r.category, novelty=r.novelty)

    de = pd.read_csv(TAB / "emtab8322_DE_active_vs_remission.tsv", sep="\t", index_col=0)
    for g, r in de[(de.padj < 0.05) & (de.log2FoldChange > 0.5) & ~de.index.str.match(r"^(ENSG|MT-$)")].iterrows():
        edge(node(f"gene:{g}", "Gene", g), ra, "up_in_active_RA_macrophages", "04", "log2FC", r.log2FoldChange, padj=float(r.padj))

    node(M6, "CellState", "M-6 STAT1+ CXCL10+ macrophage")
    for g, r in pr[pr.L2_tissue_M6 > 0].iterrows():
        edge(f"gene:{g}", M6, "marks_cell_state", "22", "M-6 enrichment", 1)

    gw = pd.read_csv(TAB / "gwas_enrichment.tsv", sep="\t")
    gw = gw[(gw.trait == "RA") & (gw.window_kb == 25) & (gw.fdr < 0.1)]
    for _, r in gw.iterrows():
        edge(node(r.state, "CellState"), ra, "accessible_chromatin_enriched_for_RA_risk_loci", "10", "enrichment",
             r.enrichment, fdr=float(r.fdr))

    cv = pd.read_csv(TAB / "gwas_candidate_variants.tsv", sep="\t")
    for _, r in cv.drop_duplicates(["rs_id", "state"]).iterrows():
        v = node(f"snp:{r.rs_id}", "Variant", r.rs_id, location=r.location)
        edge(v, node(r.state, "CellState"), "in_state_specific_peak", "10", "GWAS p", r.p_value)
    for _, r in cv.drop_duplicates("rs_id").iterrows():
        edge(f"snp:{r.rs_id}", node(f"gene:{r.nearest_protein_coding_gene}", "Gene", r.nearest_protein_coding_gene),
             "nearest_protein_coding_gene", "10")
        edge(f"snp:{r.rs_id}", ra, "associated_with", "09", "GWAS p", r.p_value)

    vm = pd.read_csv(TAB / "variant_motif_effects.tsv", sep="\t")
    vm = vm[vm.tf_expressed_in_state & ~vm.tf.str.upper().str.match(r"^(ZNF|ZSCAN|ZKSCAN|ZBTB|ZBED|ZFP)")]
    for _, r in vm.drop_duplicates(["rs_id", "tf"]).iterrows():
        edge(f"snp:{r.rs_id}", node(f"tf:{r.tf}", "TF", r.tf), "risk_allele_" + str(r.risk_allele_effect).replace(" ", "_"),
             "16", "delta log-odds (alt - ref)", r.delta_alt_minus_ref)

    cb = pd.read_csv(TAB / "chromvar_M6_vs_macrophages.tsv", sep="\t", index_col=0).head(10)
    for k, r in cb.iterrows():
        edge(node(f"tf:{r.motif}", "TF", r.motif), M6, "motif_more_accessible_in", "14", "mean delta z", r.mean_diff_z,
             n_donors=int(r.n_donors))

    li = pd.read_csv(TAB / "liana_to_M6_specific.tsv", sep="\t").head(20)
    ifng = pd.read_csv(TAB / "liana_results.tsv.gz", sep="\t")
    ifng = ifng[(ifng.ligand_complex == "IFNG") & (ifng.target == M6)].sort_values("specificity_rank").head(3)
    cx = pd.read_csv(TAB / "liana_results.tsv.gz", sep="\t")
    cx = cx[(cx.source == M6) & cx.ligand_complex.isin(["CXCL9", "CXCL10", "CXCL11"]) & (cx.receptor_complex == "CXCR3")] \
        .sort_values("magnitude_rank").head(3)
    for _, r in pd.concat([li, ifng, cx]).iterrows():
        edge(node(r.source, "CellState"), node(r.target, "CellState"), "predicted_ligand_receptor_interaction", "13", "magnitude rank",
             r.magnitude_rank, ligand=r.ligand_complex, receptor=r.receptor_complex)

    gs = pd.read_csv(TAB / "gsea_results.tsv", sep="\t")
    h = gs[(gs.library == "hallmark") & (gs.fdr < 0.05)]
    for _, r in h.iterrows():
        p = node(f"pw:{r.pathway}", "Pathway", r.pathway.replace("HALLMARK_", ""))
        if r.contrast.startswith("A") and r.nes > 0:
            edge(p, ra, "up_in_active_RA_macrophages", "12", "NES", r.nes, fdr=float(r.fdr))
        if r.contrast.startswith("B") and r.nes < 0:
            edge(p, node("anti-TNF good responders", "Disease", "Anti-TNF cohort: good responders"), "lower_post_than_pre_in_good_responders", "12", "NES", r.nes, fdr=float(r.fdr))
    le = [set(str(x).split(";")) for x in gs[(gs.pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE") & (gs.library == "hallmark")]["Lead_genes"]]
    for g in sorted(set.intersection(*le)) if len(le) == 2 else []:
        if f"gene:{g}" in nodes:
            edge(f"gene:{g}", "pw:HALLMARK_INTERFERON_GAMMA_RESPONSE", "leading_edge_member", "12")

    mine = pd.read_csv(TAB / "target_drug_mining.tsv", sep="\t")
    for _, r in mine.iterrows():
        for field, rel in (("late_stage_inhibitors", "late_stage_inhibitor_of"), ("RA_drugs", "tested_or_approved_in_RA_targeting")):
            if isinstance(r[field], str) and r[field]:
                for d in r[field].split("; ")[:4]:
                    name = d.split(" (")[0]
                    edge(node(f"drug:{name}", "Drug", name), f"gene:{r.gene}", rel, "15", "", None,
                         stage=d.split(" (")[1].rstrip(")") if "(" in d else "")

    cox = pd.read_csv(TAB / "nsclc_cox.tsv", sep="\t")
    for _, r in cox[cox.score == "polarity_CXCL9_minus_SPP1"].iterrows():
        for g in ("CXCL9", "SPP1"):
            edge(node(f"gene:{g}", "Gene", g), lc, f"macrophage_polarity_axis_{r.cancer}", "19", "HR per SD", r.HR_per_SD,
                 p=float(r.p))

    cpath = TAB / "eqtl_coloc.tsv"
    if cpath.exists():
        co = pd.read_csv(cpath, sep="\t").dropna(subset=["PP_H4"])
        for _, r in co[co.PP_H4 > 0.5].iterrows():
            edge(node(f"gene:{r.gene}", "Gene", r.gene), ra, "eQTL_colocalises_with_RA_GWAS", "25", "PP.H4", r.PP_H4,
                 eqtl_dataset=r.label, variant=r.top_shared_variant)
    vpath = TAB / "visium_sections.tsv"
    if vpath.exists():
        v = pd.read_csv(vpath, sep="\t")
        edge(s2, ra,
             "S2_IFN_high_spots_associated_with_RA_tissue", "24", "median % IFN-high spots RA",
             v.loc[v.condition == "RA", "pct_spots_ifn_high"].median())
    fpath = TAB / "finemap_enrichment.tsv"
    if fpath.exists():
        fe = pd.read_csv(fpath, sep="\t")
        for _, r in fe[(fe.study != "GCST90018961") & (fe.fdr < 0.1)].iterrows():
            edge(node(r.state, "CellState"), ra, "chromatin_enriched_for_fine_mapped_RA_variants", "29", "PIP enrichment",
                 r.enrichment, fdr=float(r.fdr), gwas=r.label)
        fc = pd.read_csv(TAB / "finemap_candidates.tsv", sep="\t").sort_values("pip", ascending=False)
        for _, r in fc.drop_duplicates(["variant", "state"]).iterrows():
            rs = str(r.rsid).split(";")[-1] if isinstance(r.rsid, str) and r.rsid else r.variant
            v = node(f"snp:{rs}", "Variant", rs, location=r.variant)
            edge(v, node(r.state, "CellState"), "fine_mapped_variant_in_state_specific_peak", "29", "PIP", r.pip, gwas=r.study)
            edge(v, node(f"gene:{r.nearest_protein_coding_gene}", "Gene", r.nearest_protein_coding_gene),
                 "nearest_protein_coding_gene", "29")
    spath = TAB / "gse89408_tests.tsv"
    if spath.exists():
        st = pd.read_csv(spath, sep="\t").set_index(["score", "test"])
        r = st.loc[("S2", "H2 RA vs OA")]
        if r.p_vs_random < 0.05:
            edge(s2, ra, "S2_IFN_score_associated_with_RA_vs_OA_synovium",
                 "30", "t (adjusted)", r.t, p_vs_random=float(r.p_vs_random))
    tpath = TAB / "myeloid_paga.tsv"
    if tpath.exists():
        names = {s.split(":")[0]: s for s in pd.read_csv(TAB / "signature_by_myeloid_subset.tsv", sep="\t").subset}
        pg = pd.read_csv(tpath, sep="\t", index_col=0).drop(index="M-6")
        for col, rel in (("rna_paga_to_M6", "nearest_state_on_RNA_graph"), ("atac_paga_to_M6", "nearest_state_on_chromatin_graph")):
            k = pg[col].idxmax()
            edge(M6, node(names.get(k, k), "CellState"), rel, "31", "PAGA connectivity", pg.loc[k, col])
    dpath = TAB / "coloc_direction.tsv"
    if dpath.exists():
        dr = pd.read_csv(dpath, sep="\t")
        for _, r in dr[dr.risk_allele_effect.isin(["raises expression", "lowers expression"]) & (dr.immune_datasets >= 2)].iterrows():
            edge(node(f"gene:{r.gene}", "Gene", r.gene), ra, "risk_allele_" + r.risk_allele_effect.replace(" ", "_") + "_in_immune_cells",
                 "32", "share of immune QTL datasets where risk raises level", r.frac_risk_raises_expression,
                 datasets=int(r.immune_datasets))
    spath2 = TAB / "sldsc_results.tsv"
    if spath2.exists():
        sl = pd.read_csv(spath2, sep="\t")
        safe = lambda x: "".join(ch if ch.isalnum() else "_" for ch in x).strip("_")      # as in step 33
        by_safe = {safe(n): n for n, v in nodes.items() if v["type"] == "CellState"}
        for _, r in sl[sl.gwas.str.startswith("RA") & (sl.fdr_tau < 0.1)].iterrows():
            target = by_safe.get(r.annotation[3:], r.annotation[3:])
            edge(node(target, "CellState"), ra, "heritability_enriched_S_LDSC", "33", "tau z", r.tau_z, fdr=float(r.fdr_tau), gwas=r.gwas)
    tpath2 = TAB / "tf_M6_drivers.tsv"
    if tpath2.exists():
        tf = pd.read_csv(tpath2, sep="\t", index_col=0)
        for t, r in tf[tf.both_layers & ~tf.complex].iterrows():
            edge(node(f"tf:{t}", "TF", t), M6, "regulon_and_motif_more_active_in", "34", "RNA activity difference",
                 r.rna_activity_diff, motif_diff_z=float(r.motif_diff_z))
    mpath = TAB / "mr_primary.tsv"
    if mpath.exists():
        mr = pd.read_csv(mpath, sep="\t")
        mr = mr[(mr.n_datasets >= 2) & (mr.share_same_sign >= 0.75) & ((mr.OR_lo > 1) | (mr.OR_hi < 1))]
        for _, r in mr.iterrows():
            rel = "higher_level_associated_with_higher_RA_risk_Wald" if r.OR_per_SD > 1 else "higher_level_associated_with_lower_RA_risk_Wald"
            edge(node(f"gene:{r.gene}", "Gene", r.gene), ra, rel, "37", "OR per SD", r.OR_per_SD, gwas=r.gwas, dataset=r.dataset)
    dpath2 = TAB / "scdrs_groups.tsv"
    if dpath2.exists():
        sd = pd.read_csv(dpath2, sep="\t")
        sd = sd[sd.trait.str.startswith("RA") & (sd.level == "subtype")]
        both = sd.groupby("group").assoc_mcp.max()
        for g, pmax in both[both <= 0.001].items():
            z = sd[sd.group == g].assoc_mcz.min()
            edge(node(g, "CellState"), ra, "expresses_RA_risk_genes_scDRS", "36", "min MC z over RA GWAS", z)
    spath3 = TAB / "tph_spatial.tsv"
    if spath3.exists():
        tp = pd.read_csv(spath3, sep="\t")
        if len(tp) and (tp.p_B < 0.05).all():
            edge(node("CXCL13-positive T cells (Xenium)", "CellState"), node("B cells (Xenium)", "CellState"), "spatial_neighbour_association",
                 "35", "fraction sections empirical P <0.05", float((tp.p_B<.05).mean()), sections=len(tp), evidence_type="global label shuffle; not causal")
    try:
        xt = pd.read_csv(TAB / "xenium_tests.tsv", sep="\t")
        h1 = xt.iloc[0]
        edge(node("IFNG-positive T/NK cells (Xenium)", "CellState"), node("S2 IFN-scored macrophages (Xenium)", "CellState"), "spatial_neighbour_association",
             "20", "fraction sections empirical P <0.05", float(h1.n_perm_p_lt_0_05/h1.n_sections), sections=int(h1.n_sections), evidence_type="global label shuffle; not causal")
    except FileNotFoundError:
        pass

    N = pd.DataFrame(nodes.values())
    E = pd.DataFrame(edges)
    E = E[E.source.isin(N.id) & E.target.isin(N.id)].drop_duplicates(["source", "target", "relation"])
    N.to_csv(OUT / "nodes.tsv", sep="\t", index=False)
    E.to_csv(OUT / "edges.tsv", sep="\t", index=False)
    G = nx.MultiDiGraph()
    for _, r in N.iterrows():
        G.add_node(r.id, **{k: ("" if pd.isna(v) else v) for k, v in r.items() if k != "id"})
    for _, r in E.iterrows():
        G.add_edge(r.source, r.target, **{k: ("" if pd.isna(v) else v) for k, v in r.items() if k not in ("source", "target")})
    nx.write_graphml(G, OUT / "ra_multiomics_kg.graphml")
    print(f"Knowledge graph: {len(N)} nodes ({N.type.value_counts().to_dict()}), {len(E)} edges "
          f"({E.relation.value_counts().head(8).to_dict()})")

    colours = {"Gene": "#d62728", "CellState": "#1f77b4", "Variant": "#9467bd", "TF": "#ff7f0e", "Pathway": "#2ca02c",
               "Score": "#CC79A7", "Drug": "#8c564b", "Disease": "#111111"}
    connected = set(E.source) | set(E.target)                    # the viewer hides isolated nodes; TSV/GraphML keep them
    vis_nodes = [dict(id=r.id, label=r.label, group=r.type, color=colours.get(r.type, "#999"),
                      title=f"{r.type}: {r.label}" + (f"<br>evidence score {r.get('evidence_score')}" if r.type == "Gene" and pd.notna(r.get("evidence_score")) else ""),
                      value=float(r.get("evidence_score", 1) if pd.notna(r.get("evidence_score", None)) else 1))
                 for _, r in N.iterrows() if r.id in connected]
    vis_edges = [dict(**{"from": r.source, "to": r.target}, label="", arrows="to",
                      title=f"{r.relation} (step {r.step}){'; ' + r.statistic + ' = ' + str(r.value) if r.statistic else ''}")
                 for _, r in E.iterrows()]
    html = (VIEWER.replace("__NODES__", json.dumps(vis_nodes)).replace("__EDGES__", json.dumps(vis_edges))
            .replace("__LEGEND__", " ".join(f'<span style="color:{c}">&#9679; {t}</span>' for t, c in colours.items())))
    vendor=ROOT/'references/vendor'
    for asset in ['vis-network-9.1.9.min.js','LICENSE-MIT.txt']:
        shutil.copy2(vendor/asset,OUT/asset)
    (OUT / "kg_viewer.html").write_text(html)
    print(f"Viewer: {OUT / 'kg_viewer.html'}")


VIEWER = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>RA Multi-omics Graph</title>
<script src="vis-network-9.1.9.min.js"></script>
<style>
:root { --bg:#ffffff; --fg:#1a1a1a; --muted:#666; --border:#ddd; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { --bg:#121212; --fg:#eee; --muted:#aaa; --border:#333; } }
:root[data-theme="dark"] { --bg:#121212; --fg:#eee; --muted:#aaa; --border:#333; }
body { margin:0; font-family: system-ui, sans-serif; background:var(--bg); color:var(--fg); }
header { padding:12px 16px; border-bottom:1px solid var(--border); }
header h1 { font-size:18px; margin:0 0 4px; } header p { margin:0; color:var(--muted); font-size:13px; }
#controls { padding:8px 16px; display:flex; gap:8px; flex-wrap:wrap; align-items:center; font-size:13px; }
#controls input { padding:4px 8px; border:1px solid var(--border); border-radius:6px; background:var(--bg); color:var(--fg); }
#net { height: calc(100vh - 150px); border-top:1px solid var(--border); }
#info { padding:6px 16px; font-size:12px; color:var(--muted); }
</style></head><body>
<header><h1>RA synovium multi-omics knowledge graph</h1>
<p>Findings of the GeNeDis 2026 RA workflow. Hover a node or edge for its evidence; every edge names the pipeline step that produced it.</p></header>
<div id="controls">__LEGEND__ &nbsp; <input id="q" placeholder="Find a node (e.g. SLAMF7)"></div>
<div id="info">Node size = integrated evidence score (genes). Drag to explore; scroll to zoom.</div>
<div id="net"></div>
<script>
const nodes = new vis.DataSet(__NODES__);
const edges = new vis.DataSet(__EDGES__);
const net = new vis.Network(document.getElementById('net'), {nodes, edges},
  {nodes:{shape:'dot', scaling:{min:6,max:30}, font:{size:12, color: getComputedStyle(document.body).color}},
   edges:{color:{opacity:0.4}, smooth:false, width:0.8},
   layout:{improvedLayout:false},
   physics:{stabilization:{iterations:400}, barnesHut:{gravitationalConstant:-6000, springLength:120}}});
net.once('stabilizationIterationsDone', () => net.setOptions({physics:false}));
function find(t) {
  t = t.trim().toLowerCase(); if (!t) return;
  const hit = nodes.get().find(n => n.label.toLowerCase() === t) || nodes.get().find(n => n.label.toLowerCase().includes(t));
  if (!hit) { document.getElementById('info').textContent = 'No node matches "' + t + '".'; return; }
  const nb = net.getConnectedNodes(hit.id);
  net.selectNodes([hit.id, ...nb]);
  net.focus(hit.id, {scale:1.6, animation:{duration:600}});
  document.getElementById('info').textContent = hit.label + ': ' + nb.length + ' connected nodes: ' +
    nb.map(id => nodes.get(id).label).slice(0, 25).join(', ');
}
const q = document.getElementById('q');
q.addEventListener('keydown', e => { if (e.key === 'Enter') find(q.value); });
q.addEventListener('change', () => find(q.value));
</script></body></html>
"""

if __name__ == "__main__":
    main()
