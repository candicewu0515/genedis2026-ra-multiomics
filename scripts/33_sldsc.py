"""Stratified LD score regression (S-LDSC): do cell-state-specific accessible regions carry more RA heritability
than expected, given the baselineLD model and accessible chromatin in general?

Annotations (GRCh38, peak +/- 500 bp):
  ALLPEAKS        all 131,644 non-MHC synovial peaks (background for every model)
  broad types     top 10 % specific peaks per author cell type (13,164 peaks; same specificity score as step 10)
  subtypes        the top-4,000 specific peaks per subtype used in steps 10 and 29
Reference: 1000 Genomes EUR GRCh38 plink files, baselineLD v2.2 and HapMap3 weights (Broad, Zenodo 10515792);
LD scores with a 1 cM window, printed for HapMap3 SNPs without MHC.
GWAS: Ishigaki 2022 RA EUR (GCST90132223, effective N 68,838); FinnGen R12 RA (M13_RHEUMA, effective N 62,044);
control trait FinnGen R12 type 2 diabetes (T2D, effective N 275,022). Munged to HapMap3 SNPs, strand-ambiguous
SNPs removed, MAF >= 0.01 where frequencies are given.
Model per state: baselineLD v2.2 + ALLPEAKS + state (ldsc --h2 --overlap-annot). Test: the state's coefficient (tau),
one-sided p from its z score; BH within GWAS x level. Enrichment (share of h2 / share of SNPs) reported alongside.
Hypotheses fixed before analysis:
  S1 T-cell-specific peaks carry RA heritability beyond baselineLD + all peaks (Ishigaki), replicated in FinnGen.
  S2 Within lineages (exploratory): which subtypes; are STAT1+ CXCL10+ macrophage peaks (M-6) enriched? (M-6 has
     65 nuclei, below the 100-nucleus minimum for a peak set, so it cannot be tested.)
  S3 control: T2D shows no immune-state enrichment.
LDSC: Python 3 port in polyfun (commit 8b17b4e) run in the separate environment of environment_ldsc.yml
(LDSC_PYTHON, default the genedis_ldsc env).

Output: data/interim/ldsc/ (annotations, LD scores, munged sumstats), results/tables/sldsc_results.tsv,
        sldsc_allpeaks_summary.tsv,
        results/figures/sldsc.png
"""
import gzip
import importlib
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

from config import ROOT, RAW, INTERIM, FIG, TAB

step10 = importlib.import_module("10_gwas_peak_enrichment")
# the LDSC environment is expected next to the active one (…/envs/genedis_ldsc); override with LDSC_PYTHON
LDSC_PY = os.environ.get("LDSC_PYTHON", str(Path(sys.executable).resolve().parents[2] / "genedis_ldsc" / "bin" / "python"))
LDSC = ROOT / "tools" / "polyfun" / "ldsc.py"
REF = RAW / "ldsc" / "GRCh38"
HM3 = RAW / "ldsc" / "hm3_no_MHC.list.txt"
W = INTERIM / "ldsc"
EXT = 500
TOP_FRAC = 0.10
N_JOBS = 4                 # ldsc itself is multithreaded
MHC = ("chr6", 28_510_120, 33_480_577)
GWAS = {"RA_Ishigaki_EUR": dict(file=RAW / "gwas_sumstats" / "GCST90132223.h.tsv.gz", fmt="gwascat", ncas=22350, ncon=74823),
        "RA_FinnGen_R12": dict(file=RAW / "gwas_sumstats" / "finngen_R12_M13_RHEUMA.gz", fmt="finngen", ncas=16314, ncon=315115),
        "T2D_FinnGen_R12_control": dict(file=RAW / "gwas_sumstats" / "finngen_R12_T2D.gz", fmt="finngen", ncas=82878, ncon=403489)}


def safe(name):
    return "".join(ch if ch.isalnum() else "_" for ch in name).strip("_")


def peak_sets():
    """ALLPEAKS, top-10 % broad types, and the step-10 subtype sets, as (chrom, start, end) arrays."""
    atac = sc.read_h5ad(INTERIM / "gse243917_atac.h5ad")
    peaks = atac.var[["chrom", "start", "end"]].copy()
    ok = ~((peaks.chrom == MHC[0]) & (peaks.end >= MHC[1]) & (peaks.start <= MHC[2])).values
    sets = step10.specific_sets(atac, ok)
    pk = peaks[ok].reset_index(drop=True)
    out = {"ALLPEAKS": pk}
    # broad types: same specificity score as step 10, but top 10 % instead of 4,000
    C = atac.layers["counts"].tocsr()[:, ok]
    ct = atac.obs["author_ct"].astype(str).values
    cts = [c for c in pd.unique(ct) if (ct == c).sum() >= step10.MIN_STATE]
    cpm = step10.pseudobulk_cpm(C, ct, cts)
    n_top = int(TOP_FRAC * C.shape[1])
    for i, c in enumerate(cts):
        others = np.delete(cpm, i, axis=0).mean(axis=0)
        score = np.log2((cpm[i] + 1) / (others + 1))
        score[cpm[i] < np.quantile(cpm[i], 0.5)] = -np.inf
        out[f"CT_{safe(c)}"] = pk.iloc[np.argsort(score)[::-1][:n_top]]
    for (level, state), idx in sets.items():
        if level == "subtype":
            out[f"ST_{safe(state)}"] = pk.iloc[idx]
    return out


def write_annotations(sets):
    (W / "annot").mkdir(parents=True, exist_ok=True)
    names = list(sets)
    for chrom in range(1, 23):
        f = W / "annot" / f"states.{chrom}.annot.gz"
        bim = pd.read_csv(REF / "plink_files" / f"1000G.EUR.hg38.{chrom}.bim", sep="\t", header=None,
                          names=["chr", "snp", "cm", "bp", "a1", "a2"])
        bp = bim.bp.values
        cols = {}
        for n in names:
            s = sets[n][sets[n].chrom == f"chr{chrom}"]
            starts, ends = np.sort(s.start.values - EXT), np.sort(s.end.values + EXT)
            # a SNP is inside some interval if #starts <= bp exceeds #ends < bp (intervals may overlap)
            cols[n] = ((np.searchsorted(starts, bp, side="right") - np.searchsorted(ends, bp, side="left")) > 0).astype(np.int8)
        pd.DataFrame(cols).to_csv(f, sep="\t", index=False, compression="gzip")
    return names


def run(cmd, log):
    env = dict(os.environ, OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="4", MKL_NUM_THREADS="4")   # 4 jobs x 4 threads
    with open(log, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=LDSC.parent, env=env)
    return r.returncode


def ld_scores():
    (W / "l2").mkdir(parents=True, exist_ok=True)
    snps = W / "hm3_snps.txt"
    if not snps.exists():
        pd.read_csv(HM3, header=None)[0].to_csv(snps, index=False, header=False)
    jobs = []
    for chrom in range(1, 23):
        out = W / "l2" / f"states.{chrom}"
        if (W / "l2" / f"states.{chrom}.l2.ldscore.gz").exists():
            continue
        jobs.append(([LDSC_PY, str(LDSC), "--l2", "--bfile", str(REF / "plink_files" / f"1000G.EUR.hg38.{chrom}"),
                      "--ld-wind-cm", "1", "--annot", str(W / "annot" / f"states.{chrom}.annot.gz"), "--thin-annot",
                      "--out", str(out), "--print-snps", str(snps)], str(out) + ".runlog"))
    with ThreadPoolExecutor(N_JOBS) as ex:
        codes = list(ex.map(lambda j: run(*j), jobs))
    if any(codes):
        raise RuntimeError(f"ldsc --l2 failed for {sum(c != 0 for c in codes)} chromosomes; see {W / 'l2'}/*.runlog")


def split_ld(names):
    """One LD-score file set per annotation so each model can include ALLPEAKS + one state."""
    (W / "split").mkdir(parents=True, exist_ok=True)
    for chrom in range(1, 23):
        l2 = pd.read_csv(W / "l2" / f"states.{chrom}.l2.ldscore.gz", sep="\t")
        M = np.loadtxt(W / "l2" / f"states.{chrom}.l2.M", ndmin=1)
        M550 = np.loadtxt(W / "l2" / f"states.{chrom}.l2.M_5_50", ndmin=1)
        ann = pd.read_csv(W / "annot" / f"states.{chrom}.annot.gz", sep="\t")      # thin: rows follow the .bim / .frq
        for k, n in enumerate(names):
            l2[["CHR", "SNP", "BP", f"{n}L2"]].to_csv(W / "split" / f"{n}.{chrom}.l2.ldscore.gz", sep="\t", index=False,
                                                      compression="gzip")
            (W / "split" / f"{n}.{chrom}.l2.M").write_text(f"{M[k]}\n")
            (W / "split" / f"{n}.{chrom}.l2.M_5_50").write_text(f"{M550[k]}\n")
            ann[[n]].to_csv(W / "split" / f"{n}.{chrom}.annot.gz", sep="\t", index=False, compression="gzip")


def munge(name, spec):
    out = W / "sumstats" / f"{name}.sumstats.gz"
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    hm3 = set(pd.read_csv(HM3, header=None)[0])
    neff = 4 / (1 / spec["ncas"] + 1 / spec["ncon"])
    parts = []
    if spec["fmt"] == "gwascat":
        it = pd.read_csv(spec["file"], sep="\t", chunksize=1_000_000, dtype=str)
        for d in it:
            rs = d["rsid"] if "rsid" in d else d["variant_id"]
            d = d.assign(SNP=rs, CHR=d.chromosome, BP=d.base_pair_location, A1=d.effect_allele, A2=d.other_allele,
                         B=pd.to_numeric(d.beta, errors="coerce"),
                         SE=pd.to_numeric(d.standard_error, errors="coerce"),
                         F=pd.to_numeric(d.get("effect_allele_frequency"), errors="coerce"))
            parts.append(d[["SNP", "CHR", "BP", "A1", "A2", "B", "SE", "F"]][d.SNP.isin(hm3)])
    else:
        it = pd.read_csv(spec["file"], sep="\t", chunksize=1_000_000, dtype=str)
        for d in it:
            d = d.assign(SNP=d.rsids.str.split(",").str[0], CHR=d["#chrom"], BP=d.pos, A1=d.alt, A2=d.ref,
                         B=pd.to_numeric(d.beta, errors="coerce"),
                         SE=pd.to_numeric(d.sebeta, errors="coerce"), F=pd.to_numeric(d.af_alt, errors="coerce"))
            parts.append(d[["SNP", "CHR", "BP", "A1", "A2", "B", "SE", "F"]][d.SNP.isin(hm3)])
    s = pd.concat(parts).dropna(subset=["B", "SE"])
    s = s[(s.SE > 0) & s.A1.str.len().eq(1) & s.A2.str.len().eq(1)]
    amb = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
    s = s[[(a, b) not in amb for a, b in zip(s.A1.str.upper(), s.A2.str.upper())]]
    s = s[s.F.isna() | s.F.between(0.01, 0.99)]
    s = s.drop_duplicates("SNP", keep=False)
    s = pd.DataFrame({"SNP": s.SNP, "CHR": s.CHR, "BP": s.BP, "A1": s.A1.str.upper(), "A2": s.A2.str.upper(),
                      "Z": s.B / s.SE, "N": round(neff)})          # polyfun's LDSC also requires CHR and BP
    s.to_csv(out, sep="\t", index=False, compression="gzip")
    print(f"  munged {name}: {len(s):,} HapMap3 SNPs, N_eff {neff:,.0f}, mean chi2 {np.mean(s.Z ** 2):.3f}")
    return out


def regress(gw, ss, n):
    out = W / "h2" / f"{gw}__{n}"
    if Path(str(out) + ".results").exists():
        return 0
    (W / "h2").mkdir(parents=True, exist_ok=True)
    ref = ",".join([str(REF / "baselineLD_v2.2" / "baselineLD."), str(W / "split" / "ALLPEAKS."), str(W / "split" / f"{n}.")])
    cmd = [LDSC_PY, str(LDSC), "--h2", str(ss), "--ref-ld-chr", ref,
           "--w-ld-chr", str(REF / "weights" / "weights.hm3_noMHC."), "--overlap-annot",
           "--frqfile-chr", str(REF / "plink_files" / "1000G.EUR.hg38."), "--print-coefficients", "--out", str(out)]
    return run(cmd, str(out) + ".runlog")


def validate_complete_results(res, tests):
    """Keep the fixed 3-GWAS × 65-annotation family before BH adjustment."""
    expected={(gw,n) for gw,_,n in tests}
    if len(tests)!=195 or len(expected)!=195:
        raise RuntimeError('Expected 195 unique S-LDSC models for the fixed atlas')
    if len(res)!=195 or res.duplicated(['gwas','annotation']).any() or set(zip(res.gwas,res.annotation))!=expected:
        raise RuntimeError('Incomplete or duplicate S-LDSC results; no partial-family FDR is allowed')
    numeric=['prop_snps','prop_h2','enrichment','tau','tau_se','tau_z','p_tau']
    if not np.isfinite(res[numeric].to_numpy(dtype=float)).all():
        raise RuntimeError('Non-finite S-LDSC statistics; no partial-family FDR is allowed')
    if not res.tau_se.gt(0).all() or not res.p_tau.between(0,1).all() or not res.prop_snps.gt(0).all():
        raise RuntimeError('Invalid S-LDSC standard errors, P values or annotation sizes')


def main():
    W.mkdir(parents=True, exist_ok=True)
    if not (W / "annot" / "states.22.annot.gz").exists():
        sets = peak_sets()
        names = write_annotations(sets)
        pd.DataFrame({"annotation": names, "n_peaks": [len(sets[n]) for n in names]}).to_csv(W / "annotations.tsv", sep="\t", index=False)
    ann = pd.read_csv(W / "annotations.tsv", sep="\t")
    names = ann.annotation.tolist()
    print(f"Annotations: {len(names)} ({sum(n.startswith('CT_') for n in names)} broad types, "
          f"{sum(n.startswith('ST_') for n in names)} subtypes, ALLPEAKS)")
    # Value-verified native storage cache; statistical inputs are unchanged.
    subprocess.run([LDSC_PY, str(Path(__file__).with_name("prepare_ldsc_parquet.py"))], check=True)
    ld_scores()
    if not (W / "split" / f"{names[-1]}.22.annot.gz").exists():
        split_ld(names)
    sums = {gw: munge(gw, spec) for gw, spec in GWAS.items()}
    tests = [(gw, ss, n) for gw, ss in sums.items() for n in names if n != "ALLPEAKS"]
    with ThreadPoolExecutor(N_JOBS) as ex:
        codes = list(ex.map(lambda t: regress(*t), tests))
    if any(codes):
        raise RuntimeError(f"{sum(c != 0 for c in codes)} S-LDSC regressions failed; see {W / 'h2'}/*.runlog")

    rows, allp, h2 = [], [], []
    for gw, ss, n in tests:
        f = W / "h2" / f"{gw}__{n}.results"
        if not f.exists():
            raise RuntimeError("Missing S-LDSC result: "+str(f))
        r = pd.read_csv(f, sep="\t").set_index("Category")
        a = r.loc[[c for c in r.index if c.startswith("ALLPEAKS")][0]]
        allp.append(dict(gwas=gw, enrichment=a["Enrichment"], enrichment_p=a["Enrichment_p"], tau_z=a["Coefficient_z-score"]))
        lg = (W / "h2" / f"{gw}__{n}.log").read_text()
        h2.append(dict(gwas=gw, h2=lg.split("Total Observed scale h2: ")[1].split()[0] if "Total Observed" in lg else "",
                       intercept=lg.split("Intercept: ")[1].split()[0] if "Intercept: " in lg else ""))
        key = [c for c in r.index if c.startswith(n)][0]
        x = r.loc[key]
        rows.append(dict(gwas=gw, level="cell type" if n.startswith("CT_") else "subtype", annotation=n,
                         prop_snps=x["Prop._SNPs"], prop_h2=x["Prop._h2"], enrichment=x["Enrichment"],
                         enrichment_p=x["Enrichment_p"], tau=x["Coefficient"], tau_se=x["Coefficient_std_error"],
                         tau_z=x["Coefficient_z-score"], p_tau=norm.sf(x["Coefficient_z-score"])))
    res = pd.DataFrame(rows)
    validate_complete_results(res, tests)
    res["fdr_tau"] = np.nan
    for _, g in res.groupby(["gwas", "level"]):
        res.loc[g.index, "fdr_tau"] = multipletests(g.p_tau, method="fdr_bh")[1]
    res = res.sort_values(["gwas", "level", "p_tau"])
    res.to_csv(TAB / "sldsc_results.tsv", sep="\t", index=False, float_format="%.4g")
    pd.set_option("display.width", 220)
    for gw in GWAS:
        print(f"\n== {gw}")
        print(res[res.gwas == gw].groupby("level").head(6)[["level", "annotation", "prop_snps", "prop_h2", "enrichment",
                                                             "tau_z", "p_tau", "fdr_tau"]].round(4).to_string(index=False))
    ap = pd.DataFrame(allp).groupby("gwas").median()
    hh = pd.DataFrame(h2).drop_duplicates("gwas").set_index("gwas")
    summ = ap.join(hh)
    summ.to_csv(TAB / "sldsc_allpeaks_summary.tsv", sep="\t", float_format="%.4g")
    print("\nAll synovial peaks (median over models) and total h2:\n", summ.to_string())
    t = res[res.annotation == "CT_Tcell"].set_index("gwas")
    print("\nS1 T cells: " + "; ".join(f"{g} tau z {r.tau_z:.2f}, p {r.p_tau:.3g}, FDR {r.fdr_tau:.3g}, enrichment {r.enrichment:.1f}"
                                       for g, r in t.iterrows()))
    m6 = res[res.annotation.str.contains("M_6")]
    print("S2 M-6 peaks: " + ("; ".join(f"{r.gwas} tau z {r.tau_z:.2f}, p {r.p_tau:.3g}" for r in m6.itertuples()) or
                              f"not testable (M-6 has fewer than {step10.MIN_STATE} nuclei, so step 10 defines no peak set)"))
    ctl = res[res.gwas.str.startswith("T2D")]
    print(f"S3 T2D control: annotations with FDR < 0.05: {int((ctl.fdr_tau < 0.05).sum())} of {len(ctl)}")

    fig, axes = plt.subplots(1, 2, figsize=(16, 9), gridspec_kw={"width_ratios": [1, 2.4]})
    cols = {"RA_Ishigaki_EUR": "C3", "RA_FinnGen_R12": "C1", "T2D_FinnGen_R12_control": "C7"}
    for ax, level in zip(axes, ("cell type", "subtype")):
        p = res[res.level == level].pivot_table(index="annotation", columns="gwas", values="tau_z")
        f = res[res.level == level].pivot_table(index="annotation", columns="gwas", values="fdr_tau")
        p = p.sort_values("RA_Ishigaki_EUR")
        y = np.arange(len(p)); h = 0.8 / len(GWAS)
        for k, gw in enumerate(GWAS):
            off = h * (len(GWAS) / 2 - 0.5 - k)
            ax.barh(y + off, p[gw], height=h, color=cols[gw], label=gw)
            for i, a in enumerate(p.index):
                if f.loc[a, gw] < 0.05:
                    ax.text(max(p.loc[a, gw], 0) + 0.1, i + off, "*", va="center", fontsize=9)
        ax.axvline(0, color="k", lw=0.8); ax.axvline(norm.isf(0.05), color="k", lw=0.6, ls="--")
        ax.set_yticks(y, [a[3:] for a in p.index], fontsize=7 if level == "subtype" else 9)
        ax.set_xlabel("tau z (state beyond baselineLD + all peaks)")
        ax.set_title(f"S-LDSC, {level} (* FDR < 0.05; dashed: one-sided p 0.05)", fontsize=10)
        ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "sldsc.png", dpi=150)


if __name__ == "__main__":
    from pipeline_step_lock import step_lock
    with step_lock("33"):
        main()
