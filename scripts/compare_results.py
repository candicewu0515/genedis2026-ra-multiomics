"""Compare two runs of the pipeline table by table (reproducibility check).

Usage: python scripts/compare_results.py <baseline_tables_dir> [<current_tables_dir>] [--out results/reproducibility_check.tsv]
Each TSV present in either folder is classed as
  identical      same shape, same text, numbers equal to 1e-9 (relative)
  float noise    numbers differ by < 1e-3 relative (thread / BLAS order), text identical
  differs        larger numeric differences or changed text / shape (max differences reported)
  missing        present in one run only
Row order is ignored when the tables only differ by order.
"""
import sys
import hashlib,gzip
from pathlib import Path

import numpy as np
import pandas as pd


def content_hash(p):
    h=hashlib.sha256()
    opener=gzip.open if p.suffix=='.gz' else open
    with opener(p,'rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.digest()


def load(p):
    return pd.read_csv(p, sep="\t", low_memory=False)


def compare(a, b):
    if a.shape != b.shape or list(a.columns) != list(b.columns):
        return "differs", f"shape {a.shape} vs {b.shape}", np.nan
    for attempt in (0, 1):
        num = a.select_dtypes("number").columns.intersection(b.select_dtypes("number").columns)
        txt = [c for c in a.columns if c not in num]
        x, y = a[num].to_numpy(float), b[num].to_numpy(float)
        both_nan = np.isnan(x) & np.isnan(y)
        diff = np.where(both_nan, 0, np.abs(x - y))
        rel = np.where(both_nan, 0, diff / np.maximum(np.maximum(np.abs(x), np.abs(y)), 1e-12))
        nan_mismatch = int((np.isnan(x) ^ np.isnan(y)).sum())
        txt_mismatch = int((a[txt].astype(str).to_numpy() != b[txt].astype(str).to_numpy()).sum()) if txt else 0
        max_rel = float(np.nanmax(rel)) if rel.size else 0.0
        if max_rel <= 1e-9 and nan_mismatch == 0 and txt_mismatch == 0:
            return "identical", "", 0.0
        if max_rel < 1e-3 and nan_mismatch == 0 and txt_mismatch == 0:
            return "float noise", "", max_rel
        if attempt == 0:                                  # retry ignoring row order
            key = list(a.columns)
            a = a.sort_values(key, kind="mergesort").reset_index(drop=True)
            b = b.sort_values(key, kind="mergesort").reset_index(drop=True)
    col = num[int(np.nanargmax(rel.max(axis=0)))] if rel.size and np.isfinite(rel).any() else ""
    return "differs", f"text cells {txt_mismatch}, NaN mismatches {nan_mismatch}, largest change in '{col}'", max_rel


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    base = Path(args[0]); cur = Path(args[1]) if len(args) > 1 else Path(__file__).resolve().parents[1] / "results" / "tables"
    out = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv else None
    names = sorted({p.name for p in base.glob("*.tsv*")} | {p.name for p in cur.glob("*.tsv*")})
    rows = []
    for n in names:
        pa, pb = base / n, cur / n
        if not (pa.exists() and pb.exists()):
            rows.append(dict(table=n, status="missing", detail="baseline only" if pa.exists() else "current only", max_rel_diff=np.nan))
            continue
        if content_hash(pa)==content_hash(pb):
            status,detail,mr='identical','equal decompressed SHA-256',0.0
        else:
            status, detail, mr = compare(load(pa), load(pb))
        rows.append(dict(table=n, status=status, detail=detail, max_rel_diff=mr))
    res = pd.DataFrame(rows)
    print(res.status.value_counts().to_string())
    print(res[res.status.isin(["differs", "missing"])].to_string(index=False))
    if out:
        res.to_csv(out, sep="\t", index=False, float_format="%.3g")


if __name__ == "__main__":
    main()
