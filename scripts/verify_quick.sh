#!/usr/bin/env bash
# Non-mutating package check: hashes; small-input verification; steps 11/30/23; all six figures.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="${PY:-python}"
"$PY" "$ROOT/scripts/verify_manifest.py" "$ROOT"
TMP="$(mktemp -d)"
echo "Isolated verification workspace: $TMP"
mkdir -p "$TMP/project/data/raw" "$TMP/project/results"
cp -R "$ROOT/scripts" "$TMP/project/scripts"
mkdir -p "$TMP/project/references"
cp -R "$ROOT/references/vendor" "$TMP/project/references/vendor"
cp -R "$ROOT/data/processed" "$TMP/project/data/processed"
cp -R "$ROOT/results/tables" "$TMP/project/results/tables"
cp "$ROOT/data/raw/SHA256SUMS.txt" "$TMP/project/data/raw/SHA256SUMS.txt"
cd "$TMP/project"
get() { mkdir -p "$(dirname "$1")"; curl -f -sS -L --retry 3 -o "$1" "$2"; }
for cohort in GSE89408 GSE198520; do
  series="${cohort:0:${#cohort}-3}nnn"
  get "data/raw/$cohort/${cohort}_series_matrix.txt.gz" "https://ftp.ncbi.nlm.nih.gov/geo/series/$series/$cohort/matrix/${cohort}_series_matrix.txt.gz"
done
get data/raw/GSE89408/GSE89408_GEO_count_matrix_rename.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE89nnn/GSE89408/suppl/GSE89408_GEO_count_matrix_rename.txt.gz
get data/raw/GSE198520/GSE198520_Raw_gene_count_matrix.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198520/suppl/GSE198520_Raw_gene_count_matrix.txt.gz
"$PY" - <<'CHECK'
from pathlib import Path
import hashlib
root=Path('data/raw')
checks={p.removeprefix('./'):s for s,p in (x.split(maxsplit=1) for x in (root/'SHA256SUMS.txt').read_text().splitlines())}
for p in list((root/'GSE89408').glob('*'))+list((root/'GSE198520').glob('*')):
    name=str(p.relative_to(root));assert name in checks,'Missing input checksum: '+name
    assert hashlib.sha256(p.read_bytes()).hexdigest()==checks[name],'Input checksum differs: '+name
print('PASS: four small raw inputs, including metadata')
CHECK
for s in 11_gse198520_response 30_gse89408_disease_stage 23_knowledge_graph; do
  (cd scripts && "$PY" -u "$s.py" > "$TMP/$s.log" 2>&1) || { tail -20 "$TMP/$s.log"; exit 1; }
  echo "$s: completed"
done
"$PY" scripts/compare_results.py "$ROOT/results/tables" results/tables --out "$TMP/check_tables.tsv"
"$PY" scripts/compare_results.py "$ROOT/data/processed" data/processed --out "$TMP/check_processed.tsv"
"$PY" - "$TMP" <<'CHECK'
import sys,pandas as pd
from pathlib import Path
for name in ['check_tables.tsv','check_processed.tsv']:
    d=pd.read_csv(Path(sys.argv[1])/name,sep='\t')
    assert not d.status.isin(['differs','missing']).any(),name+' has unresolved differences'
print('PASS: isolated table comparisons')
CHECK
"$PY" scripts/39_main_figures.py > "$TMP/figures.log" 2>&1 || { tail -20 "$TMP/figures.log"; exit 1; }
"$PY" - <<'CHECK'
from pathlib import Path
for name in ['Fig1','Fig2','Fig3','Fig4','ExtData_Fig1','ExtData_Fig2']:
    for ext in ['pdf','svg','png']:
        p=Path('results/figures/main')/f'{name}.{ext}'
        assert p.is_file() and p.stat().st_size>1000,str(p)
print('PASS: all six figures rendered without raw single-cell or interim objects')
CHECK
"$PY" "$ROOT/scripts/verify_manifest.py" "$ROOT"
echo "PASS: package unchanged; logs and comparison evidence retained at $TMP"
