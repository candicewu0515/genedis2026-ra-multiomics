#!/usr/bin/env bash
# Run the whole pipeline in dependency order (or a subset: bash scripts/run_all.sh 28 33 36; "22:pass2" labels a re-run).
# Prerequisites: `bash scripts/00_download.sh` done; environment `genedis_ra` active (or PY=/path/to/python);
# the LDSC environment `genedis_ldsc` next to it (or LDSC_PYTHON=/path/to/python). One log per step in logs/.
#
# Order notes: step 22 runs twice. Its first pass defines the genes with genetic evidence that step 25 tests, its
# second pass adds the colocalisation and fine-mapping layers of steps 25, 29 and 32; step 32 is re-run after the
# second pass only to refresh the prioritisation columns of its direction table; steps 23 and 39 use everything.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PY:-python}"
mkdir -p ../logs
ALL=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 24 25 26 27 28 29 30 31 32 22 32 33 34 35 36 37 38 40 41 23 39)
STEPS=("$@")
[ ${#STEPS[@]} -eq 0 ] && STEPS=("${ALL[@]}")
[ -f ../logs/run_all_timing.tsv ] || printf "step\tstarted\tseconds\tstatus\n" > ../logs/run_all_timing.tsv
seen=" "
for arg in "${STEPS[@]}"; do
  n="${arg%%:*}"                                   # "22:pass2" runs step 22 and labels its log as a second pass
  script=$(ls ${n}_*.py | head -1)
  name="${script%.py}"
  label="$name"
  [[ "$arg" == *:* ]] && label="$name.${arg#*:}"
  [[ "$arg" != *:* && "$seen" == *" $n "* ]] && label="$name.pass2"
  seen="$seen$n "
  start=$(date +%s); stamp=$(date "+%Y-%m-%d %H:%M:%S")
  echo "[$stamp] $label"
  if "$PY" -u "$script" > "../logs/$label.log" 2>&1; then status=ok; else status=FAILED; fi
  printf "%s\t%s\t%s\t%s\n" "$label" "$stamp" "$(( $(date +%s) - start ))" "$status" >> ../logs/run_all_timing.tsv
  if [ "$status" = FAILED ]; then
    echo "  FAILED; last lines of logs/$label.log:"; tail -20 "../logs/$label.log"; exit 1
  fi
done
echo "Done."
