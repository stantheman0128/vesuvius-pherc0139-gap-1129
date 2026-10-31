#!/usr/bin/env bash
# 7-point sheet traces from the w035 seed at C240, thresholds 0.5 (checkpoint default) and 0.3 (docstring's best).
set -u
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 EXTRA_SITE=${EXTRA_SITE:-}
PY=${PY:-python3}
STEPS=${STEPS:-8}
for th in 0.5 0.3; do
  n="sheet7_s${STEPS}_t${th/./}"
  /usr/bin/time -v $PY -I scripts/run_sheet.py "out/${n}_record.json" -- \
    --checkpoint downloads/sheet_tracer_7pt_gap60/last.pt --volume out/ct_view_l3_as_l0 \
    --seed-patch out/inputs/w035_C240_seed_l3 --line row:22 --range 1338:1484 --direction + \
    --steps "$STEPS" --threshold "$th" --out "out/$n" --device cpu > "logs/$n.log" 2>&1
  echo "exit=$?" >> "logs/$n.log"
  $PY -I scripts/trace_to_npz.py "out/$n"
  grep -E "steps x|status|Elapsed|Maximum res" "logs/$n.log"
done
