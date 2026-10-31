#!/usr/bin/env bash
# Sheet-tracer part of the CT ablation. run_ablate.sh wrote the traces into the CT view directories themselves
# (out/ablate/sheet_*), which the tracer refused (FileExistsError, nothing written); traces go to out/ablate/trace_*.
set -u
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 EXTRA_SITE=${EXTRA_SITE:-}
PY=${PY:-python3}
A=out/ablate
for v in const shift_p shift_m; do
  n=trace_$v
  echo "$(date +%T) start $n"
  /usr/bin/time -v $PY -I scripts/run_sheet.py "$A/${n}_record.json" -- \
    --checkpoint downloads/sheet_tracer_7pt_gap60/last.pt --volume "$A/sheet_$v" \
    --seed-patch out/inputs/w035_C240_seed_l3 --line row:22 --range 1338:1484 --direction + \
    --steps 8 --threshold 0.3 --out "$A/$n" --device cpu > "logs/abl_$n.log" 2>&1
  echo "exit=$?" >> "logs/abl_$n.log"
  $PY -I scripts/trace_to_npz.py "$A/$n" >> "logs/abl_$n.log" 2>&1
  grep -E "steps x|status|Elapsed|Maximum res|exit=" "logs/abl_$n.log"
done
echo "$(date +%T) SHEET_DONE"
