#!/usr/bin/env bash
# Calibration runs: sheet tracer at each place of out/calib/plans.json on the real, +1 and -1 voxel CT views.
set -u
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 EXTRA_SITE=${EXTRA_SITE:-}
PY=${PY:-python3}
$PY -I -c "import json; [print(p['name'], p['line'], p['range'], p['direction']) for p in json.load(open('out/calib/plans.json'))]" > out/calib/args.txt
while read -r name line range dir; do
  for v in real p m; do
    n=out/calib/$name/trace_$v
    L=logs/calib_${name}_$v.log
    echo "$(date +%T) start $name $v"
    /usr/bin/time -v $PY -I scripts/run_sheet_calib.py "${n}_record.json" -- \
      --checkpoint downloads/sheet_tracer_7pt_gap60/last.pt --volume "out/calib/$name/view_$v" \
      --seed-patch "out/calib/$name/seed_l3" --line "$line" --range "$range" --direction "$dir" \
      --steps 8 --threshold 0.3 --out "$n" --device cpu > "$L" 2>&1 < /dev/null
    echo "exit=$?" >> "$L"
    $PY -I scripts/trace_to_npz.py "$n" >> "$L" 2>&1 < /dev/null
    grep -E "steps x|status|Elapsed|Maximum res|exit=|chunks read" "$L"
  done
done < out/calib/args.txt
echo "$(date +%T) CALIB_DONE"
