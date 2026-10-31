#!/usr/bin/env bash
# CT ablation at C240 (views built by ablate_build.py). One process at a time (the copy model peaks at 17 to 21 GiB).
# Order: the runs the coordinator asked for first, then the opposite-sign shift, then (c) shift with TTA.
set -u
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 EXTRA_SITE=${EXTRA_SITE:-}
PY=${PY:-python3}
A=out/ablate

copy() {  # name volume_dir extra_args...
  local n=$1 vol=$2; shift 2
  echo "$(date +%T) start $n"
  /usr/bin/time -v $PY -I scripts/run_copy.py "$A/$n" -- \
    --tifxyz-path out/inputs/w035_C240_patch --volume-path "$vol" \
    --checkpoint-path out/ckpt/copy_displacement_latest_modelonly.pth --device cpu --volume-scale 2 --no-compile \
    --batch-size 1 --tta-batch-size 1 --bbox-band-workers 1 --crop-input-workers 1 --tifxyz-voxel-size-um 1.129 \
    --merge-temp-dir "$A/tmp_merge_$n" --out-dir "$A/$n" --verbose "$@" > "logs/abl_$n.log" 2>&1
  echo "exit=$?" >> "logs/abl_$n.log"
  grep -E "Elapsed|Maximum res|exit=" "logs/abl_$n.log"
}

sheet() {  # name view_dir
  local n=$1 vol=$2
  echo "$(date +%T) start $n"
  /usr/bin/time -v $PY -I scripts/run_sheet.py "$A/${n}_record.json" -- \
    --checkpoint downloads/sheet_tracer_7pt_gap60/last.pt --volume "$vol" \
    --seed-patch out/inputs/w035_C240_seed_l3 --line row:22 --range 1338:1484 --direction + \
    --steps 8 --threshold 0.3 --out "$A/$n" --device cpu > "logs/abl_$n.log" 2>&1
  echo "exit=$?" >> "logs/abl_$n.log"
  $PY -I scripts/trace_to_npz.py "$A/$n" >> "logs/abl_$n.log" 2>&1
  grep -E "steps x|status|Elapsed|Maximum res|exit=" "logs/abl_$n.log"
}

copy copy_const $A/ct_const_l2/vol.zarr --no-tta
copy copy_shift_p $A/ct_shift_l2_p/vol.zarr --no-tta
sheet sheet_const $A/sheet_const
sheet sheet_shift_p $A/sheet_shift_p
sheet sheet_shift_m $A/sheet_shift_m
copy copy_shift_m $A/ct_shift_l2_m/vol.zarr --no-tta
copy copy_shift_p_tta $A/ct_shift_l2_p/vol.zarr
echo "$(date +%T) ALL_DONE"
