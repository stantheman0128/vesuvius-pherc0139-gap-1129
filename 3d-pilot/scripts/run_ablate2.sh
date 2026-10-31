#!/usr/bin/env bash
# Second ablation batch, started after run_ablate.sh prints ALL_DONE (one tool process at a time):
# sheet-tracer runs (run_ablate_sheet.sh), then the small-shift views and two copy runs without TTA (+/-2 level-2
# voxels = +/-8 level-0 voxels = 9.03 um).
set -u
cd "$(dirname "$0")/.."
until grep -q ALL_DONE logs/abl_driver.log; do sleep 5; done
bash scripts/run_ablate_sheet.sh
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=${PY:-python3}
A=out/ablate
$PY -I scripts/ablate_build_small.py > logs/abl_build_small.log 2>&1
echo "build_small exit=$?"
for tag in p2 m2; do
  n=copy_shift_$tag
  echo "$(date +%T) start $n"
  /usr/bin/time -v $PY -I scripts/run_copy.py "$A/$n" -- \
    --tifxyz-path out/inputs/w035_C240_patch --volume-path "$A/ct_shift_l2_$tag/vol.zarr" \
    --checkpoint-path out/ckpt/copy_displacement_latest_modelonly.pth --device cpu --volume-scale 2 --no-compile \
    --batch-size 1 --tta-batch-size 1 --bbox-band-workers 1 --crop-input-workers 1 --tifxyz-voxel-size-um 1.129 \
    --merge-temp-dir "$A/tmp_merge_$n" --out-dir "$A/$n" --verbose --no-tta > "logs/abl_$n.log" 2>&1
  echo "exit=$?" >> "logs/abl_$n.log"
  grep -E "Elapsed|Maximum res|exit=" "logs/abl_$n.log"
done
echo "$(date +%T) ALL2_DONE"
