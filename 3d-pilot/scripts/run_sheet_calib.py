"""run_sheet.py plus two records for the calibration runs: every read_block call (origin, shape) so the analysis
can check that no read touched a shift-fill voxel, and the raw rows/confidences returned by trace() (all columns,
before connected_keep masks them) so the stop reason can be read off exactly (inference.py:148 vs :152).
Usage (WSL venv): python -I run_sheet_calib.py <record_json> -- <inference.py args...>
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_sheet as R  # noqa: E402  (villa_env, 4 threads, weights_only, chunk read tracking)
import numpy as np  # noqa: E402

T = R.T
READS, RAW = [], {}
_orig_rb, _orig_trace = T.read_block, T.trace


def _rb(volume, pos, pad, level):
    blk, origin = _orig_rb(volume, pos, pad, level)
    READS.append({"origin_zyx": [int(v) for v in origin], "shape_zyx": [int(v) for v in blk.shape],
                  "pad": float(pad)})
    return blk, origin


def _tr(*a, **k):
    rows, confs = _orig_trace(*a, **k)
    RAW["rows"], RAW["confs"] = rows, confs
    return rows, confs


T.read_block, T.trace = _rb, _tr

if __name__ == "__main__":
    base = sys.argv[1].replace("_record.json", "")
    try:
        R.main()
    finally:
        json.dump(READS, open(base + "_reads.json", "w"))
        if RAW:
            np.savez(base + "_raw.npz", rows=RAW["rows"], confs=RAW["confs"])
