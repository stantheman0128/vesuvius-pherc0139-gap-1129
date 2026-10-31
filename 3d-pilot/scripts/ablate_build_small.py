"""Small-shift CT views for the copy model at C240: level 2 shifted by +/-2 voxels (9.03 um = 8 level-0 voxels)
along the same axis as ablate_build.py (read from out/ablate/build.json). Builds only these two views, so the views
in use by other runs are not touched. numpy only (runs in the WSL venv).
Usage: python -I ablate_build_small.py      Writes out/ablate/ct_shift_l2_{p2,m2}/ and out/ablate/build_small.json
"""
import json
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "shim"))
from netlog import ROOT  # noqa: E402
import vc  # noqa: E402

A = os.path.join(ROOT, "out", "ablate")
SRC = os.path.join(ROOT, "downloads", "ct", "vol.zarr")
CH = 128
b = json.load(open(os.path.join(A, "build.json")))
ax, sgn = b["axis_zyx_index"], b["outward_sign_on_axis"]
n_mean = np.array(b["n_mean_xyz_outward"])
cmn, cmx = b["copy_l2"]["crop_zyx_min"], b["copy_l2"]["crop_zyx_max"]
med = b["copy_l2"]["median"]

meta = json.load(open(os.path.join(SRC, "2", ".zarray")))
keys = [(int(z), int(y), int(x)) for z in os.listdir(os.path.join(SRC, "2")) if z.isdigit()
        for y in os.listdir(os.path.join(SRC, "2", z)) for x in os.listdir(os.path.join(SRC, "2", z, y))]
k = np.array(keys)
lo, hi = k.min(0), k.max(0)
assert len(keys) == int(np.prod(hi - lo + 1)) == b["copy_l2"]["box_chunks"]
arr = np.empty(tuple((hi - lo + 1) * CH), np.uint8)
for (z, y, x) in keys:
    o = (np.array((z, y, x)) - lo) * CH
    arr[o[0]:o[0] + CH, o[1]:o[1] + CH, o[2]:o[2] + CH] = np.fromfile(
        os.path.join(SRC, "2", str(z), str(y), str(x)), np.uint8).reshape(CH, CH, CH)

info = {}
for s, tag in ((1, "p2"), (-1, "m2")):
    d = 2 * sgn * s
    out = np.full_like(arr, med)
    src, dst = [slice(None)] * 3, [slice(None)] * 3
    if d > 0:
        dst[ax], src[ax] = slice(d, None), slice(None, -d)
    else:
        dst[ax], src[ax] = slice(None, d), slice(-d, None)
    out[tuple(dst)] = arr[tuple(src)]
    lvl = os.path.join(A, f"ct_shift_l2_{tag}", "vol.zarr", "2")
    if os.path.exists(lvl):
        shutil.rmtree(lvl)
    os.makedirs(lvl)
    json.dump(meta, open(os.path.join(lvl, ".zarray"), "w"))
    for (z, y, x) in keys:
        o = (np.array((z, y, x)) - lo) * CH
        dd = os.path.join(lvl, str(z), str(y))
        os.makedirs(dd, exist_ok=True)
        np.ascontiguousarray(out[o[0]:o[0] + CH, o[1]:o[1] + CH, o[2]:o[2] + CH]).tofile(os.path.join(dd, str(x)))
    t_zyx = [0, 0, 0]
    t_zyx[ax] = d
    t_um = np.array(t_zyx[::-1], np.float64) * 4 * 1.129
    va, vb = vc.Volume.open(os.path.dirname(lvl)), vc.Volume.open(SRC)
    rb = va.read_zyx(cmn, [128, 384, 384], level=2)
    ref = vb.read_zyx((np.array(cmn) - np.array(t_zyx)).tolist(), [128, 384, 384], level=2)
    rel_lo = np.array(cmn) - lo * CH
    rel_hi = np.array(cmx) - lo * CH
    strip = (0, d) if d > 0 else (arr.shape[ax] + d, arr.shape[ax])
    info[f"copy_shift_{tag}"] = {"d_level2_vox": d, "t_zyx_level2": t_zyx, "t_um_xyz": t_um.round(2).tolist(),
                                 "t_dot_n_mean_um": round(float(t_um @ n_mean), 2),
                                 "readback_equal": bool(np.array_equal(rb, ref)),
                                 "crop_avoids_fill": bool(rel_hi[ax] <= strip[0] or rel_lo[ax] >= strip[1])}
json.dump(info, open(os.path.join(A, "build_small.json"), "w"), indent=1)
print(json.dumps(info, indent=1))
