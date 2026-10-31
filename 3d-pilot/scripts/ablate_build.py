"""Build the CT-ablation volumes for C240 from the CT chunks already downloaded (no new downloads).

(a) const: every voxel of the box set to the box's median value.
(b) shift: CT content translated by an integer number of voxels along the zyx axis closest to w035's mean normal at
    C240, so that shifted[p] = original[p - t]. Copy model: level 2, 7 voxels (31.6 um = 28 level-0 voxels). Sheet
    tracer: level 3 served as its level 0, 3 voxels (27.1 um; 3.5 would be 31.6 um). Both signs are built.
Voxels whose source lies outside the downloaded box are filled with the median; the script checks that no read
region of the tools touches them.
Usage: python -I ablate_build.py      Writes out/ablate/<view>/... and out/ablate/build.json
"""
import json
import os
import shutil
import site
import sys

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
import tifffile  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "shim"))
from netlog import ROOT  # noqa: E402
import vc  # noqa: E402  (the numpy reader, used here as an independent read-back)

OUT = os.path.join(ROOT, "out", "ablate")
SRC = os.path.join(ROOT, "downloads", "ct", "vol.zarr")
UM0 = 1.129
CH = 128


def unit(v):
    return v / np.linalg.norm(v)


# --- w035 normal at C240, oriented outward (away from the umbilicus, toward w036) ---
pd = os.path.join(ROOT, "out", "inputs", "w035_C240_patch")
g = np.stack([tifffile.imread(os.path.join(pd, f"{c}.tif")).astype(np.float64) for c in "xyz"], -1)
n = np.cross(np.gradient(g, axis=1), np.gradient(g, axis=0))
n /= np.linalg.norm(n, axis=-1, keepdims=True)
umb = np.array(json.load(open(os.path.join(ROOT, "out", "locate.json")))["umb11_xyz"], np.float64)
radial = g[..., :2] - umb[None, None, :2]
if np.mean(np.sum(n[..., :2] * radial, -1)) < 0:
    n = -n
n_mean = unit(n.reshape(-1, 3).mean(0))
ax_xyz = int(np.argmax(np.abs(n_mean)))
sgn = int(np.sign(n_mean[ax_xyz]))
ax_zyx = 2 - ax_xyz
info = {"n_mean_xyz_outward": n_mean.round(4).tolist(), "axis_xyz": "xyz"[ax_xyz], "axis_zyx_index": ax_zyx,
        "outward_sign_on_axis": sgn,
        "angle_normal_to_axis_deg": round(float(np.degrees(np.arccos(abs(n_mean[ax_xyz])))), 1)}


def load_level(level):
    meta = json.load(open(os.path.join(SRC, str(level), ".zarray")))
    keys = []
    for z in os.listdir(os.path.join(SRC, str(level))):
        if not z.isdigit():
            continue
        for y in os.listdir(os.path.join(SRC, str(level), z)):
            for x in os.listdir(os.path.join(SRC, str(level), z, y)):
                keys.append((int(z), int(y), int(x)))
    k = np.array(keys)
    lo, hi = k.min(0), k.max(0)
    if len(keys) != int(np.prod(hi - lo + 1)):
        raise SystemExit(f"level {level}: downloaded chunks are not a full box")
    arr = np.empty(tuple((hi - lo + 1) * CH), np.uint8)
    for (z, y, x) in keys:
        b = np.fromfile(os.path.join(SRC, str(level), str(z), str(y), str(x)), np.uint8).reshape(CH, CH, CH)
        o = (np.array((z, y, x)) - lo) * CH
        arr[o[0]:o[0] + CH, o[1]:o[1] + CH, o[2]:o[2] + CH] = b
    return meta, lo, arr, keys


def write_view(dst_level_dir, meta, lo, arr, keys):
    if os.path.exists(dst_level_dir):
        shutil.rmtree(dst_level_dir)
    os.makedirs(dst_level_dir)
    json.dump(meta, open(os.path.join(dst_level_dir, ".zarray"), "w"))
    for (z, y, x) in keys:
        o = (np.array((z, y, x)) - lo) * CH
        d = os.path.join(dst_level_dir, str(z), str(y))
        os.makedirs(d, exist_ok=True)
        np.ascontiguousarray(arr[o[0]:o[0] + CH, o[1]:o[1] + CH, o[2]:o[2] + CH]).tofile(os.path.join(d, str(x)))


def shifted(arr, d, axis, fill):
    out = np.full_like(arr, fill)
    src = [slice(None)] * 3
    dst = [slice(None)] * 3
    if d > 0:
        dst[axis], src[axis] = slice(d, None), slice(None, -d)
    else:
        dst[axis], src[axis] = slice(None, d), slice(-d, None)
    out[tuple(dst)] = arr[tuple(src)]
    return out


def check_region(name, lo_box, hi_box, lo_chunk, shape, d, axis):
    """Read region [lo_box, hi_box) in level voxels must avoid the |d|-voxel fill strip of the shifted box."""
    rel_lo = np.array(lo_box) - lo_chunk * CH
    rel_hi = np.array(hi_box) - lo_chunk * CH
    strip = (0, d) if d > 0 else (shape[axis] + d, shape[axis])
    ok = bool(rel_hi[axis] <= strip[0] or rel_lo[axis] >= strip[1]) and bool(np.all(rel_lo >= 0)) and bool(
        np.all(rel_hi <= np.array(shape)))
    return {"region": name, "rel_lo": rel_lo.tolist(), "rel_hi": rel_hi.tolist(), "fill_strip_on_axis": list(strip),
            "avoids_fill_and_inside_box": ok}


views = {}
# --- copy model, level 2 ---
meta2, lo2, a2, keys2 = load_level(2)
plan = json.load(open(os.path.join(ROOT, "out", "plan_ct_l2.json")))
(cmn, cmx) = plan["boxes"][0]
rel = np.array(cmn) - lo2 * CH
crop = a2[rel[0]:rel[0] + 128, rel[1]:rel[1] + 384, rel[2]:rel[2] + 384]
assert crop.shape == (128, 384, 384)
med2 = int(np.median(crop))
info["copy_l2"] = {"crop_zyx_min": cmn, "crop_zyx_max": cmx, "median": med2, "box_chunks": len(keys2)}
write_view(os.path.join(OUT, "ct_const_l2", "vol.zarr", "2"), meta2, lo2, np.full_like(a2, med2), keys2)
views["ct_const_l2"] = None
for s, tag in ((1, "p"), (-1, "m")):
    d = 7 * sgn * s
    write_view(os.path.join(OUT, f"ct_shift_l2_{tag}", "vol.zarr", "2"), meta2, lo2, shifted(a2, d, ax_zyx, med2), keys2)
    t_zyx = [0, 0, 0]
    t_zyx[ax_zyx] = d
    t_um_xyz = (np.array(t_zyx[::-1], np.float64) * 4 * UM0)
    info[f"copy_shift_{tag}"] = {"d_level2_vox": d, "t_zyx_level2": t_zyx, "t_um_xyz": t_um_xyz.round(2).tolist(),
                                 "t_dot_n_mean_um": round(float(t_um_xyz @ n_mean), 2),
                                 "t_dot_n_per_vertex_mean_um": round(float((n.reshape(-1, 3) @ t_um_xyz).mean()), 2),
                                 "check": check_region("copy crop", cmn, cmx, lo2, a2.shape, d, ax_zyx)}
    # independent read-back through the shim: shifted crop == original crop read at (crop - t)
    va, vb = vc.Volume.open(os.path.join(OUT, f"ct_shift_l2_{tag}", "vol.zarr")), vc.Volume.open(SRC)
    rb = va.read_zyx(cmn, [128, 384, 384], level=2)
    ref = vb.read_zyx((np.array(cmn) - np.array(t_zyx)).tolist(), [128, 384, 384], level=2)
    info[f"copy_shift_{tag}"]["readback_equal"] = bool(np.array_equal(rb, ref))
    info[f"copy_shift_{tag}"]["readback_differs_from_unshifted"] = bool(
        not np.array_equal(rb, vb.read_zyx(cmn, [128, 384, 384], level=2)))
va = vc.Volume.open(os.path.join(OUT, "ct_const_l2", "vol.zarr"))
rb = va.read_zyx(cmn, [128, 384, 384], level=2)
info["copy_const_readback_unique"] = np.unique(rb).tolist()

# --- sheet tracer, level 3 served as level 0 ---
meta3, lo3, a3, keys3 = load_level(3)
sp = json.load(open(os.path.join(ROOT, "out", "inputs", "sheet_plan.json")))
rlo = np.floor(np.array(sp["read_lo_xyz_l3"])[::-1]).astype(int)
rhi = np.ceil(np.array(sp["read_hi_xyz_l3"])[::-1]).astype(int) + 1
r0 = rlo - lo3 * CH
r1 = rhi - lo3 * CH
med3 = int(np.median(a3[r0[0]:r1[0], r0[1]:r1[1], r0[2]:r1[2]]))
info["sheet_l3"] = {"read_zyx_lo": rlo.tolist(), "read_zyx_hi": rhi.tolist(), "median": med3, "box_chunks": len(keys3)}
write_view(os.path.join(OUT, "sheet_const", "0"), meta3, lo3, np.full_like(a3, med3), keys3)
for s, tag in ((1, "p"), (-1, "m")):
    d = 3 * sgn * s
    write_view(os.path.join(OUT, f"sheet_shift_{tag}", "0"), meta3, lo3, shifted(a3, d, ax_zyx, med3), keys3)
    t_zyx = [0, 0, 0]
    t_zyx[ax_zyx] = d
    t_um_xyz = np.array(t_zyx[::-1], np.float64) * 8 * UM0
    info[f"sheet_shift_{tag}"] = {"d_level3_vox": d, "t_um_xyz": t_um_xyz.round(2).tolist(),
                                  "t_dot_n_mean_um": round(float(t_um_xyz @ n_mean), 2),
                                  "check": check_region("sheet plan read box", rlo, rhi, lo3, a3.shape, d, ax_zyx)}
os.makedirs(OUT, exist_ok=True)
json.dump(info, open(os.path.join(OUT, "build.json"), "w"), indent=1)
print(json.dumps(info, indent=1))
