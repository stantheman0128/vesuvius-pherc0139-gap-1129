"""Displacement of the ablation outputs relative to the real-CT runs, along w035's outward normal.

Copy model: vertex-by-vertex on the shared 25 x 25 grid (front and back), normal = input w035 patch normal at the same
vertex. Sheet tracer: point-by-point at the same (step, column) where both runs kept the point (seed row excluded),
normal = w035 band normal at the nearest band vertex. Outward = away from the umbilicus (toward w036).
For the shifts, each point's translation component along its normal (t . n) is reported beside the measured mean.
Usage: python -I ablate_compare.py      Writes out/ablate/compare.json
"""
import json
import os
import site
import sys

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
import tifffile  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import ROOT  # noqa: E402

UM0 = 1.129
OUT = os.path.join(ROOT, "out")
A = os.path.join(OUT, "ablate")
build = json.load(open(os.path.join(A, "build.json")))
if os.path.exists(os.path.join(A, "build_small.json")):
    build.update(json.load(open(os.path.join(A, "build_small.json"))))
umb = np.array(json.load(open(os.path.join(OUT, "locate.json")))["umb11_xyz"], np.float64)


def grid(d):
    g = np.stack([tifffile.imread(os.path.join(d, f"{c}.tif")).astype(np.float64) for c in "xyz"], -1)
    return g, np.all(g != -1, -1) & np.all(np.isfinite(g), -1)


def outward_normals(g):
    n = np.cross(np.gradient(g, axis=1), np.gradient(g, axis=0))
    n /= np.linalg.norm(n, axis=-1, keepdims=True) + 1e-12
    flip = np.sum(n[..., :2] * (g[..., :2] - umb[:2]), -1) < 0
    n[flip] *= -1
    return n


def summ(v):
    v = np.asarray(v, np.float64)
    if v.size == 0:
        return {"n": 0}
    return {"n": int(v.size), "mean": round(float(v.mean()), 2), "median": round(float(np.median(v)), 2),
            "p10": round(float(np.percentile(v, 10)), 2), "p90": round(float(np.percentile(v, 90)), 2)}


res = {}
pg, pv = grid(os.path.join(OUT, "inputs", "w035_C240_patch"))
pn = outward_normals(pg)
copy_runs = {"copy_const": ("copy_notta", None), "copy_shift_p": ("copy_notta", "copy_shift_p"),
             "copy_shift_m": ("copy_notta", "copy_shift_m"), "copy_shift_p_tta": ("copy_tta", "copy_shift_p"),
             "copy_shift_p2": ("copy_notta", "copy_shift_p2"), "copy_shift_m2": ("copy_notta", "copy_shift_m2"),
             "copy_real_repeat": ("copy_notta", None)}
for name, (ref, tkey) in copy_runs.items():
    d = os.path.join(A, name)
    if not os.path.isdir(os.path.join(d, "w035_C240_patch_front")):
        res[name] = "not run"
        continue
    r = {"reference_run": ref}
    t = np.array(build[tkey]["t_um_xyz"]) if tkey else None
    for side in ("front", "back"):
        ga, va = grid(os.path.join(d, f"w035_C240_patch_{side}"))
        gr, vr = grid(os.path.join(OUT, ref, f"w035_C240_patch_{side}"))
        ok = va & vr & pv
        disp = (ga - gr)[ok] * UM0
        nn = pn[ok]
        r[side] = {"valid_both": int(ok.sum()), "signed_along_normal_um": summ(np.sum(disp * nn, -1)),
                   "abs_disp_um": summ(np.linalg.norm(disp, axis=-1))}
        if t is not None:
            tn = nn @ t
            r[side]["t_dot_n_um"] = summ(tn)
            r[side]["along_t_direction_um"] = summ(disp @ (t / np.linalg.norm(t)))
            r[side]["ratio_mean_signed_to_mean_t_dot_n"] = round(float(np.sum(disp * nn, -1).mean() / tn.mean()), 3)
    res[name] = r

# sheet tracer
band = np.load(os.path.join(ROOT, "downloads", "mesh_bands", "w035_band.npz"))["g"].astype(np.float64)
bv = np.all(band != -1, -1) & np.all(np.isfinite(band), -1)
bn = outward_normals(band)
tree = cKDTree(band[bv])
bn_valid = bn[bv]
ref = np.load(os.path.join(OUT, "sheet7_s8_t03", "trace.npz"))
for name, tkey in (("trace_const", None), ("trace_shift_p", "sheet_shift_p"), ("trace_shift_m", "sheet_shift_m")):
    p = os.path.join(A, name, "trace.npz")
    if not os.path.exists(p):
        res[name] = "not run"
        continue
    a = np.load(p)
    rows = min(a["x"].shape[0], ref["x"].shape[0])
    ga = np.stack([a[c][:rows] for c in "xyz"], -1).astype(np.float64) * 8
    gr = np.stack([ref[c][:rows] for c in "xyz"], -1).astype(np.float64) * 8
    ok = np.all(ga >= 0, -1) & np.all(gr >= 0, -1)
    ok[0] = False
    disp = (ga - gr)[ok] * UM0
    _, idx = tree.query(gr[ok])
    nn = bn_valid[idx]
    r = {"rows_compared": rows - 1, "points_both_kept": int(ok.sum()), "rows_in_run": int(a["x"].shape[0]),
         "signed_along_normal_um": summ(np.sum(disp * nn, -1)), "abs_disp_um": summ(np.linalg.norm(disp, axis=-1)),
         "per_row_mean_signed_um": [round(float(np.sum(((ga - gr)[i][ok[i]] * UM0) * bn_valid[tree.query(gr[i][ok[i]])[1]], -1).mean()), 1)
                                    if ok[i].any() else None for i in range(1, rows)]}
    if tkey:
        t = np.array(build[tkey]["t_um_xyz"])
        r["t_dot_n_um"] = summ(nn @ t)
        r["ratio_mean_signed_to_mean_t_dot_n"] = round(float(np.sum(disp * nn, -1).mean() / (nn @ t).mean()), 3)
    res[name] = r
json.dump(res, open(os.path.join(A, "compare.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
