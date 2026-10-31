"""Compare a sheet-tracer trace seeded on w035 at C240 with the published w034/w035/w036.

The trace is in level-3 voxels (x8 = 1.129 um level-0 voxels). Per traced row (step) it reports the median
distance of the kept points to each mesh and which mesh is nearest; a trace that stays on w035's sheet should be
nearest w035 on every row, with w034 about 110 um and w036 about 130 um away.
Usage: python -I analyze_sheet.py <trace_dir> <tag>
"""
import json
import os
import site
import sys

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
import tifffile  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import ROOT  # noqa: E402
from meshgeom import densify, cut  # noqa: E402

UM = 1.129
trace_dir, tag = sys.argv[1], sys.argv[2]
OUT = os.path.join(ROOT, "out")
PNG = os.path.join(OUT, "png")

_t = np.load(os.path.join(trace_dir, "trace.npz"))  # from trace_to_npz.py (tifs are LZW)
g = np.stack([_t[c].astype(np.float64) for c in "xyz"], -1)
keep = np.all(g != -1, -1)
g0 = np.where(keep[..., None], g * 8.0, np.nan)   # level-3 -> level-0 voxels
conf = _t["conf"]
trees = {}
for w in ("w034", "w035", "w036"):
    b = np.load(os.path.join(ROOT, "downloads", "mesh_bands", f"{w}_band.npz"))["g"].astype(np.float64)
    trees[w] = cKDTree(densify(b, np.all(b != -1, -1)))
rows = []
for r in range(g.shape[0]):
    m = keep[r]
    if not m.any():
        rows.append({"row": r, "kept": 0})
        continue
    Q = g0[r][m]
    dd = {w: trees[w].query(Q)[0] * UM for w in trees}
    near = np.argmin(np.stack([dd[w] for w in ("w034", "w035", "w036")]), 0)
    rows.append({"row": r, "kept": int(m.sum()), "of": int(m.size), "z_l0_median": round(float(np.median(Q[:, 2])), 1),
                 **{f"med_to_{w}_um": round(float(np.median(dd[w])), 1) for w in dd},
                 **{f"max_to_{w}_um": round(float(np.max(dd[w])), 1) for w in ("w035",)},
                 "nearest": {w: int((near == i).sum()) for i, w in enumerate(("w034", "w035", "w036"))}})
allQ = g0[keep]
dd_all = {w: trees[w].query(allQ)[0] * UM for w in trees}
near_all = np.argmin(np.stack([dd_all[w] for w in ("w034", "w035", "w036")]), 0)
res = {"tag": tag, "shape_rows_cols": list(g.shape[:2]), "kept_points": int(keep.sum()), "of": int(keep.size),
       "kept_after_seed_row": int(keep[1:].sum()), "of_after_seed_row": int(keep[1:].size),
       "median_to_um": {w: round(float(np.median(dd_all[w])), 1) for w in dd_all},
       "p90_to_um": {w: round(float(np.percentile(dd_all[w], 90)), 1) for w in dd_all},
       "nearest_counts": {w: int((near_all == i).sum()) for i, w in enumerate(("w034", "w035", "w036"))},
       "conf_mean_by_row": [round(float(c), 3) for c in conf.mean(1)], "rows": rows}
json.dump(res, open(os.path.join(OUT, f"{tag}_metrics.json"), "wt"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "rows"}, indent=1))
for row in rows:
    print(row)

# figure: CT level-3 slice through the trace's middle row height, with mesh cuts and traced points near it
zmid = np.nanmedian(g0[g.shape[0] // 2][:, 2]) if keep[g.shape[0] // 2].any() else np.nanmedian(allQ[:, 2])
zc3 = int(round(zmid / 8.0))
plan = json.load(open(os.path.join(OUT, "inputs", "sheet_plan.json")))
lo, hi = np.array(plan["read_lo_xyz_l3"]).astype(int), np.array(plan["read_hi_xyz_l3"]).astype(int)
img = np.zeros((hi[1] - lo[1], hi[0] - lo[0]), np.uint8)
root = os.path.join(ROOT, "downloads", "ct", "vol.zarr", "3")
for cy in range(lo[1] // 128, (hi[1] - 1) // 128 + 1):
    for cx in range(lo[0] // 128, (hi[0] - 1) // 128 + 1):
        fn = os.path.join(root, str(zc3 // 128), str(cy), str(cx))
        if not os.path.exists(fn):
            continue
        a = np.memmap(fn, dtype=np.uint8, mode="r", shape=(128, 128, 128))[zc3 % 128]
        ys, xs = max(lo[1], cy * 128), max(lo[0], cx * 128)
        ye, xe = min(hi[1], cy * 128 + 128), min(hi[0], cx * 128 + 128)
        img[ys - lo[1]:ye - lo[1], xs - lo[0]:xe - lo[0]] = a[ys - cy * 128:ye - cy * 128, xs - cx * 128:xe - cx * 128]
f = 8 * UM / 1000
ext = [lo[0] * f, hi[0] * f, hi[1] * f, lo[1] * f]
fig, axs = plt.subplots(1, 2, figsize=(15, 6.5))
cols = {"w034": "#1f77b4", "w035": "#2ca02c", "w036": "#d62728"}

for k, ax in enumerate(axs):
    ax.imshow(img, cmap="gray", extent=ext, vmin=0, vmax=255)
    for w in cols:
        b = np.load(os.path.join(ROOT, "downloads", "mesh_bands", f"{w}_band.npz"))["g"].astype(np.float64)
        c = cut(b, np.all(b != -1, -1), zc3 * 8.0)
        ax.plot(c[:, 0] * UM / 1000, c[:, 1] * UM / 1000, ".", ms=3, color=cols[w], label=f"{w} (published)")
    gz = np.where(keep[..., None], g0, -1.0)
    c = cut(gz, keep, zc3 * 8.0)
    ax.plot(c[:, 0] * UM / 1000, c[:, 1] * UM / 1000, "x", ms=7, mew=1.5, color="#ff7f0e", label="sheet trace (kept)")
    ax.set_xlabel("x (mm, 1.129 um frame)")
    ax.set_ylabel("y (mm)")
axs[0].set_xlim(ext[0], ext[1])
axs[0].set_ylim(ext[2], ext[3])
cx_, cy_ = np.nanmedian(allQ[:, 0]) * UM / 1000, np.nanmedian(allQ[:, 1]) * UM / 1000
axs[1].set_xlim(cx_ - 0.6, cx_ + 0.6)
axs[1].set_ylim(cy_ + 0.4, cy_ - 0.4)
axs[0].legend(loc="lower left", fontsize=8)
axs[0].set_title(f"sheet tracer seeded on w035 at C240; slice z = {zc3 * 8} (1.129 um frame), CT level 3 (9.0 um)")
axs[1].set_title("zoom")
fig.tight_layout()
fig.savefig(os.path.join(PNG, f"{tag}_slice.png"), dpi=110)
print("figure", os.path.join(PNG, f"{tag}_slice.png"))
