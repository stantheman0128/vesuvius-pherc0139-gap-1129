"""Fetch row bands of the on-1.129 tifxyz copies of w034/w035/w036 around C240 and build the inputs.

- downloads/mesh_bands/<w>_band.npz: rows [r-B, r+B] (all columns) of x, y, z, plus the row offset.
- out/inputs/w035_C240_patch/: a (2P+1) x (2P+1) tifxyz crop of w035 centred on the vertex nearest C240
  (x.tif, y.tif, z.tif float32, meta.json copied from the original with bbox/uuid updated).
Usage: python -I make_inputs.py [P=12] [B=48]
"""
import json
import os
import site
import sys

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
import tifffile  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import ROOT  # noqa: E402
import tif_rows  # noqa: E402

P = int(sys.argv[1]) if len(sys.argv) > 1 else 12
BAND = int(sys.argv[2]) if len(sys.argv) > 2 else 48
loc = json.load(open(os.path.join(ROOT, "out", "locate.json")))
P11 = np.array(loc["P11_xyz"])
LAY = os.path.join(ROOT, "downloads", "mesh_layout")
BD = os.path.join(ROOT, "downloads", "mesh_bands")
os.makedirs(BD, exist_ok=True)

bands = {}
for w, m in loc["meshes"].items():
    fn = os.path.join(BD, f"{w}_band.npz")
    r0 = max(0, m["coarse_row"] - BAND)
    r1 = min(m["H"], m["coarse_row"] + BAND + 1)
    if os.path.exists(fn):
        d = np.load(fn)
        g, r0 = d["g"], int(d["row0"])
    else:
        lays = {c: json.load(open(os.path.join(LAY, f"{w}_{c}.json"))) for c in "xyz"}
        g = np.stack([tif_rows.rows(lays[c], r0, r1) for c in "xyz"], -1)
        np.savez(fn, g=g, row0=r0)
    bands[w] = (g, r0)
    valid = np.all(g != -1, -1) & np.all(np.isfinite(g), -1)
    d = np.linalg.norm(g - P11, axis=-1)
    d[~valid] = np.inf
    i, j = np.unravel_index(np.argmin(d), d.shape)
    print(f"{w}: band rows {r0}..{r0 + len(g) - 1}, nearest vertex row {r0 + i} col {j}, "
          f"{d[i, j]:.1f} vox = {d[i, j] * 1.129:.1f} um, xyz {np.round(g[i, j], 1)}")
    if w == "w035":
        ci, cj = i, j

g, r0 = bands["w035"]
assert P <= ci < len(g) - P, "patch rows outside band"
patch = g[ci - P:ci + P + 1, cj - P:cj + P + 1]
pv = np.all(patch != -1, -1)
print(f"w035 patch: rows {r0 + ci - P}..{r0 + ci + P}, cols {cj - P}..{cj + P}, valid {pv.sum()}/{pv.size}")
# grid axes vs world axes
dr = np.nanmean(np.diff(patch, axis=0).reshape(-1, 3), 0)
dc = np.nanmean(np.diff(patch, axis=1).reshape(-1, 3), 0)
print("mean step along rows (xyz):", np.round(dr, 2), " along cols:", np.round(dc, 2))
lo, hi = patch[pv].min(0), patch[pv].max(0)
print("patch bbox xyz lo", np.round(lo, 1), "hi", np.round(hi, 1), "extent um", np.round((hi - lo) * 1.129, 0))

out = os.path.join(ROOT, "out", "inputs", "w035_C240_patch")
os.makedirs(out, exist_ok=True)
for k, c in enumerate("xyz"):
    tifffile.imwrite(os.path.join(out, f"{c}.tif"), np.ascontiguousarray(patch[..., k].astype(np.float32)))
meta = json.load(open(os.path.join(ROOT, "downloads", "meta", "w035-1.129-meta.json")))
meta["bbox"] = [lo.tolist(), hi.tolist()]
meta["uuid"] = "w035_C240_patch"
meta.pop("area_vx2", None)
json.dump(meta, open(os.path.join(out, "meta.json"), "wt"), indent=4)
json.dump({"P": P, "row0": int(r0 + ci - P), "col0": int(cj - P), "center_row": int(r0 + ci), "center_col": int(cj),
           "bbox_lo_xyz": lo.tolist(), "bbox_hi_xyz": hi.tolist()},
          open(os.path.join(ROOT, "out", "inputs", "w035_C240_patch.json"), "wt"), indent=1)
print("wrote", out)
