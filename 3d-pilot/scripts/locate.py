"""Locate the C240 control (gap-1129 README: z 17,995 in the 2.399 um frame, 240 deg about the umbilicus,
w035 crossing at r 8382.6 um) in the 1.129 um frame, and find it on the on-1.129 tifxyz copies.

Outputs out/locate.json: P24, P11, nearest grid vertex per mesh (from a coarse row sample), and the tif layouts.
Never touches the gap sector (327 to 85 deg); asserts the place's angle is outside it.
"""
import json
import os
import site
import sys

sys.path.append(site.getusersitepackages())  # python -I drops the user site, where numpy lives
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import ROOT  # noqa: E402
import tif_rows  # noqa: E402

BUCKET = "https://vesuvius-challenge-open-data.s3.amazonaws.com"
SEG = {
    "w034": "20260303000000-w034_2026030317",
    "w035": "20260317000000-w035_2026031718",
    "w036": "20260316000000-w036_2026031607",
}
META = os.path.join(ROOT, "downloads", "meta")
LAY = os.path.join(ROOT, "downloads", "mesh_layout")
os.makedirs(LAY, exist_ok=True)
UM24 = 2.399

Z24, ANG, R_W035_UM = 17995.0, 240.0, 8382.6
assert not (ANG >= 327 or ANG <= 85), "place inside gap sector"

T = np.array(json.load(open(os.path.join(META, "zarr", "transform.json")))["transformation_matrix"], float)
A = np.eye(4)
A[:3] = T          # 1.129 (moving) -> 2.399 (fixed), xyz
Ainv = np.linalg.inv(A)

cp = sorted(json.load(open(os.path.join(META, "umbilicus.json")))["control_points"], key=lambda c: c["z"])
uz = np.array([c["z"] for c in cp], float)
ux = np.interp(Z24, uz, [c["x"] for c in cp])
uy = np.interp(Z24, uz, [c["y"] for c in cp])
r_vox = R_W035_UM / UM24
P24 = np.array([ux + r_vox * np.cos(np.radians(ANG)), uy + r_vox * np.sin(np.radians(ANG)), Z24])
P11 = (Ainv @ np.r_[P24, 1.0])[:3]
U11 = (Ainv @ np.r_[ux, uy, Z24, 1.0])[:3]
print("umbilicus 2.399:", round(ux, 1), round(uy, 1), " P24:", np.round(P24, 1), " P11 (xyz):", np.round(P11, 1))

res = {"P24_xyz": P24.tolist(), "P11_xyz": P11.tolist(), "umb24_xy": [ux, uy], "umb11_xyz": U11.tolist(),
       "meshes": {}}
STRIDE = 16
for w, sid in SEG.items():
    ts = sid.split("-")[0]
    base = f"{BUCKET}/PHerc0139/segments/{sid}/mesh/{ts}-on-20260413113053-1.129um.tifxyz"
    lays = {c: tif_rows.layout(f"{base}/{c}.tif", os.path.join(LAY, f"{w}_{c}.json")) for c in "xyz"}
    H, W = lays["x"]["height"], lays["x"]["width"]
    cache = os.path.join(LAY, f"{w}_coarse{STRIDE}.npz")
    if os.path.exists(cache):
        g = np.load(cache)["g"]
    else:
        idx = list(range(0, H, STRIDE))
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(8) as ex:
            g = np.stack([np.concatenate(list(ex.map(lambda r, c=c: tif_rows.rows(lays[c], r, r + 1), idx)))
                          for c in "xyz"], -1)
        np.savez(cache, g=g)
    valid = np.all(g != -1, -1) & np.all(np.isfinite(g), -1)
    d = np.linalg.norm(g - P11[None, None], axis=-1)
    d[~valid] = np.inf
    i, j = np.unravel_index(np.argmin(d), d.shape)
    res["meshes"][w] = {"H": H, "W": W, "coarse_row": int(i * STRIDE), "col": int(j), "dist_vox": float(d[i, j]),
                        "xyz": g[i, j].tolist(), "base": base}
    print(w, "H x W", H, W, "nearest coarse vertex row", i * STRIDE, "col", j, "dist", round(float(d[i, j]), 1),
          "vox =", round(float(d[i, j]) * 1.129, 1), "um", np.round(g[i, j], 1))
json.dump(res, open(os.path.join(ROOT, "out", "locate.json"), "wt"), indent=1)
