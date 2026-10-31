"""Numbers for PILOT-RULES-3D-DRAFT.md from sheet traces already on disk (no new data).

Per step, at the control points only (dense columns 0, 6, ..., 36: the 37 columns are Hermite samples between the
7 control points, inference.py:113), the number kept and the number whose nearest published mesh is w035.
Per trace, the nearest-mesh margin m = d(w035) / min(d(w034), d(w036)) at every kept dense point (m < 1 means
assigned to w035), and the local w035-to-neighbour spacings at the w035 vertex nearest each point.
Binomial table: smallest k of n with P(X >= k | n, p0) <= alpha, for the p0/alpha options in the draft.
Usage: python -I rules_numbers.py      Writes out/rules_numbers.json
"""
import json
import math
import os
import site
import sys

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import ROOT  # noqa: E402

UM = 1.129
OUT = os.path.join(ROOT, "out")


def densify(g, v, k=8):
    q = v[:-1, :-1] & v[1:, :-1] & v[:-1, 1:] & v[1:, 1:]
    ii, jj = np.nonzero(q)
    t = (np.arange(k) + 0.5) / k
    a, b = np.meshgrid(t, t, indexing="ij")
    a, b = a.ravel()[None, :, None], b.ravel()[None, :, None]
    p00, p10, p01, p11 = (g[ii, jj][:, None], g[ii + 1, jj][:, None], g[ii, jj + 1][:, None],
                          g[ii + 1, jj + 1][:, None])
    pts = (1 - a) * (1 - b) * p00 + a * (1 - b) * p10 + (1 - a) * b * p01 + a * b * p11
    return np.concatenate([pts.reshape(-1, 3), g[v]])


trees, verts = {}, {}
for w in ("w034", "w035", "w036"):
    g = np.load(os.path.join(ROOT, "downloads", "mesh_bands", f"{w}_band.npz"))["g"].astype(np.float64)
    v = np.all(g != -1, -1) & np.all(np.isfinite(g), -1)
    trees[w] = cKDTree(densify(g, v))
    verts[w] = g[v]
vtree35 = cKDTree(verts["w035"])
res = {}
for tag, d in (("real_t03", "sheet7_s8_t03"), ("shift_out", "ablate/trace_shift_p"), ("shift_in", "ablate/trace_shift_m")):
    z = np.load(os.path.join(OUT, d, "trace.npz"))
    P = np.stack([z[c] for c in "xyz"], -1).astype(np.float64) * 8
    kept = np.all(P >= 0, -1)
    steps = []
    for r in range(1, P.shape[0]):
        cp = [c for c in range(0, P.shape[1], 6) if kept[r, c]]
        if cp:
            q = P[r, cp]
            dd = np.stack([trees[w].query(q)[0] for w in ("w034", "w035", "w036")])
            k = int((np.argmin(dd, 0) == 1).sum())
        else:
            k = 0
        steps.append({"step": r, "control_kept": len(cp), "control_nearest_w035": k})
    q = P[1:][kept[1:]]
    d = {w: trees[w].query(q)[0] * UM for w in ("w034", "w035", "w036")}
    m = d["w035"] / np.minimum(d["w034"], d["w036"])
    _, vi = vtree35.query(q)
    s36 = trees["w036"].query(verts["w035"][vi])[0] * UM
    s34 = trees["w034"].query(verts["w035"][vi])[0] * UM
    res[tag] = {"steps": steps,
                "control_points_pooled": [sum(s["control_nearest_w035"] for s in steps),
                                          sum(s["control_kept"] for s in steps)],
                "margin_m": {"median": round(float(np.median(m)), 3), "p90": round(float(np.percentile(m, 90)), 3),
                             "max": round(float(m.max()), 3), "share_ge_0.5": round(float((m >= 0.5).mean()), 3)},
                "spacing_w035_w036_um": [round(float(np.percentile(s36, p)), 1) for p in (10, 50, 90)],
                "spacing_w035_w034_um": [round(float(np.percentile(s34, p)), 1) for p in (10, 50, 90)]}


def kmin(n, p0, alpha):
    for k in range(n + 1):
        if sum(math.comb(n, j) * p0 ** j * (1 - p0) ** (n - j) for j in range(k, n + 1)) <= alpha:
            return k
    return None


res["binomial_kmin"] = {f"p0={p0},alpha={a}": {n: kmin(n, p0, a) for n in range(1, 8)}
                        for p0 in (0.5, 0.8) for a in (0.05, 0.01)}
res["binomial_kmin_pooled"] = {f"p0={p0},alpha={a}": {n: kmin(n, p0, a) for n in (28, 42, 56)}
                               for p0 in (0.5, 0.8, 0.9) for a in (0.05, 0.01)}
json.dump(res, open(os.path.join(OUT, "rules_numbers.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
