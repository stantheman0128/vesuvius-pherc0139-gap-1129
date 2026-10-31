"""Mesh rows only (already fetched, no downloads): per height, the per-vertex spacing from w035 to w034 and w036
over the whole fetched band, and for each gap edge the closest region (same geometry as the edge controls of
pick_places.py) whose spacing never drops below 21 um. Answers "how far from the edge must an edge control move to
survive the exclusion rule". Usage: python -I edge_scan.py      Writes out/places/edge_scan.json
"""
import json
import os
import site
import sys

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pick_places as PP  # noqa: E402

ROOT = PP.ROOT
P = json.load(open(os.path.join(ROOT, "out", "places", "places.json")))
out = {}
for Zh, hp in PP.HEIGHTS.items():
    G = np.load(os.path.join(ROOT, "downloads", "mesh_rows", f"w035_z{Zh}.npz"))["g"].astype(np.float64)
    v = np.all(G != -1, -1)
    a, q = PP.angle24(G)
    z24 = q[..., 2]
    S = {}
    for w in ("w034", "w036"):
        g = np.load(os.path.join(ROOT, "downloads", "mesh_rows", f"{w}_z{Zh}.npz"))["g"].astype(np.float64)
        gv = np.all(g != -1, -1)
        tree = cKDTree(PP.densify(g, gv, k=6))
        s = np.full(v.shape, np.nan)
        s[v] = tree.query(G[v])[0] * PP.UM
        S[w] = s
    smin = np.fmin(S["w034"], S["w036"])
    ctl = P["heights"][str(Zh)]["controls"]
    nadv = ctl["window"]["advance_cols"]
    W = G.shape[1]
    res = {}
    for kind, sgn in (("low_edge", +1), ("high_edge", -1)):
        c = ctl[kind]
        hits = []
        for cs in range(1, W - 1):
            ce = cs + sgn * nadv
            if ce < 1 or ce >= W - 1 or not v[:, cs].any():
                continue
            rc = int(np.argmin(np.where(v[:, cs], np.abs(z24[:, cs] - Zh), np.inf)))
            r0, r1 = rc - PP.SEED_ROWS // 2, rc + PP.SEED_ROWS // 2
            if r0 < 0 or r1 > G.shape[0]:
                continue
            c0, c1 = min(cs, ce), max(cs, ce) + 1
            reg = (slice(r0, r1), slice(c0, c1))
            if not v[reg].all():
                continue
            ra = a[reg]
            gap_lo, gap_hi = hp["gap"]
            if kind == "low_edge" and ra.max() > gap_lo - 5:
                continue
            if kind == "high_edge" and ra.min() < gap_hi + 5:
                continue
            m = float(np.nanmin(smin[reg]))
            hits.append({"seed_col": cs, "angle_min": round(float(ra.min()), 1), "angle_max": round(float(ra.max()), 1),
                         "min_spacing_um": round(m, 1),
                         "p10_w034": round(float(np.nanpercentile(S["w034"][reg], 10)), 1),
                         "p10_w036": round(float(np.nanpercentile(S["w036"][reg], 10)), 1)})
        hits.sort(key=lambda h: -h["angle_max"] if kind == "low_edge" else h["angle_min"])
        first_ok = next((h for h in hits if h["min_spacing_um"] >= 21.0), None)
        edge = hp["gap"][0] if kind == "low_edge" else hp["gap"][1]
        res[kind] = {"edge_deg": edge, "control_as_planned": hits[0] if hits else None,
                     "closest_surviving": first_ok,
                     "degrees_further_than_planned": None if not (first_ok and hits) else round(
                         (hits[0]["angle_max"] - first_ok["angle_max"]) if kind == "low_edge"
                         else (first_ok["angle_min"] - hits[0]["angle_min"]), 1)}
    out[Zh] = res
    print(Zh, json.dumps(res))
json.dump(out, open(os.path.join(ROOT, "out", "places", "edge_scan.json"), "w"), indent=1)
