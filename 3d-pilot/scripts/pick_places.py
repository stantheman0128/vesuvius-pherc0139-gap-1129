"""Controls and calibration places for PILOT-RULES-3D-DRAFT.md, from mesh rows only (no CT is read here).

Step A (fetch): at each of the four heights (2.399 um frame z 12,400 / 14,856 / 17,000 / 19,500), the rows of w035
(every 2nd row) and w034 / w036 (every 3rd row) within about 2.2 mm of the height, and only the columns whose angle
(coarse 16-row sample, every coarse row of the window) lies on the non-gap arc [high edge + 2 deg, low edge - 2 deg].
So no gap-row column of any mesh is fetched. Cached in downloads/mesh_rows/.
Step B (places): planned traces are w035 column lines (seed of SEED_ROWS rows centred on the height) advanced along
the arc by 8 steps of 20 level-3 voxels (1.445 mm). The planned trace region is the w035 grid rows of the seed by the
columns from the seed column to the end column.
  - 12 controls (draft section 5): per height one centred in the control window and one at each gap edge, every
    vertex of the region at least 5 deg outside the edge, as close to the edge as that allows. Direction: toward
    the gap (low edge +, high edge -; window toward the nearer edge).
  - Calibration candidates: region centres every 5 deg on the non-gap arc, direction toward the nearer gap edge,
    eligible if (1) every level-3 CT chunk of its box has all 8 corners at least 30 deg outside the 327..85 deg
    envelope (per-height umbilicus, 2.399 frame), (2) its CT box contains no vertex of any control region,
    (3) no spacing below 21 um (exclusion rule), (4) w035 valid over the region.
    Pick: the two narrowest (p10 to the nearer neighbour) at two different heights, then at each remaining height
    the candidate whose p10 is closest to that height's median candidate p10.
Spacing at a w035 vertex = distance to the w034 / w036 surface (fetched rows, each grid cell upsampled 8 x 8).
Usage: python -I pick_places.py           Writes out/places/places.json
"""
import json
import math
import os
import site
import sys
import threading
from concurrent.futures import ThreadPoolExecutor

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import netlog  # noqa: E402
import tif_rows  # noqa: E402

ROOT = netlog.ROOT
UM = 1.129
HEIGHTS = {12400: {"gap": (328, 77), "win": (130, 170)}, 14856: {"gap": (334, 80), "win": (130, 170)},
           17000: {"gap": (338, 79), "win": (150, 190)}, 19500: {"gap": (346, 83), "win": (170, 210)}}
SEED_ROWS = 86           # band rows (w035 fetched every 2nd row); trimmed below to the C240 seed length (364.9 level-3 voxels gives 7 points at spacing 60)
SEED_L3_TARGET = (370.0, 410.0)
STEP_L3, STEPS, PAD_L3 = 20.0, 8, 64.0
ADV_UM = STEP_L3 * STEPS * 8 * UM            # 1445 um of arc
DZ_WIN = 900.0                               # 2.399-frame voxels around the height (2.16 mm)
OUTD = os.path.join(ROOT, "out", "places")
CACHE = os.path.join(ROOT, "downloads", "mesh_rows")
os.makedirs(OUTD, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

T = np.array(json.load(open(os.path.join(ROOT, "downloads", "meta", "zarr", "transform.json")))
             ["transformation_matrix"], float)
cp = sorted(json.load(open(os.path.join(ROOT, "downloads", "meta", "umbilicus.json")))["control_points"],
            key=lambda c: c["z"])
UZ = np.array([c["z"] for c in cp], float)
UX = np.array([c["x"] for c in cp], float)
UY = np.array([c["y"] for c in cp], float)


def to24(p):
    return p @ T[:, :3].T + T[:, 3]


def angle24(p11):
    q = to24(p11)
    ux, uy = np.interp(q[..., 2], UZ, UX), np.interp(q[..., 2], UZ, UY)
    return np.degrees(np.arctan2(q[..., 1] - uy, q[..., 0] - ux)) % 360, q


def on_arc(a, lo, hi):
    """a on the arc from lo increasing to hi (lo < hi, both on the non-gap side)."""
    return (a >= lo) & (a <= hi)


def in_envelope(a):
    return (a >= 327) | (a <= 85)


# ---- cap check without re-reading the whole ledger on every request ----
_BASE = netlog.total_fetched()
_ACC = [0]
_LK = threading.Lock()
netlog.total_fetched = lambda: _BASE + _ACC[0]
_orig_fetch = tif_rows.fetch


def _counted_fetch(*a, **k):
    data = _orig_fetch(*a, **k)
    with _LK:
        _ACC[0] += len(data)
    return data


def fetch_rows(w, Zh):
    fn = os.path.join(CACHE, f"{w}_z{Zh}.npz")
    if os.path.exists(fn):
        d = np.load(fn)
        return d["g"], d["rows"], d["colmask"]
    gap_lo, gap_hi = HEIGHTS[Zh]["gap"]
    arc = (gap_hi + 2, gap_lo - 2)
    g16 = np.load(os.path.join(ROOT, "downloads", "mesh_layout", f"{w}_coarse16.npz"))["g"].astype(np.float64)
    v16 = np.all(g16 != -1, -1) & np.all(np.isfinite(g16), -1)
    a16, p16 = angle24(g16)
    near = v16 & (np.abs(p16[..., 2] - Zh) <= DZ_WIN) & on_arc(a16, *arc)
    rr = np.nonzero(near.any(1))[0]
    r0, r1 = int(rr.min() * 16 - 16), int(rr.max() * 16 + 17)
    win = slice(max(rr.min() - 1, 0), rr.max() + 2)
    vv, aa = v16[win], a16[win]
    colmask = vv.any(0) & np.all(~vv | on_arc(aa, *arc), 0)
    lays = {c: json.load(open(os.path.join(ROOT, "downloads", "mesh_layout", f"{w}_{c}.json"))) for c in "xyz"}
    H, W = lays["x"]["height"], lays["x"]["width"]
    assert lays["x"]["rows_per_strip"] == H, "expected a single strip"
    r0, r1 = max(r0, 0), min(r1, H)
    step = 2 if w == "w035" else 3
    rows = np.arange(r0, r1, step)
    runs, c = [], 0
    while c < W:
        if colmask[c]:
            s = c
            while c < W and colmask[c]:
                c += 1
            runs.append((s, c))
        else:
            c += 1
    g = np.full((len(rows), W, 3), -1.0, np.float32)
    jobs = [(i, r, k, s, e) for i, r in enumerate(rows) for k in range(3) for (s, e) in runs]

    def one(job):
        i, r, k, s, e = job
        lay = lays["xyz"[k]]
        start = lay["offsets"][0] + (int(r) * W + s) * 4
        buf = tif_rows.fetch(lay["url"], start, start + (e - s) * 4 - 1, tag=f"tif-cols-{w}-z{Zh}")
        g[i, s:e, k] = np.frombuffer(buf, "<f4")
        return len(buf)

    tif_rows.fetch = _counted_fetch
    with ThreadPoolExecutor(16) as ex:
        nbytes = sum(ex.map(one, jobs))
    np.savez(fn, g=g, rows=rows, colmask=colmask)
    print(f"  fetched {w} z{Zh}: rows {r0}..{r1} step {step} ({len(rows)}), {int(colmask.sum())} of {W} cols in "
          f"{len(runs)} runs, {nbytes} bytes, {len(jobs)} requests", flush=True)
    return g, rows, colmask


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


def pct(v):
    return [round(float(np.percentile(v, p)), 1) for p in (10, 50, 90)]


def chunk_plan(seed11, end11):
    pts = np.concatenate([seed11, end11]) / 8.0
    lo, hi = pts.min(0) - PAD_L3, pts.max(0) + PAD_L3
    ch = [[z, y, x] for z in range(int(lo[2] // 128), int(hi[2] // 128) + 1)
          for y in range(int(lo[1] // 128), int(hi[1] // 128) + 1)
          for x in range(int(lo[0] // 128), int(hi[0] // 128) + 1)]
    return lo, hi, ch


def chunk_guard(ch):
    """Smallest angular distance (deg) of any chunk corner to the 327..85 envelope; negative if inside."""
    worst = 1e9
    for z, y, x in ch:
        corners = np.array([[(x + dx) * 1024.0, (y + dy) * 1024.0, (z + dz) * 1024.0]
                            for dz in (0, 1) for dy in (0, 1) for dx in (0, 1)])
        a, _ = angle24(corners)
        for v in a:
            dist = -1.0 if in_envelope(v) else min(v - 85, 327 - v)
            worst = min(worst, dist)
    return worst


def build_region(G, rows, Zh, a_center=None, edge=None, direction=None):
    """Region on w035 for a seed column; returns dict or None. a_center: centre angle of the region; edge: (kind,
    angle) to place the region as close to that edge as the 5 deg rule allows."""
    v = np.all(G != -1, -1)
    a, q = angle24(G.astype(np.float64))
    z24 = q[..., 2]
    colstep_um = np.nanmedian(np.where(v[:, 1:] & v[:, :-1],
                                       np.linalg.norm(np.diff(G, axis=1), axis=-1), np.nan)) * UM
    nadv = int(round(ADV_UM / colstep_um))
    best = None
    W = G.shape[1]
    cands = range(1, W - 1)
    out = []
    for cs in cands:
        ce = cs + direction * nadv
        if ce < 1 or ce >= W - 1:
            continue
        # seed row: the row of column cs whose z24 is closest to the height
        col_ok = v[:, cs]
        if not col_ok.any():
            continue
        rc = int(np.argmin(np.where(col_ok, np.abs(z24[:, cs] - Zh), np.inf)))
        r0, r1 = rc - SEED_ROWS // 2, rc + SEED_ROWS // 2
        if r0 < 0 or r1 > G.shape[0]:
            continue
        c0, c1 = min(cs, ce), max(cs, ce) + 1
        reg = (slice(r0, r1), slice(c0, c1))
        if not v[reg].all() or not v[r0:r1, cs - direction].all():
            continue
        ra = a[reg]
        if edge is not None:
            kind, ea = edge
            if kind == "low" and ra.max() > ea - 5:
                continue
            if kind == "high" and ra.min() < ea + 5:
                continue
        out.append((cs, ce, rc, r0, r1, c0, c1, float(np.median(ra))))
    if not out:
        return None
    if edge is not None:
        kind, ea = edge
        pick = max(out, key=lambda o: o[7]) if kind == "low" else min(out, key=lambda o: o[7])
    else:
        pick = min(out, key=lambda o: abs(o[7] - a_center))
    cs, ce, rc, r0, r1, c0, c1, amed = pick
    return finish_region(G, rows, a, cs, ce, r0, r1, c0, c1, direction, colstep_um, nadv)


def finish_region(G, rows, a, cs, ce, r0, r1, c0, c1, direction, colstep_um, nadv):
    seed11 = G[r0:r1, cs].astype(np.float64)
    # trim the seed to the target length in level-3 voxels (7 front points at spacing 60)
    seg = np.linalg.norm(np.diff(seed11 / 8.0, axis=0), axis=1)
    while seg.sum() > SEED_L3_TARGET[1] and len(seed11) > 10:
        seed11, seg = seed11[1:-1], seg[1:-1]
        r0, r1 = r0 + 1, r1 - 1
    end11 = G[r0:r1, ce].astype(np.float64)
    lo, hi, ch = chunk_plan(seed11, end11)
    ra = a[r0:r1, c0:c1]
    return {"seed_col": int(cs), "end_col": int(ce), "direction": "+" if direction > 0 else "-",
            "band_rows": [int(r0), int(r1)], "global_rows": [int(rows[r0]), int(rows[r1 - 1]) + 1],
            "cols": [int(c0), int(c1)], "col_step_um": round(float(colstep_um), 2), "advance_cols": nadv,
            "seed_len_l3": round(float(seg.sum()), 1), "seed_len_ok_7pts": bool(360.0 <= seg.sum() < 420.0),
            "angle_min": round(float(ra.min()), 2), "angle_max": round(float(ra.max()), 2),
            "angle_seed": round(float(np.median(a[r0:r1, cs])), 2),
            "z24_range": [round(float(v), 1) for v in (to24(G[r0:r1, c0:c1].astype(np.float64))[..., 2].min(),
                                                       to24(G[r0:r1, c0:c1].astype(np.float64))[..., 2].max())],
            "read_lo_xyz_l3": lo.round(1).tolist(), "read_hi_xyz_l3": hi.round(1).tolist(),
            "chunks": ch, "n_chunks": len(ch)}


def spacing(G35, reg, nb):
    """p10/median/p90/min of the distance (um) from the region's w035 vertices to each neighbour surface."""
    r0, r1 = reg["band_rows"]
    c0, c1 = reg["cols"]
    P = G35[r0:r1, c0:c1].reshape(-1, 3).astype(np.float64)
    a35, _ = angle24(P)
    res = {}
    for w, (Gn, an) in nb.items():
        vn = np.all(Gn != -1, -1)
        # neighbour vertices within 10 deg of the region (enough for a nearest-point search at ~0.1 mm)
        lo, hi = reg["angle_min"] - 10, reg["angle_max"] + 10
        sel = vn & (an >= lo) & (an <= hi)
        rr = np.nonzero(sel.any(1))[0]
        cc = np.nonzero(sel.any(0))[0]
        if len(rr) < 2 or len(cc) < 2:
            res[w] = {"present": False}
            continue
        sub = Gn[rr.min():rr.max() + 1, cc.min():cc.max() + 1].astype(np.float64)
        sv = np.all(sub != -1, -1)
        tree = cKDTree(densify(sub, sv))
        d, idx = tree.query(P)
        d_um = d * UM
        # neighbour actually beside the vertex: nearest point within 2 deg of the vertex's angle
        qa, _ = angle24(tree.data[idx])
        beside = np.abs((qa - a35 + 180) % 360 - 180) <= 2.0
        res[w] = {"present": True, "p10_med_p90_um": pct(d_um), "min_um": round(float(d_um.min()), 1),
                  "share_beside": round(float(beside.mean()), 3), "n_vertices": int(len(P))}
    return res


def box_contains(reg, pts11):
    lo, hi = np.array(reg["read_lo_xyz_l3"]), np.array(reg["read_hi_xyz_l3"])
    ch = np.array(reg["chunks"])
    clo = ch.min(0)[::-1] * 128.0
    chi = (ch.max(0)[::-1] + 1) * 128.0
    p = pts11 / 8.0
    inside = np.all((p >= clo) & (p < chi), -1)
    return int(inside.sum())


def main():
    out = {"rule": __doc__, "heights": {}}
    for Zh, hp in HEIGHTS.items():
        print("height", Zh, flush=True)
        data = {w: fetch_rows(w, Zh) for w in ("w034", "w035", "w036")}
        G35, rows35, _ = data["w035"]
        nb = {}
        for w in ("w034", "w036"):
            Gn = data[w][0]
            an, _ = angle24(Gn.astype(np.float64))
            nb[w] = (Gn, an)
        gap_lo, gap_hi = hp["gap"]
        ctl = {}
        ctl["window"] = build_region(G35, rows35, Zh, a_center=sum(hp["win"]) / 2, direction=-1)
        ctl["low_edge"] = build_region(G35, rows35, Zh, edge=("low", gap_lo), direction=+1)
        ctl["high_edge"] = build_region(G35, rows35, Zh, edge=("high", gap_hi), direction=-1)
        for k, reg in ctl.items():
            if reg is None:
                continue
            reg["spacing"] = spacing(G35, reg, nb)
            if k == "window":
                reg["inside_window"] = bool(reg["angle_min"] >= hp["win"][0] and reg["angle_max"] <= hp["win"][1])
        # control vertices (for the calibration box check)
        cpts = np.concatenate([G35[r["band_rows"][0]:r["band_rows"][1], r["cols"][0]:r["cols"][1]].reshape(-1, 3)
                               for r in ctl.values() if r is not None]).astype(np.float64)
        # calibration candidates: centres every 5 deg on the non-gap arc
        cands = []
        for ac in range(gap_hi + 10, gap_lo - 9, 5):
            direction = -1 if (ac - gap_hi) <= (gap_lo - ac) else +1
            reg = build_region(G35, rows35, Zh, a_center=ac, direction=direction)
            if reg is None or abs(0.5 * (reg["angle_min"] + reg["angle_max"]) - ac) > 2.5:
                continue
            reg["centre_grid_deg"] = ac
            reg["guard_min_deg_from_envelope"] = round(chunk_guard(reg["chunks"]), 1)
            reg["control_vertices_in_box"] = box_contains(reg, cpts)
            if reg["guard_min_deg_from_envelope"] >= 30 and reg["control_vertices_in_box"] == 0:
                reg["spacing"] = spacing(G35, reg, nb)
            cands.append(reg)
        out["heights"][Zh] = {"gap_rows": hp["gap"], "control_window": hp["win"], "controls": ctl,
                              "calib_candidates": cands}
        json.dump(out, open(os.path.join(OUTD, "places.json"), "w"), indent=1)
    # calibration selection
    elig = []
    for Zh, h in out["heights"].items():
        for i, c in enumerate(h["calib_candidates"]):
            sp = c.get("spacing")
            if not sp or not all(sp[w].get("present") for w in ("w034", "w036")) or not c["seed_len_ok_7pts"]:
                continue
            if min(sp[w]["min_um"] for w in ("w034", "w036")) < 21.0:
                continue
            c["p10_nearer_um"] = min(sp[w]["p10_med_p90_um"][0] for w in ("w034", "w036"))
            elig.append((c["p10_nearer_um"], Zh, i))
    elig.sort()
    pick = []
    if elig:
        pick.append(("narrow", elig[0]))
        rest = [e for e in elig if e[1] != elig[0][1]]
        if rest:
            pick.append(("narrow", rest[0]))
        used = {p[1][1] for p in pick}
        for Zh in out["heights"]:
            if Zh in used:
                continue
            here = [e for e in elig if e[1] == Zh]
            if not here:
                continue
            med = float(np.median([e[0] for e in here]))
            pick.append(("typical", min(here, key=lambda e: abs(e[0] - med))))
    out["calibration"] = [{"kind": k, "height": e[1], "candidate_index": e[2], "p10_nearer_um": e[0],
                           "n_chunks": out["heights"][e[1]]["calib_candidates"][e[2]]["n_chunks"]} for k, e in pick]
    out["eligible_counts"] = {Zh: sum(1 for e in elig if e[1] == Zh) for Zh in out["heights"]}
    out["eligible_le_80um"] = sum(1 for e in elig if e[0] <= 80.0)
    json.dump(out, open(os.path.join(OUTD, "places.json"), "w"), indent=1)
    print("calibration picks:", out["calibration"])
    print("eligible per height:", out["eligible_counts"], "eligible with p10 <= 80 um:", out["eligible_le_80um"])
    print("CT bytes for picks:", sum(c["n_chunks"] for c in out["calibration"]) * 2097152)
    print("fetched this run:", _ACC[0], "ledger total:", netlog.total_fetched())


if __name__ == "__main__":
    main()
