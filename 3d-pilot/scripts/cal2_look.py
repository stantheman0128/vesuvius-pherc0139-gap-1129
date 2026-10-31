"""cal2 look (draft rules C3/C4, option "hold the freeze"): CT intensity profiles along the local outward normal at
every kept dense point of the real-CT cal2 trace, steps 1 to 8, with the crossings of w034 / w035 / w036 on the same
lines. Asks: is the trace on a CT sheet, is w035 on one, and are the w034 and w035 sheets separate maxima.

Only cal2's own 30 level-3 chunks are read (out/calib/plans.json "chunks"; the envelope guard and the control check
of calib_prep.plan are re-asserted). Coordinates as in calib_compare.py: trace.npz x 8 = level-0 voxels of the
1.129 um frame, mesh rows in level-0, level-3 voxel index = level-0 / 8.
Line: o + s n, s in [-S_HALF, +S_HALF] um; n = w035 outward normal (calib_prep.outward_normals), averaged over the 4
nearest w035 vertices (inverse-distance weights), so + is outward, toward w036. CT sampled trilinearly every DS um.
Ridges: scipy find_peaks on the raw profile (no smoothing); every local maximum is kept with its prominence. A
maximum counts as a ridge (a sheet candidate) if its prominence is at least the cut; see CUTS below for why.
"On a sheet": the point's intensity is at least half-way from the line's gap level (10th percentile of the line) to
the crest of its nearest ridge, i.e. inside that ridge's half-maximum width. No um tolerance is needed.
Crossings: line against the mesh triangles (Moller-Trumbore, as ray_check.py), every crossing with |s| <= S_HALF.
Row 0 (the seed, on w035 by construction) is measured too, as step 0, for reference.
Usage (WSL venv): python -I cal2_look.py      Writes out/cal2_look/profiles.json and profiles_raw.npz
"""
import json
import os
import site
import sys

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402
from scipy.signal import find_peaks  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pick_places as PP  # noqa: E402
from calib_prep import outward_normals, umb11_xy, CH, DST  # noqa: E402
from meshgeom import densify  # noqa: E402

ROOT = PP.ROOT
NAME = "cal2_z14856_a245"
CAL = os.path.join(ROOT, "out", "calib")
OUT = os.path.join(ROOT, "out", "cal2_look")
UM = 1.129                 # level-0 voxel, um
L3UM = 8 * UM              # level-3 voxel, 9.032 um
S_HALF = 200.0             # um each side. The brief says about 120. At 150 um, w036 fell outside the line at 12 to 18
#                            of 30 to 37 points per step at steps 5 to 8, and a ridge crest beyond about 140 um cannot
#                            be found, so the w035 crossings at +100 to +113 um (columns 14 to 19, step 8) had none.
DS = 1.0                   # um between samples (about 9 per level-3 voxel; trilinear, no smoothing)
HALF = 0.5                 # "on a sheet" level, as a fraction from gap level to ridge crest
MESHES = ("w034", "w035", "w036")
CUT_MAIN = "noise3"
CUTS_WHY = {
    "noise3": "3 x sqrt(2) x sigma: a ridge's prominence is a difference of two noisy samples (std sigma x sqrt 2), "
              "sigma = the larger of the two in-plane voxel noise estimates. A noise-only maximum reaches this in "
              "about 0.1% of cases. Main cut: it asks only that a dip is not noise.",
    "noise2": "2 x sqrt(2) x sigma, a lenient version of the same rule.",
    "null_p90": "90th percentile of the prominences of all local maxima on null lines (arc tangent, same length "
                "and sampling, through the trace point and through the w035 crossing).",
    "null_p95": "95th percentile of the same null. Kept as an upper, conservative cut: the cross-sections show "
                "the tangent lines do cross thin sheets (sheets here are wavy and bundled), so this null contains "
                "real sheet edges and overstates noise."}


def load_plan():
    p = [x for x in json.load(open(os.path.join(CAL, "plans.json"))) if x["name"] == NAME][0]
    keys = [tuple(int(v) for v in c) for c in p["chunks"]]
    assert p["control_vertices_in_chunks"] == 0
    assert PP.chunk_guard(keys) >= 30, "a cal2 chunk is within 30 deg of the gap envelope"
    return p, keys


def load_ct(keys):
    """The cal2 box from downloads/ct/vol.zarr/3 (the same bytes the real-CT trace read). Returns the zyx box as
    float32 and its origin in level-3 voxels (zyx)."""
    k = np.array(keys)
    lo = k.min(0)
    shape = tuple((k.max(0) - lo + 1) * CH)
    assert len(keys) == int(np.prod(k.max(0) - lo + 1)), "cal2 chunks are not a full box"
    arr = np.zeros(shape, np.uint8)
    for key in keys:
        fn = os.path.join(DST, "3", *map(str, key))
        assert os.path.getsize(fn) == CH ** 3, fn
        o = (np.array(key) - lo) * CH
        arr[o[0]:o[0] + CH, o[1]:o[1] + CH, o[2]:o[2] + CH] = np.fromfile(fn, np.uint8).reshape(CH, CH, CH)
    return arr.astype(np.float32), lo * CH


class CT:
    def __init__(self, arr, lo_zyx):
        self.a, self.lo = arr, np.asarray(lo_zyx, float)
        self.hi = np.array(arr.shape, float) - 1

    def _c(self, p0):
        p = np.asarray(p0, float).reshape(-1, 3)
        return np.stack([p[:, 2] / 8 - self.lo[0], p[:, 1] / 8 - self.lo[1], p[:, 0] / 8 - self.lo[2]])

    def sample(self, p0):
        """Trilinear CT at level-0 xyz points (..., 3); raises if any sample leaves the box."""
        c = self._c(p0)
        if np.any(c < 0) or np.any(c > self.hi[:, None]):
            raise ValueError("sample outside the cal2 box")
        return map_coordinates(self.a, c, order=1).reshape(np.shape(p0)[:-1])

    def sample_masked(self, p0):
        """As sample, but NaN where a sample would leave the box (figures only)."""
        c = self._c(p0)
        bad = np.any(c < 0, 0) | np.any(c > self.hi[:, None], 0)
        v = map_coordinates(self.a, np.clip(c, 0, self.hi[:, None]), order=1)
        v[bad] = np.nan
        return v.reshape(np.shape(p0)[:-1])


class Mesh:
    """Triangles of a mesh-row sub grid (level-0), split as in ray_check.hits."""

    def __init__(self, sub):
        v = np.all(sub != -1, -1)
        q = v[:-1, :-1] & v[1:, :-1] & v[:-1, 1:] & v[1:, 1:]
        tri = [np.stack([a[q], b[q], c[q]], 1) for a, b, c in
               ((sub[:-1, :-1], sub[1:, :-1], sub[1:, 1:]), (sub[:-1, :-1], sub[1:, 1:], sub[:-1, 1:]))]
        self.tri = np.concatenate(tri)
        cen = self.tri.mean(1)
        self.tree = cKDTree(cen)
        self.rmax = float(np.linalg.norm(self.tri - cen[:, None], axis=-1).max())
        self.dense = cKDTree(densify(sub, v))
        self.sub, self.v = sub, v

    def cross(self, o, d, tmax):
        """Line parameters t (level-0 voxels, d unit) of every crossing with |t| <= tmax, sorted by |t|."""
        idx = self.tree.query_ball_point(o, tmax + self.rmax)
        if not idx:
            return []
        a, b, c = (self.tri[idx, k] for k in range(3))
        e1, e2 = b - a, c - a
        h = np.cross(d, e2)
        det = np.einsum("ij,ij->i", e1, h)
        ok = np.abs(det) > 1e-9
        f = np.where(ok, 1.0 / np.where(ok, det, 1), 0)
        s = o - a
        u = f * np.einsum("ij,ij->i", s, h)
        qq = np.cross(s, e1)
        w = f * (qq @ d)
        t = f * np.einsum("ij,ij->i", e2, qq)
        m = ok & (u >= -1e-9) & (w >= -1e-9) & (u + w <= 1 + 1e-9) & (np.abs(t) <= tmax)
        out = []
        for x in sorted(np.round(t[m], 3).tolist(), key=abs):   # a line through a shared edge hits both triangles
            if all(abs(x - y) > 0.05 for y in out):
                out.append(x)
        return out


def load_meshes(p):
    """Sub grids within 15 deg of the place (calib_compare.py), as Mesh; plus a w035 normal function."""
    lo_a, hi_a = p["angle_range"][0] - 15, p["angle_range"][1] + 15
    res, normal_at = {}, None
    for w in MESHES:
        g = np.load(os.path.join(ROOT, "downloads", "mesh_rows", f"{w}_z{p['height']}.npz"))["g"].astype(np.float64)
        v = np.all(g != -1, -1)
        a, _ = PP.angle24(g)
        sel = v & (a >= lo_a) & (a <= hi_a)
        rr, cc = np.nonzero(sel.any(1))[0], np.nonzero(sel.any(0))[0]
        sub = g[rr.min():rr.max() + 1, cc.min():cc.max() + 1]
        res[w] = Mesh(sub)
        if w == "w035":
            sv = np.all(sub != -1, -1)
            nn = outward_normals(sub, umb11_xy)[sv]
            tree = cKDTree(sub[sv])

            def normal_at(q, tree=tree, nn=nn):
                d, j = tree.query(q, k=4)
                n = (nn[j] * (1.0 / (d + 1.0))[..., None]).sum(-2)
                return n / np.linalg.norm(n, axis=-1, keepdims=True)
    return res, normal_at


def load_trace(v="real"):
    z = np.load(os.path.join(CAL, NAME, f"trace_{v}", "trace.npz"))
    P = np.stack([z[c] for c in "xyz"], -1).astype(np.float64) * 8
    return P, np.all(P >= 0, -1), z["conf"]


def ribbon(ct, normal_at, pts, du, ext, s_um, masked=False):
    """CT on a ribbon: points every du (level-0) along the polyline pts (level-0, extended by ext at both ends along
    the end segments), each with the CT along its own w035 normal at offsets s_um. Returns u (level-0, 0 at pts[0]),
    the polyline arc position of every vertex, the base points, the normals and the (len(u), len(s)) image."""
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    cum = np.r_[0.0, np.cumsum(seg)]
    u = np.arange(-ext, cum[-1] + ext + du / 2, du)
    base = np.stack([np.interp(u, cum, pts[:, k]) for k in range(3)], -1)
    d0, d1 = (pts[1] - pts[0]) / seg[0], (pts[-1] - pts[-2]) / seg[-1]
    lo, hi = u < 0, u > cum[-1]
    base[lo] = pts[0] + u[lo, None] * d0
    base[hi] = pts[-1] + (u[hi] - cum[-1])[:, None] * d1
    n = normal_at(base)
    q = base[:, None, :] + (np.asarray(s_um) / UM)[None, :, None] * n[:, None, :]
    img = ct.sample_masked(q) if masked else ct.sample(q)
    return u, cum, base, n, img


def ridges(prof, cut):
    idx, pr = find_peaks(prof, prominence=(None, None))
    keep = pr["prominences"] >= cut
    return idx[keep], pr["prominences"][keep]


def describe(prof, s, cut, cr):
    """Ridges on one line and where the trace (s = 0) and each mesh crossing sit relative to them."""
    ridx, rprom = ridges(prof, cut)
    floor = float(np.percentile(prof, 10))
    rs = s[ridx]
    r = {"gap_level": round(floor, 1), "ridges_s_um": [round(float(x), 1) for x in rs],
         "ridge_heights": [round(float(prof[k]), 1) for k in ridx],
         "ridge_prominences": [round(float(x), 1) for x in rprom]}

    def at(x):
        if x is None or ridx.size == 0:
            return None
        k = int(np.argmin(np.abs(rs - x)))
        i = int(np.argmin(np.abs(s - x)))
        lev = (prof[i] - floor) / max(prof[ridx[k]] - floor, 1e-6)
        return {"s_um": round(float(x), 1), "ridge": k, "dist_to_ridge_um": round(float(x - rs[k]), 1),
                "intensity": round(float(prof[i]), 1), "level": round(float(lev), 2), "on_sheet": bool(lev >= HALF)}

    pos = {"trace": 0.0, **{w: (cr[w][0] if cr[w] else None) for w in MESHES}}
    for k, x in pos.items():
        r[k] = at(x)
    # width of the trace's ridge at half level (contiguous run above the half level around the crest)
    if r["trace"] is not None:
        k0 = ridx[r["trace"]["ridge"]]
        half = floor + HALF * (prof[k0] - floor)
        a = k0
        while a > 0 and prof[a - 1] >= half:
            a -= 1
        b = k0
        while b < len(prof) - 1 and prof[b + 1] >= half:
            b += 1
        r["trace_ridge_width_um"] = round(float(s[b] - s[a]), 1)
        r["trace_ridge_width_clipped"] = bool(a == 0 or b == len(prof) - 1)

    def pair(a, b):
        """Can the CT tell the two meshes' sheets apart on this line? Judged only where both crossings sit on a
        sheet: "separate" = two ridges with a dip of at least cut between them, "merged" = one ridge holds both."""
        A, B = r[a], r[b]
        if A is None or B is None:
            return None, None
        if not (A["on_sheet"] and B["on_sheet"]):
            return "not judged (a mesh is off any sheet)", None
        if A["ridge"] == B["ridge"]:
            return "merged", None
        ka, kb = sorted((A["ridge"], B["ridge"]))
        ia, ib = ridx[ka], ridx[kb]
        dip = float(min(prof[ia], prof[ib]) - prof[ia:ib + 1].min())
        return ("separate" if dip >= cut else "merged"), round(dip, 1)

    for a, b in (("w034", "w035"), ("w035", "w036")):
        r[f"{a}_{b}"], r[f"{a}_{b}_dip"] = pair(a, b)
    # n ridges strictly between the w034 and w035 crossings (a third sheet the meshes do not have would show here)
    if pos["w034"] is not None and pos["w035"] is not None:
        lo_, hi_ = sorted((pos["w034"], pos["w035"]))
        r["ridges_between_w034_w035"] = int(np.sum((rs > lo_ + L3UM) & (rs < hi_ - L3UM)))
    # what the trace is doing relative to w035 (and w034)
    T, M5, M4 = r["trace"], r["w035"], r["w034"]
    m4_on = M4 is not None and M4["on_sheet"]
    if T is None or M5 is None:
        cat = "no ridge or no w035 crossing"
    elif not T["on_sheet"]:
        cat = "trace off any sheet"
    elif T["ridge"] == M5["ridge"] and M5["on_sheet"]:
        cat = ("trace and w035 on one ridge, w034 on it too (merged)" if (m4_on and M4["ridge"] == T["ridge"])
               else "trace and w035 on one ridge")
    elif T["ridge"] == M5["ridge"]:
        cat = "trace on a ridge, w035 off it and off any other sheet"
    elif m4_on and T["ridge"] == M4["ridge"]:
        cat = ("trace on w034's ridge, w035 on another ridge" if M5["on_sheet"] else
               "trace on w034's ridge, w035 off any sheet")
    else:
        cat = "trace on a ridge that neither w034 nor w035 sits on"
    r["category"] = cat
    return r


def main():
    os.makedirs(OUT, exist_ok=True)
    p, keys = load_plan()
    arr, lo = load_ct(keys)
    ct = CT(arr, lo)
    meshes, normal_at = load_meshes(p)
    P, kept, conf = load_trace("real")
    s = np.arange(-S_HALF, S_HALF + DS / 2, DS)
    tmax = S_HALF / UM

    # voxel noise: neighbour differences along x and z (in the sheets' plane here; the normal is 7.5 deg off -y),
    # over the y rows the profiles use (box y 100..200 = level-3 y 1508..1608)
    sig = {}
    for ax, nm in ((2, "x"), (0, "z")):
        dd = np.diff(arr[:, 100:200, :], axis=ax).ravel()
        sig[nm] = float(1.4826 * np.median(np.abs(dd - np.median(dd))) / np.sqrt(2))
    sigma = max(sig.values())

    pts, raw, null_prom = [], [], []
    adv = P[-1] - P[0]
    trace_n = outward_normals(np.where(kept[..., None], P, np.nan), umb11_xy)
    for i in range(0, P.shape[0]):           # row 0 = the seed (on w035 by construction), kept as a reference
        for c in range(P.shape[1]):
            if not kept[i, c]:
                continue
            o = P[i, c]
            n = normal_at(o)
            prof = ct.sample(o[None] + (s / UM)[:, None] * n[None])
            cr = {w: [round(t * UM, 2) for t in meshes[w].cross(o, n, tmax)] for w in MESHES}
            near = {w: round(float(meshes[w].dense.query(o)[0] * UM), 2) for w in MESHES}
            t = adv[c] - (adv[c] @ n) * n
            t /= np.linalg.norm(t)
            for off in [0.0] + ([cr["w035"][0]] if cr["w035"] else []):
                nl = ct.sample(o[None] + (off / UM) * n[None] + (s / UM)[:, None] * t[None])
                null_prom += find_peaks(nl, prominence=(None, None))[1]["prominences"].tolist()
            tn = trace_n[i, c]
            ang = (None if not np.all(np.isfinite(tn)) else
                   round(float(np.degrees(np.arccos(min(1.0, abs(tn @ n) / np.linalg.norm(tn))))), 1))
            idx, pr = find_peaks(prof, prominence=(None, None))
            pts.append({"step": i, "col": c, "control": c % 6 == 0, "conf": round(float(conf[i, c]), 3),
                        "xyz_l0": np.round(o, 2).tolist(), "normal_xyz": np.round(n, 4).tolist(),
                        "angle_trace_normal_to_line_deg": ang, "crossings_s_um": cr, "nearest_surface_um": near,
                        "maxima_all": [[round(float(s[k]), 1), round(float(prof[k]), 1), round(float(q), 1)]
                                       for k, q in zip(idx, pr["prominences"])]})
            raw.append(prof.astype(np.float32))
    null_prom = np.asarray(null_prom)
    cuts = {"noise3": 3 * np.sqrt(2) * sigma, "noise2": 2 * np.sqrt(2) * sigma,
            "null_p90": float(np.percentile(null_prom, 90)), "null_p95": float(np.percentile(null_prom, 95))}
    per_cut = {cn: [describe(pr, s, cut, pt["crossings_s_um"]) for pt, pr in zip(pts, raw)]
               for cn, cut in cuts.items()}
    for k, pt in enumerate(pts):
        pt.update(per_cut[CUT_MAIN][k])
        pt["category_other_cuts"] = {cn: per_cut[cn][k]["category"] for cn in cuts if cn != CUT_MAIN}
    out = {"place": NAME, "height": p["height"], "angle_range": p["angle_range"], "direction": p["direction"],
           "chunks_read": [list(k) for k in keys], "ct_box_lo_l3_zyx": lo.tolist(),
           "conventions": {
               "s_um": "along the w035 outward normal (4 nearest vertices); 0 = trace point; + outward (toward w036)",
               "line_half_length_um": S_HALF, "sample_step_um": DS, "interp": "trilinear, no smoothing",
               "level3_voxel_um": L3UM, "on_sheet": f"level >= {HALF}: (I - gap level) / (crest - gap level), gap "
                                                     "level = 10th percentile of the line, crest = nearest ridge",
               "maxima_all": "[s_um, intensity, prominence] of every local maximum, no cut",
               "coords": "trace x 8 and mesh rows in level-0 (1.129 um); level-3 index = level-0 / 8"},
           "noise": {"sigma_voxel": {k: round(v, 2) for k, v in sig.items()},
                     "how": "1.4826 x MAD of neighbour differences / sqrt 2, level-3 voxels"},
           "cuts": {k: round(float(v), 2) for k, v in cuts.items()}, "cut_main": CUT_MAIN, "cuts_why": CUTS_WHY,
           "null": {"n_maxima": int(null_prom.size),
                    "prominence_percentiles": {str(q): round(float(np.percentile(null_prom, q)), 2)
                                               for q in (50, 75, 90, 95, 99)}},
           "summary_by_cut": {cn: summarise(pts, per_cut[cn]) for cn in cuts},
           "separation_by_mesh_spacing": {cn: {"steps_1_to_8": by_spacing(pts, per_cut[cn], range(1, 9)),
                                               "steps_5_to_8": by_spacing(pts, per_cut[cn], range(5, 9))}
                                          for cn in cuts},
           "ct_ridges": ct_ridges(raw, s, cuts),
           "stacked_means": stacked(pts, raw, s),
           "points": pts}
    json.dump(out, open(os.path.join(OUT, "profiles.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(OUT, "profiles_raw.npz"), s_um=s, prof=np.stack(raw),
                        step_col=np.array([[pt["step"], pt["col"]] for pt in pts]))
    report(out)


def tally(v):
    out = {}
    for x in v:
        out[str(x)] = out.get(str(x), 0) + 1
    return out


def by_spacing(pts, res, steps):
    """Over the given steps: how often the w034 / w035 (and w035 / w036) sheets are separate ridges, by the meshes'
    own spacing on that line. Only lines where both crossings sit on a sheet are judged."""
    bins = [(0, 60), (60, 80), (80, 100), (100, 1e9)]
    out = {}
    for a, b in (("w034", "w035"), ("w035", "w036")):
        rows = []
        for lo_, hi_ in bins:
            sep = mer = off = 0
            for pt, r in zip(pts, res):
                if pt["step"] not in steps or r[f"{a}_{b}"] is None:
                    continue
                d = abs(pt["crossings_s_um"][b][0] - pt["crossings_s_um"][a][0])
                if not lo_ <= d < hi_:
                    continue
                v = r[f"{a}_{b}"]
                sep += v == "separate"
                mer += v == "merged"
                off += v.startswith("not judged")
            rows.append({"mesh_spacing_um": [lo_, None if hi_ > 1e8 else hi_], "separate": sep, "merged": mer,
                         "not_judged_mesh_off_sheet": off})
        out[f"{a}_{b}"] = rows
    return out


def summarise(pts, res):
    out = []
    for i in sorted({pt["step"] for pt in pts}):
        ks = [k for k, pt in enumerate(pts) if pt["step"] == i]

        def cnt(fn, ks=ks):
            v = [fn(res[k]) for k in ks]
            v = [x for x in v if x is not None]
            return [int(sum(v)), len(v)]

        def med(fn, ks=ks):
            v = [fn(k) for k in ks]
            v = [x for x in v if x is not None]
            return round(float(np.median(v)), 1) if v else None

        cats = {}
        for k in ks:
            cats[res[k]["category"]] = cats.get(res[k]["category"], 0) + 1
        out.append({
            "step": i, "n": len(ks),
            "trace_on_sheet": cnt(lambda r: r["trace"]["on_sheet"] if r["trace"] else None),
            "w035_on_sheet": cnt(lambda r: r["w035"]["on_sheet"] if r["w035"] else None),
            "w034_on_sheet": cnt(lambda r: r["w034"]["on_sheet"] if r["w034"] else None),
            "w036_on_sheet": cnt(lambda r: r["w036"]["on_sheet"] if r["w036"] else None),
            "trace_level_median": med(lambda k: res[k]["trace"]["level"] if res[k]["trace"] else None),
            "w035_level_median": med(lambda k: res[k]["w035"]["level"] if res[k]["w035"] else None),
            "w034_level_median": med(lambda k: res[k]["w034"]["level"] if res[k]["w034"] else None),
            "trace_abs_dist_to_ridge_median_um": med(lambda k: abs(res[k]["trace"]["dist_to_ridge_um"])
                                                     if res[k]["trace"] else None),
            "w035_abs_dist_to_ridge_median_um": med(lambda k: abs(res[k]["w035"]["dist_to_ridge_um"])
                                                    if res[k]["w035"] else None),
            "trace_ridge_width_median_um": med(lambda k: res[k].get("trace_ridge_width_um")),
            "w034_w035": tally([res[k]["w034_w035"] for k in ks]),
            "w035_w036": tally([res[k]["w035_w036"] for k in ks]),
            "ridges_between_w034_w035": cnt(lambda r: None if "ridges_between_w034_w035" not in r else
                                            r["ridges_between_w034_w035"] > 0),
            "s_w034_median_um": med(lambda k: pts[k]["crossings_s_um"]["w034"][0]
                                    if pts[k]["crossings_s_um"]["w034"] else None),
            "s_w035_median_um": med(lambda k: pts[k]["crossings_s_um"]["w035"][0]
                                    if pts[k]["crossings_s_um"]["w035"] else None),
            "s_w036_median_um": med(lambda k: pts[k]["crossings_s_um"]["w036"][0]
                                    if pts[k]["crossings_s_um"]["w036"] else None),
            "nearest_surface_w035_median_um": med(lambda k: pts[k]["nearest_surface_um"]["w035"]),
            "categories": dict(sorted(cats.items(), key=lambda kv: -kv[1]))})
    return out


def ct_ridges(raw, s, cuts):
    """What the CT itself shows, without the meshes: spacing between adjacent ridges on every line, and each ridge's
    width at half level (gap level = 10th percentile of the line). Ridges whose half-level run reaches the line's end
    are left out of the widths."""
    out = {}
    for cn, cut in cuts.items():
        pitch, width = [], []
        for prof in raw:
            ridx, _ = ridges(prof, cut)
            pitch += np.diff(s[ridx]).tolist()
            floor = np.percentile(prof, 10)
            for k0 in ridx:
                half = floor + HALF * (prof[k0] - floor)
                a = b = k0
                while a > 0 and prof[a - 1] >= half:
                    a -= 1
                while b < len(prof) - 1 and prof[b + 1] >= half:
                    b += 1
                if a > 0 and b < len(prof) - 1:
                    width.append(s[b] - s[a])
        pitch, width = np.array(pitch), np.array(width)
        out[cn] = {"adjacent_ridge_spacing_um": {"n": int(pitch.size),
                                                 **{f"p{q}": round(float(np.percentile(pitch, q)), 1)
                                                    for q in (5, 10, 25, 50, 75)},
                                                 "share_below_60": round(float(np.mean(pitch < 60)), 3),
                                                 "share_below_40": round(float(np.mean(pitch < 40)), 3)},
                   "ridge_half_level_width_um": {"n": int(width.size),
                                                 **{f"p{q}": round(float(np.percentile(width, q)), 1)
                                                    for q in (10, 25, 50, 75, 90)}}}
    return out


def stacked(pts, raw, s):
    """Per step, the mean profile aligned at the trace point, at the w035 crossing and at the w034 crossing
    (+-100 um around the alignment point), every 2 um, with the local maxima of each mean."""
    W = int(round(100 / DS))
    rel = np.arange(-W, W + 1) * DS
    res = []
    for i in sorted({pt["step"] for pt in pts}):
        row = {"step": i, "rel_um_every_2": rel[::2].tolist()}
        for tag in ("trace", "w035", "w034"):
            segs = []
            for pt, prof in zip(pts, raw):
                if pt["step"] != i:
                    continue
                x = 0.0 if tag == "trace" else (pt["crossings_s_um"][tag][0] if pt["crossings_s_um"][tag] else None)
                if x is None:
                    continue
                k = int(np.argmin(np.abs(s - x)))
                if W <= k < len(s) - W:
                    segs.append(prof[k - W:k + W + 1])
            if not segs:
                row[tag] = None
                continue
            m = np.mean(segs, 0)
            idx, pr = find_peaks(m, prominence=(None, None))
            row[tag] = {"n": len(segs), "mean": np.round(m[::2], 1).tolist(),
                        "maxima": [[round(float(rel[k]), 1), round(float(m[k]), 1), round(float(q), 1)]
                                   for k, q in zip(idx, pr["prominences"])]}
        res.append(row)
    return res


def report(out):
    print("noise", out["noise"]["sigma_voxel"], "cuts", out["cuts"], "null", out["null"])
    print("ct ridges", json.dumps(out["ct_ridges"]))
    for cn, summ in out["summary_by_cut"].items():
        print(f"--- cut {cn} = {out['cuts'][cn]}")
        for r in summ:
            print(f"  step {r['step']} n {r['n']}: on sheet trace {r['trace_on_sheet']} w035 {r['w035_on_sheet']} "
                  f"w034 {r['w034_on_sheet']} w036 {r['w036_on_sheet']} | level t/35/34 {r['trace_level_median']}/"
                  f"{r['w035_level_median']}/{r['w034_level_median']} | |d ridge| t {r['trace_abs_dist_to_ridge_median_um']}"
                  f" w035 {r['w035_abs_dist_to_ridge_median_um']} | width {r['trace_ridge_width_median_um']} | between "
                  f"{r['ridges_between_w034_w035']} | s34/35/36 {r['s_w034_median_um']}/{r['s_w035_median_um']}/"
                  f"{r['s_w036_median_um']}")
            print("       34/35", r["w034_w035"], "35/36", r["w035_w036"])
            print("      ", r["categories"])
        print("  separation by mesh spacing", json.dumps(out["separation_by_mesh_spacing"][cn]))


if __name__ == "__main__":
    main()
