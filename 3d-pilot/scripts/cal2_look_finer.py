"""Apply out/cal2_look/FINER-RULE.md (fixed before any level-1/2 voxel was read) to the cal2 seed.

Levels 1 and 2: the 14 chunks of finer_ct_estimate.json (downloads/ct/vol.zarr/<L>/z/y/x). Level 3: cal2's box.
Writes out/cal2_look/finer_look.json and out/cal2_look/png/cal2_col15_L1.png; console to finer_console.txt.
Usage (WSL venv): python -I cal2_look_finer.py
"""
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cal2_look as L  # noqa: E402
import numpy as np  # noqa: E402
from fetch_ct import DST  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402
from scipy.signal import find_peaks  # noqa: E402

COLS = list(range(11, 18))
ROWS = (0, 1, 2)
S = np.arange(-130.0, 60.0 + 1e-9, 1.0)
W_PAD, UP_TOL, SEP_MIN, LO_TOL = 15.0, 15.0, 12.0, 15.0
STEP = 4.516
REACH = STEP + 2.0
MATCH = 3.0
RULE = os.path.join(L.OUT, "FINER-RULE.md")


class Box:
    def __init__(self, lv, keys):
        self.lv = lv
        k = np.array(keys)
        self.lo = k.min(0)
        shape = tuple((k.max(0) - self.lo + 1) * L.CH)
        self.a = np.full(shape, np.nan, np.float32)
        diffs = {"x": [], "z": []}
        for key in keys:
            fn = os.path.join(DST, str(lv), *map(str, key))
            assert os.path.getsize(fn) == L.CH ** 3, fn
            c = np.fromfile(fn, np.uint8).reshape(L.CH, L.CH, L.CH)
            o = (np.array(key) - self.lo) * L.CH
            self.a[o[0]:o[0] + L.CH, o[1]:o[1] + L.CH, o[2]:o[2] + L.CH] = c
            c = c.astype(np.int16)
            for ax, nm in ((2, "x"), (0, "z")):
                d = np.diff(c, axis=ax)
                both = (np.delete(c, -1, axis=ax) > 0) & (np.delete(c, 0, axis=ax) > 0)
                diffs[nm].append(d[both].ravel())
        self.sig_vox = {}
        for nm, v in diffs.items():
            v = np.concatenate(v).astype(np.float64)
            self.sig_vox[nm] = float(1.4826 * np.median(np.abs(v - np.median(v))) / np.sqrt(2))
        self.sigma_vox = max(self.sig_vox.values())
        self.lo_vox = self.lo * L.CH

    def sample(self, p0):
        p = np.asarray(p0, float).reshape(-1, 3) / (2 ** self.lv)
        c = np.stack([p[:, 2] - self.lo_vox[0], p[:, 1] - self.lo_vox[1], p[:, 0] - self.lo_vox[2]])
        return map_coordinates(self.a, c, order=1, mode="constant", cval=np.nan,
                               prefilter=False).reshape(np.shape(p0)[:-1])


class Box3:
    """Level-3 cal2 box with the same sample() interface."""

    def __init__(self, keys):
        arr, lo = L.load_ct(keys)
        self.lv, self.a, self.lo_vox = 3, arr, np.asarray(lo)

    def sample(self, p0):
        p = np.asarray(p0, float).reshape(-1, 3) / 8
        c = np.stack([p[:, 2] - self.lo_vox[0], p[:, 1] - self.lo_vox[1], p[:, 0] - self.lo_vox[2]])
        return map_coordinates(self.a, c, order=1, mode="constant", cval=np.nan,
                               prefilter=False).reshape(np.shape(p0)[:-1])


def frame(n, d):
    t1 = d - np.dot(d, n) * n
    t1 /= np.linalg.norm(t1)
    return t1, np.cross(n, t1)


def profile(box, base, n, t1, t2):
    """Mean over the k x k in-plane square (one level-3 voxel wide) at each s; also d(s) = half1 - half2."""
    k = 2 ** (3 - box.lv)
    vox = L.UM * 2 ** box.lv
    o = (np.arange(k) - (k - 1) / 2) * vox
    q = (base[None, None, None, :] + (S[:, None, None, None] * n + o[None, :, None, None] * t1
                                      + o[None, None, :, None] * t2) / L.UM)
    v = box.sample(q)                                   # (len(S), k, k)
    prof = v.mean(axis=(1, 2))
    d = None if k == 1 else v[:, :k // 2].mean(axis=(1, 2)) - v[:, k // 2:].mean(axis=(1, 2))
    return prof, d


def ridges(prof, cut):
    idx, pr = find_peaks(prof, prominence=(None, None))
    keep = pr["prominences"] >= cut
    return S[idx[keep]], pr["prominences"][keep]


def classify(rs, s34, s35):
    if s34 is None or s35 is None:
        return "other", None, None, "no w034 or w035 crossing"
    lo_w, hi_w = s34 - W_PAD, s35 + W_PAD
    mid = (s34 + s35) / 2
    w = rs[(rs >= lo_w) & (rs <= hi_w)]
    up = [r for r in w if abs(r - s35) <= UP_TOL]
    if up:
        r_up = min(up, key=lambda r: abs(r - s35))
        low = [r for r in w if r <= r_up - SEP_MIN and r <= mid]
        if low:
            r_lo = min(low, key=lambda r: abs(r - s34))
            return "two sheets (A-like)", float(r_up), float(r_lo), f"{len(w)} ridges in W0"
    if len(w) == 1 and w[0] > mid and abs(w[0] - s34) > LO_TOL:
        return "one sheet (B-like)", float(w[0]), None, "1 ridge in W0"
    return "other", None, None, f"{len(w)} ridges in W0: {[round(float(x), 1) for x in w]}"


def call(classes):
    a = sum(c.startswith("two") for c in classes)
    b = sum(c.startswith("one") for c in classes)
    verdict = "A" if a >= 5 and b <= 1 else "B" if b >= 5 and a <= 1 else "undecided"
    return verdict, {"A-like": a, "B-like": b, "other": len(classes) - a - b}


def follow(profs, j0, s0, cut):
    cur, path, lost = s0, [s0], 0
    for j in range(j0 + 1, len(profs)):
        rs, _ = ridges(profs[j], cut)
        c = rs[np.abs(rs - cur) <= REACH]
        if c.size:
            cur = float(c[np.argmin(np.abs(c - cur))])
        else:
            win = np.abs(S - cur) <= REACH
            cur = float(S[win][np.argmax(profs[j][win])])
            lost += 1
        path.append(cur)
    return path, lost


def font(sz):
    try:
        return ImageFont.load_default(size=sz)
    except TypeError:
        return ImageFont.load_default()


def main():
    rule_sha = hashlib.sha256(open(RULE, "rb").read()).hexdigest()
    p, keys3 = L.load_plan()
    meshes, normal_at = L.load_meshes(p)
    P, kept, _ = L.load_trace("real")
    assert kept[0:3, 10:19].all()
    est = json.load(open(os.path.join(L.OUT, "finer_ct_estimate.json")))
    boxes = {lv: Box(lv, [tuple(k) for k in est["levels"][str(lv)]["keys_zyx"]]) for lv in (1, 2)}
    boxes[3] = Box3(keys3)
    prof_json = json.load(open(os.path.join(L.OUT, "profiles.json")))
    tmax = 200.0 / L.UM

    def cross(w, b, n):
        t = meshes[w].cross(b, n, tmax)
        return round(float(t[0] * L.UM), 1) if t else None

    geo = {}
    for r in ROWS:
        for c in COLS:
            b = P[r, c]
            n = normal_at(b[None])[0]
            t1, t2 = frame(n, P[r, c + 1] - P[r, c - 1])
            geo[(r, c)] = (b, n, t1, t2, cross("w034", b, n), cross("w035", b, n))
    out = {"rule_file": "out/cal2_look/FINER-RULE.md", "rule_sha256": rule_sha, "levels": {}}
    profs_all = {}
    for lv in (1, 2, 3):
        box = boxes[lv]
        profs, ds = {}, []
        for key, (b, n, t1, t2, _, _) in geo.items():
            pr, d = profile(box, b, n, t1, t2)
            profs[key] = pr
            if d is not None:
                ds.append(d[np.isfinite(d)])
        profs_all[lv] = profs
        if lv == 3:
            sig = {"sigma_prof": round(float(prof_json["noise"]["sigma_voxel"]["z"]), 2),
                   "source": "cal2_look.py (profiles.json noise)"}
            cut = float(prof_json["cuts"]["noise3"])
        else:
            dd = np.concatenate(ds)
            s_half = float(1.4826 * np.median(np.abs(dd - np.median(dd))) / 2)
            k = 2 ** (3 - lv)
            s_white = box.sigma_vox / k
            sp = max(s_half, s_white)
            sig = {"sigma_vox_x": round(box.sig_vox["x"], 2), "sigma_vox_z": round(box.sig_vox["z"], 2),
                   "sigma_half": round(s_half, 2), "sigma_white": round(s_white, 2), "sigma_prof": round(sp, 2),
                   "k": k}
            cut = 3 * np.sqrt(2) * sp
        variants = [("main", cut)] + ([("2x_cut", 2 * cut)] if lv == 1 else [])
        lev = {"noise": sig, "cut": round(float(cut), 2), "columns": {}, "calls": {}}
        for vname, cv in variants:
            classes = []
            for c in COLS:
                b, n, t1, t2, s34, s35 = geo[(0, c)]
                pr = profs[(0, c)]
                if not np.all(np.isfinite(pr)):
                    cls, r_up, r_lo, why = "other", None, None, "not checked: sample outside fetched chunks"
                    rs, prom = np.array([]), np.array([])
                else:
                    rs, prom = ridges(pr, cv)
                    cls, r_up, r_lo, why = classify(rs, s34, s35)
                classes.append(cls)
                if vname == "main":
                    ent = {}
                    for r in ROWS:
                        _, _, _, _, a34, a35 = geo[(r, c)]
                        prr = profs[(r, c)]
                        ok = bool(np.all(np.isfinite(prr)))
                        rr, pp = ridges(prr, cv) if ok else (np.array([]), np.array([]))
                        ent[f"row{r}"] = {"s34": a34, "s35": a35, "checked": ok,
                                          "ridges_s_um": [round(float(x), 1) for x in rr],
                                          "prominences": [round(float(x), 1) for x in pp]}
                    ent["row0_class"] = cls
                    ent["row0_R_up"] = r_up
                    ent["row0_R_lo"] = r_lo
                    ent["row0_why"] = why
                    lev["columns"][str(c)] = ent
            verdict, counts = call(classes)
            lev["calls"][vname] = {"verdict": verdict, "counts": counts, "classes": classes,
                                   "cut": round(float(cv), 2)}
            print(f"level {lv} {vname} cut {cv:.2f}: {verdict} {counts} {classes}")
        lev["noise"]["cut_main"] = round(float(cut), 2)
        print(f"level {lv} noise {sig}")
        out["levels"][str(lv)] = lev
    out["decision"] = {"level": 1, **out["levels"]["1"]["calls"]["main"]}

    # rows 1 and 2: followers on level 1 for A-like columns
    box, cut = boxes[1], out["levels"]["1"]["cut"]
    fol = {}
    rib = {}
    for c in COLS:
        col = out["levels"]["1"]["columns"][str(c)]
        pts = P[0:3, c]
        seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        cum = np.r_[0.0, np.cumsum(seg)]
        u = np.unique(np.r_[np.arange(0, cum[-1], STEP / L.UM), cum])
        jr = [int(np.argmin(np.abs(u - x))) for x in cum]
        profs, base_l, s34_l, s35_l = [], [], [], []
        for uu in u:
            i = min(int(np.searchsorted(cum, uu, side="right") - 1), len(cum) - 2)
            f = (uu - cum[i]) / seg[i]
            b = (1 - f) * pts[i] + f * pts[i + 1]
            nb = [(1 - f) * P[i, cc] + f * P[i + 1, cc] for cc in (c - 1, c + 1)]
            n = normal_at(b[None])[0]
            t1, t2 = frame(n, nb[1] - nb[0])
            pr, _ = profile(box, b, n, t1, t2)
            profs.append(pr)
            base_l.append(b)
            s34_l.append(cross("w034", b, n))
            s35_l.append(cross("w035", b, n))
        rib[c] = (u, jr, profs, s34_l, s35_l)
        if not col["row0_class"].startswith("two"):
            fol[str(c)] = {"applies": False, "row0_class": col["row0_class"]}
            continue
        if not all(np.all(np.isfinite(x)) for x in profs):
            fol[str(c)] = {"applies": True, "result": "not checked: path leaves fetched chunks"}
            continue
        up, l_up = follow(profs, 0, col["row0_R_up"], cut)
        lo, l_lo = follow(profs, 0, col["row0_R_lo"], cut)
        ent = {"applies": True, "lost_up": l_up, "lost_lo": l_lo, "steps": len(u) - 1,
               "path_up": [round(x, 1) for x in up], "path_lo": [round(x, 1) for x in lo],
               "u_um": [round(float(x * L.UM), 1) for x in u], "row_index": jr}
        for r in (1, 2):
            j = jr[r]
            pr = profs[j]
            rs, _ = ridges(pr, cut)
            fl = float(np.percentile(pr, 10))
            i0 = int(np.argmin(np.abs(S)))
            if rs.size == 0:
                lab, tr = "off sheet", None
            else:
                tr = float(rs[np.argmin(np.abs(rs))])
                kr = int(np.argmin(np.abs(S - tr)))
                on = (pr[i0] - fl) >= 0.5 * (pr[kr] - fl)
                if not on:
                    lab = "off sheet"
                else:
                    a, b_ = abs(up[j] - tr) <= MATCH, abs(lo[j] - tr) <= MATCH
                    lab = ("merged" if a and b_ else "on w035's side" if a else "on w034's side" if b_
                           else "neither")
            ent[f"row{r}"] = {"label": lab, "trace_ridge_s_um": tr, "up_s_um": round(up[j], 1),
                              "lo_s_um": round(lo[j], 1), "s34": s34_l[j], "s35": s35_l[j]}
        fol[str(c)] = ent
        print(f"col {c}: row1 {ent['row1']['label']} row2 {ent['row2']['label']} lost up/lo {l_up}/{l_lo}")
    out["rows_1_2_level1"] = fol
    tally = {}
    for r in (1, 2):
        t = {}
        for v in fol.values():
            if v.get("applies") and f"row{r}" in v:
                t[v[f"row{r}"]["label"]] = t.get(v[f"row{r}"]["label"], 0) + 1
        tally[f"row{r}"] = t
    out["rows_1_2_tally"] = tally
    print("rows 1-2 tally", tally)
    json.dump(out, open(os.path.join(L.OUT, "finer_look.json"), "w"), indent=1)

    # figure: column 15, level 1, rows 0 to 2
    c = 15
    u, jr, profs, s34_l, s35_l = rib[c]
    img = np.array(profs)                                    # (len(u), len(S))
    lo_g, hi_g = np.nanpercentile(img, 1), np.nanpercentile(img, 99)
    g = np.clip((img.T[::-1] - lo_g) / (hi_g - lo_g), 0, 1) * 255
    g = np.where(np.isnan(img.T[::-1]), 90, g).astype(np.uint8)
    UPX, UPY = 10, 3
    H, W = g.shape
    im = Image.fromarray(np.repeat(g[..., None], 3, -1)).resize((W * UPX, H * UPY), Image.NEAREST)
    dr = ImageDraw.Draw(im)

    def yy(s):
        return (S[-1] - s) * UPY + UPY / 2

    for j in range(W):
        for v, col_ in ((s34_l[j], (0, 200, 255)), (s35_l[j], (255, 220, 0))):
            if v is not None and S[0] <= v <= S[-1]:
                dr.rectangle((j * UPX, yy(v) - 1, j * UPX + UPX - 1, yy(v)), fill=col_)
    f = fol.get(str(c), {})
    for key, col_ in (("path_up", (255, 140, 0)), ("path_lo", (40, 220, 40))):
        for j, v in enumerate(f.get(key) or []):
            dr.rectangle((j * UPX + 2, yy(v) - 1, j * UPX + 3, yy(v) + 1), fill=col_)
    for r, j in enumerate(jr):
        x, y = j * UPX + UPX / 2, yy(0.0)
        dr.line((x - 7, y, x + 7, y), fill=(255, 40, 40), width=2)
        dr.line((x, y - 7, x, y + 7), fill=(255, 40, 40), width=2)
        dr.text((x + 5, y - 20), f"s{r}", fill=(255, 40, 40), font=font(13))
    pad_top, pad_bot = 46, 30
    canvas = Image.new("RGB", (max(im.width + 40, 900), im.height + pad_top + pad_bot), (255, 255, 255))
    canvas.paste(im, (0, pad_top))
    d2 = ImageDraw.Draw(canvas)
    d2.text((6, 4), f"cal2, level-1 CT (2.258 um voxels; mean over a 9.03 um in-plane square), column {c}, "
                    f"steps 0 to 2, s -130 to +60 um (up = outward)", fill=(0, 0, 0), font=font(14))
    xl = 6
    for col_, lab in (((0, 200, 255), "w034"), ((255, 220, 0), "w035"), ((255, 40, 40), "trace (+)"),
                      ((255, 140, 0), "follower from R_up"), ((40, 220, 40), "follower from R_lo")):
        d2.rectangle((xl, 28, xl + 10, 38), fill=col_)
        d2.text((xl + 14, 25), lab, fill=(0, 0, 0), font=font(13))
        xl += 24 + 7 * len(lab)
    y0 = pad_top + im.height + 14
    d2.line((10, y0, 10 + 50 * UPY, y0), fill=(0, 0, 0), width=3)
    d2.text((16 + 50 * UPY, y0 - 8), f"50 um along the normal (vertical); horizontal {STEP:.2f} um per "
                                      f"{UPX} px", fill=(0, 0, 0), font=font(12))
    fn = os.path.join(L.OUT, "png", "cal2_col15_L1.png")
    canvas.save(fn)
    print(fn, canvas.size, "grey window", round(float(lo_g), 1), round(float(hi_g), 1))


if __name__ == "__main__":
    main()
