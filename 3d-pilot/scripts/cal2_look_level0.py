"""Apply FINER-RULE.md unchanged at level 0 (LEVEL0-ADDENDUM.md), next to levels 1, 2 and 3.

Level 0: the 48 chunks of level0_ct_estimate.json in downloads/ct_level0/0/z/y/x. Levels 1 and 2: the 14 chunks of
finer_ct_estimate.json. Level 3: cal2's box. The rule's functions (Box, profile, ridges, classify, call, follow,
frame) are imported from cal2_look_finer.py, not copied; the per-level loop and the follower block repeat
cal2_look_finer.main for any level. Levels 1 to 3 are recomputed and checked against finer_look.json and
finer_width.json, so the level-0 numbers come from code that reproduces the earlier ones.
Outside the rule: half-height width and 10-90% edge (as cal2_look_finer_width.py), and the dip check of the addendum.
Writes out/cal2_look/level0_look.json and out/cal2_look/png/cal2_col15_L0.png.
Usage (WSL venv): python -I cal2_look_level0.py
"""
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cal2_look as L  # noqa: E402
import cal2_look_finer as F  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image, ImageDraw  # noqa: E402
from scipy.signal import find_peaks  # noqa: E402

RULE_SHA = "c7d471962f91e3f6e36e42a18ba1bf3efe3760f874ba55c13f58a23b0b2feea8"
ADD = os.path.join(L.OUT, "LEVEL0-ADDENDUM.md")
NEW_DST = os.path.join(L.ROOT, "downloads", "ct_level0")
LEVELS = (0, 1, 2, 3)


def sha(fn):
    return hashlib.sha256(open(fn, "rb").read()).hexdigest()


def box_at(lv, keys, dst):
    old = F.DST
    F.DST = dst
    try:
        return F.Box(lv, keys)
    finally:
        F.DST = old


def width(pr, rs_w):
    """cal2_look_finer_width.py's measure for one profile whose W0 holds exactly one ridge at rs_w."""
    k0 = int(np.argmin(np.abs(F.S - rs_w)))
    f10 = float(np.percentile(pr, 10))
    top = pr[k0] - f10
    a = k0
    while a > 0 and pr[a - 1] - f10 >= 0.5 * top:
        a -= 1
    z = k0
    while z < len(pr) - 1 and pr[z + 1] - f10 >= 0.5 * top:
        z += 1

    def lvl(x):
        return round(float((pr[int(np.argmin(np.abs(F.S - x)))] - f10) / top), 2)

    def rise(step):
        j, j90, j10 = k0, None, None
        while 0 <= j < len(pr):
            v = (pr[j] - f10) / top
            if j90 is None and v <= 0.9:
                j90 = j
            if v <= 0.1:
                j10 = j
                break
            j += step
        return None if j90 is None or j10 is None else float(abs(F.S[j10] - F.S[j90]))

    return {"crest_s_um": float(rs_w), "half_width_um": float(F.S[z] - F.S[a]),
            "edge_10_90_um_inner_outer": [rise(-1), rise(+1)], "half_lo_hi_um": [float(F.S[a]), float(F.S[z])]}, lvl


def dip_check(pr, s34, s35, cut):
    """Every local maximum inside W0 and between the crossings, with its prominence (no cut)."""
    idx, pp = find_peaks(pr, prominence=(None, None))
    s, prom = F.S[idx], pp["prominences"]
    lo_w, hi_w = s34 - F.W_PAD, s35 + F.W_PAD
    inw = (s >= lo_w) & (s <= hi_w)
    a, b = min(s34, s35), max(s34, s35)
    btw = (s >= a) & (s <= b)
    ws, wp = s[inw], prom[inw]
    main = int(np.argmax(wp)) if wp.size else None
    others = [float(x) for i, x in enumerate(wp) if i != main]
    second = max(others) if others else 0.0
    return {"maxima_in_W0": [[round(float(x), 1), round(float(y), 2)] for x, y in zip(ws, wp)],
            "maxima_between_crossings": [[round(float(x), 1), round(float(y), 2)] for x, y in zip(s[btw], prom[btw])],
            "main_s_um": None if main is None else float(ws[main]),
            "main_prominence": None if main is None else round(float(wp[main]), 2),
            "largest_other_prominence_in_W0": round(second, 2),
            "largest_other_over_cut": round(second / cut, 3),
            "n_other_maxima_in_W0": len(others),
            "n_other_over_cut": int(sum(x >= cut for x in others))}


def main():
    assert sha(os.path.join(L.OUT, "FINER-RULE.md")) == RULE_SHA, "FINER-RULE.md changed"
    add_rec = json.load(open(os.path.join(L.ROOT, "logs", "level0_addendum_sha256.json")))
    assert sha(ADD) == add_rec["sha256"], "LEVEL0-ADDENDUM.md changed after its hash was saved"
    p, keys3 = L.load_plan()
    meshes, normal_at = L.load_meshes(p)
    P, kept, _ = L.load_trace("real")
    assert kept[0:3, 10:19].all()
    est = json.load(open(os.path.join(L.OUT, "finer_ct_estimate.json")))
    est0 = json.load(open(os.path.join(L.OUT, "level0_ct_estimate.json")))
    keys0 = [tuple(k) for k in est0["level"]["keys_zyx"]]
    on_disk = [k for k in keys0 if os.path.exists(os.path.join(NEW_DST, "0", *map(str, k)))]
    assert len(on_disk) == len(keys0), f"{len(keys0) - len(on_disk)} level-0 chunks missing"
    boxes = {0: box_at(0, keys0, NEW_DST)}
    for lv in (1, 2):
        boxes[lv] = box_at(lv, [tuple(k) for k in est["levels"][str(lv)]["keys_zyx"]], F.DST)
    boxes[3] = F.Box3(keys3)
    prof_json = json.load(open(os.path.join(L.OUT, "profiles.json")))
    old_look = json.load(open(os.path.join(L.OUT, "finer_look.json")))
    old_width = json.load(open(os.path.join(L.OUT, "finer_width.json")))
    tmax = 200.0 / L.UM

    def cross(w, b, n):
        t = meshes[w].cross(b, n, tmax)
        return round(float(t[0] * L.UM), 1) if t else None

    geo = {}
    for r in F.ROWS:
        for c in F.COLS:
            b = P[r, c]
            n = normal_at(b[None])[0]
            t1, t2 = F.frame(n, P[r, c + 1] - P[r, c - 1])
            geo[(r, c)] = (b, n, t1, t2, cross("w034", b, n), cross("w035", b, n))
    out = {"rule_file": "out/cal2_look/FINER-RULE.md", "rule_sha256": RULE_SHA,
           "addendum_file": "out/cal2_look/LEVEL0-ADDENDUM.md", "addendum_sha256": add_rec["sha256"],
           "addendum_saved": add_rec["saved"], "levels": {}}
    profs_all = {}
    for lv in LEVELS:
        box = boxes[lv]
        profs, ds = {}, []
        for key, (b, n, t1, t2, _, _) in geo.items():
            pr, d = F.profile(box, b, n, t1, t2)
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
        variants = [("main", cut)] + ([("2x_cut", 2 * cut)] if lv in (0, 1) else [])
        lev = {"noise": sig, "cut": round(float(cut), 2), "columns": {}, "calls": {}}
        for vname, cv in variants:
            classes = []
            for c in F.COLS:
                b, n, t1, t2, s34, s35 = geo[(0, c)]
                pr = profs[(0, c)]
                if not np.all(np.isfinite(pr)):
                    cls, r_up, r_lo, why = "other", None, None, "not checked: sample outside fetched chunks"
                else:
                    rs, prom = F.ridges(pr, cv)
                    cls, r_up, r_lo, why = F.classify(rs, s34, s35)
                classes.append(cls)
                if vname == "main":
                    ent = {}
                    for r in F.ROWS:
                        _, _, _, _, a34, a35 = geo[(r, c)]
                        prr = profs[(r, c)]
                        ok = bool(np.all(np.isfinite(prr)))
                        rr, pp = F.ridges(prr, cv) if ok else (np.array([]), np.array([]))
                        ent[f"row{r}"] = {"s34": a34, "s35": a35, "checked": ok,
                                          "ridges_s_um": [round(float(x), 1) for x in rr],
                                          "prominences": [round(float(x), 1) for x in pp]}
                    ent["row0_class"] = cls
                    ent["row0_R_up"] = r_up
                    ent["row0_R_lo"] = r_lo
                    ent["row0_why"] = why
                    lev["columns"][str(c)] = ent
            verdict, counts = F.call(classes)
            lev["calls"][vname] = {"verdict": verdict, "counts": counts, "classes": classes,
                                   "cut": round(float(cv), 2)}
            print(f"level {lv} {vname} cut {cv:.2f}: {verdict} {counts} {classes}")
        lev["noise"]["cut_main"] = round(float(cut), 2)
        print(f"level {lv} noise {sig}")
        out["levels"][str(lv)] = lev

    # reproduction check: levels 1 to 3 against finer_look.json
    rep = {}
    for lv in (1, 2, 3):
        a, b = out["levels"][str(lv)], old_look["levels"][str(lv)]
        same = (a["cut"] == b["cut"] and a["noise"] == b["noise"] and a["calls"] == b["calls"]
                and a["columns"] == b["columns"])
        rep[str(lv)] = same
    out["reproduces_finer_look_json"] = rep
    print("reproduces finer_look.json:", rep)

    # rows 1 and 2: followers at level 0 for level-0 A-like columns (cal2_look_finer.main's block, level 0)
    box, cut = boxes[0], out["levels"]["0"]["cut"]
    fol, rib = {}, {}
    for c in F.COLS:
        col = out["levels"]["0"]["columns"][str(c)]
        pts = P[0:3, c]
        seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        cum = np.r_[0.0, np.cumsum(seg)]
        u = np.unique(np.r_[np.arange(0, cum[-1], F.STEP / L.UM), cum])
        jr = [int(np.argmin(np.abs(u - x))) for x in cum]
        profs, s34_l, s35_l = [], [], []
        for uu in u:
            i = min(int(np.searchsorted(cum, uu, side="right") - 1), len(cum) - 2)
            f = (uu - cum[i]) / seg[i]
            b = (1 - f) * pts[i] + f * pts[i + 1]
            nb = [(1 - f) * P[i, cc] + f * P[i + 1, cc] for cc in (c - 1, c + 1)]
            n = normal_at(b[None])[0]
            t1, t2 = F.frame(n, nb[1] - nb[0])
            pr, _ = F.profile(box, b, n, t1, t2)
            profs.append(pr)
            s34_l.append(cross("w034", b, n))
            s35_l.append(cross("w035", b, n))
        rib[c] = (u, jr, profs, s34_l, s35_l)
        path_ok = all(np.all(np.isfinite(x)) for x in profs)
        if not col["row0_class"].startswith("two"):
            fol[str(c)] = {"applies": False, "row0_class": col["row0_class"], "path_all_checked": path_ok}
            continue
        if not path_ok:
            fol[str(c)] = {"applies": True, "result": "not checked: path leaves fetched chunks"}
            continue
        up, l_up = F.follow(profs, 0, col["row0_R_up"], cut)
        lo, l_lo = F.follow(profs, 0, col["row0_R_lo"], cut)
        ent = {"applies": True, "lost_up": l_up, "lost_lo": l_lo, "steps": len(u) - 1,
               "path_up": [round(x, 1) for x in up], "path_lo": [round(x, 1) for x in lo],
               "u_um": [round(float(x * L.UM), 1) for x in u], "row_index": jr}
        for r in (1, 2):
            j = jr[r]
            pr = profs[j]
            rs, _ = F.ridges(pr, cut)
            fl = float(np.percentile(pr, 10))
            i0 = int(np.argmin(np.abs(F.S)))
            if rs.size == 0:
                lab, tr = "off sheet", None
            else:
                tr = float(rs[np.argmin(np.abs(rs))])
                kr = int(np.argmin(np.abs(F.S - tr)))
                on = (pr[i0] - fl) >= 0.5 * (pr[kr] - fl)
                if not on:
                    lab = "off sheet"
                else:
                    a_, b_ = abs(up[j] - tr) <= F.MATCH, abs(lo[j] - tr) <= F.MATCH
                    lab = ("merged" if a_ and b_ else "on w035's side" if a_ else "on w034's side" if b_
                           else "neither")
            ent[f"row{r}"] = {"label": lab, "trace_ridge_s_um": tr, "up_s_um": round(up[j], 1),
                              "lo_s_um": round(lo[j], 1), "s34": s34_l[j], "s35": s35_l[j]}
        fol[str(c)] = ent
        print(f"col {c}: row1 {ent['row1']['label']} row2 {ent['row2']['label']} lost up/lo {l_up}/{l_lo}")
    out["rows_1_2_level0"] = fol
    tally = {}
    for r in (1, 2):
        t = {}
        for v in fol.values():
            if v.get("applies") and f"row{r}" in v:
                t[v[f"row{r}"]["label"]] = t.get(v[f"row{r}"]["label"], 0) + 1
        tally[f"row{r}"] = t
    out["rows_1_2_tally_level0"] = tally
    print("rows 1-2 tally (level 0)", tally)

    # descriptive, outside the rule: width, edge, levels at the crossings, dip check (row 0)
    desc, wrep = {}, {}
    for c in F.COLS:
        b, n, t1, t2, s34, s35 = geo[(0, c)]
        ent = {}
        for lv in LEVELS:
            pr = profs_all[lv][(0, c)]
            cutl = out["levels"][str(lv)]["cut"]
            rs, _ = F.ridges(pr, cutl)
            w = rs[(rs >= s34 - F.W_PAD) & (rs <= s35 + F.W_PAD)]
            e = {"dip_check": dip_check(pr, s34, s35, cutl), "n_ridges_in_W0": int(w.size)}
            if w.size == 1:
                wd, lvl = width(pr, w[0])
                wd["level_at_w034"], wd["level_at_w035"] = lvl(s34), lvl(s35)
                e.update(wd)
            ent[f"L{lv}"] = e
            if lv in (1, 2, 3):
                o = old_width["columns"][str(c)][f"L{lv}"]
                mine = {k: e.get(k) for k in o} if o else None
                wrep[f"{c}_L{lv}"] = (mine == o)
        desc[str(c)] = ent
        l0 = ent["L0"]
        print(c, "L0 crest", l0.get("crest_s_um"), "width", l0.get("half_width_um"), "edge",
              l0.get("edge_10_90_um_inner_outer"), "levels", l0.get("level_at_w034"), l0.get("level_at_w035"),
              "| 2nd max in W0", l0["dip_check"]["largest_other_prominence_in_W0"],
              f"({l0['dip_check']['largest_other_over_cut']} x cut)",
              "| maxima between crossings", l0["dip_check"]["maxima_between_crossings"])
    out["descriptive"] = {"what": "outside FINER-RULE.md: width and edge as cal2_look_finer_width.py; dip check as "
                                  "LEVEL0-ADDENDUM.md", "columns": desc}
    for lv in LEVELS:
        ws = [e[f"L{lv}"]["half_width_um"] for e in desc.values() if "half_width_um" in e[f"L{lv}"]]
        eg = [x for e in desc.values() if "half_width_um" in e[f"L{lv}"]
              for x in e[f"L{lv}"]["edge_10_90_um_inner_outer"] if x is not None]
        out["descriptive"][f"median_half_width_L{lv}_um"] = float(np.median(ws)) if ws else None
        out["descriptive"][f"n_columns_L{lv}"] = len(ws)
        out["descriptive"][f"median_edge_10_90_L{lv}_um"] = float(np.median(eg)) if eg else None
        out["descriptive"][f"n_flanks_L{lv}"] = len(eg)
        print(f"L{lv} median half width {out['descriptive'][f'median_half_width_L{lv}_um']} um over {len(ws)} "
              f"columns; median 10-90 edge {out['descriptive'][f'median_edge_10_90_L{lv}_um']} um over {len(eg)} "
              f"flanks")
    out["reproduces_finer_width_json"] = wrep
    print("reproduces finer_width.json:", all(wrep.values()), [k for k, v in wrep.items() if not v])
    out["all_level0_profiles_checked"] = bool(all(np.all(np.isfinite(v)) for v in profs_all[0].values()))
    json.dump(out, open(os.path.join(L.OUT, "level0_look.json"), "w"), indent=1)

    # figure: column 15, level 0, rows 0 to 2 (drawn as cal2_col15_L1.png)
    c = 15
    u, jr, profs, s34_l, s35_l = rib[c]
    img = np.array(profs)
    lo_g, hi_g = np.nanpercentile(img, 1), np.nanpercentile(img, 99)
    g = np.clip((img.T[::-1] - lo_g) / (hi_g - lo_g), 0, 1) * 255
    g = np.where(np.isnan(img.T[::-1]), 90, g).astype(np.uint8)
    UPX, UPY = 10, 3
    H, W = g.shape
    im = Image.fromarray(np.repeat(g[..., None], 3, -1)).resize((W * UPX, H * UPY), Image.NEAREST)
    dr = ImageDraw.Draw(im)

    def yy(s):
        return (F.S[-1] - s) * UPY + UPY / 2

    for j in range(W):
        for v, col_ in ((s34_l[j], (0, 200, 255)), (s35_l[j], (255, 220, 0))):
            if v is not None and F.S[0] <= v <= F.S[-1]:
                dr.rectangle((j * UPX, yy(v) - 1, j * UPX + UPX - 1, yy(v)), fill=col_)
    f = fol.get(str(c), {})
    for key, col_ in (("path_up", (255, 140, 0)), ("path_lo", (40, 220, 40))):
        for j, v in enumerate(f.get(key) or []):
            dr.rectangle((j * UPX + 2, yy(v) - 1, j * UPX + 3, yy(v) + 1), fill=col_)
    for r, j in enumerate(jr):
        x, y = j * UPX + UPX / 2, yy(0.0)
        dr.line((x - 7, y, x + 7, y), fill=(255, 40, 40), width=2)
        dr.line((x, y - 7, x, y + 7), fill=(255, 40, 40), width=2)
        dr.text((x + 5, y - 20), f"s{r}", fill=(255, 40, 40), font=F.font(13))
    pad_top, pad_bot = 46, 30
    canvas = Image.new("RGB", (max(im.width + 40, 900), im.height + pad_top + pad_bot), (255, 255, 255))
    canvas.paste(im, (0, pad_top))
    d2 = ImageDraw.Draw(canvas)
    d2.text((6, 4), f"cal2, level-0 CT (1.129 um voxels; mean over a 9.03 um in-plane square), column {c}, "
                    f"steps 0 to 2, s -130 to +60 um (up = outward)", fill=(0, 0, 0), font=F.font(14))
    xl = 6
    for col_, lab in (((0, 200, 255), "w034"), ((255, 220, 0), "w035"), ((255, 40, 40), "trace (+)"),
                      ((255, 140, 0), "follower from R_up"), ((40, 220, 40), "follower from R_lo")):
        d2.rectangle((xl, 28, xl + 10, 38), fill=col_)
        d2.text((xl + 14, 25), lab, fill=(0, 0, 0), font=F.font(13))
        xl += 24 + 7 * len(lab)
    y0 = pad_top + im.height + 14
    d2.line((10, y0, 10 + 50 * UPY, y0), fill=(0, 0, 0), width=3)
    d2.text((16 + 50 * UPY, y0 - 8), f"50 um along the normal (vertical); horizontal {F.STEP:.2f} um per "
                                      f"{UPX} px", fill=(0, 0, 0), font=F.font(12))
    fn = os.path.join(L.OUT, "png", "cal2_col15_L0.png")
    canvas.save(fn)
    print(fn, canvas.size, "grey window", round(float(lo_g), 1), round(float(hi_g), 1))


if __name__ == "__main__":
    main()
