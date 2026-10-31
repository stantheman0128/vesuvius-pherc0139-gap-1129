"""Does the seed sit where two sheets touch? A check on cal2_look_ridge.py's verdicts.

The ridge verdicts start backward followers only at the last row (row 8). If the seed's ridge is two sheets in
contact (w034's and w035's) that part after a few steps, then both strands reach the seed's ridge when followed back,
and "only the trace's ridge reaches the seed" at row 8 says nothing about which strand is the seed's sheet.

Per trace column (same ribbon as cal2_look_ridge.py: CT along the local w035 normal, one sample per level-3 voxel),
for every row r = 1..8 start backward followers (cal2_look_ridge.follow, cut noise3, slope 1.0 and 0.5) from the ridge
nearest the trace (s = 0), the ridge nearest the w035 crossing and the ridge nearest the w034 crossing, each only when
that point is on a sheet (inside the half-maximum of its nearest ridge, as in cal2_look.py).
Record whether each reaches the seed's ridge (back at row 0, within 9.03 um of it).
Per column: the first row where the trace and w035 sit on different ridges (the split row) and, at that row and the
next two, whether both ridges reach the seed's ridge ("touch": the CT does not say which strand is the seed's sheet),
only the trace's ("trace keeps the seed's ridge"), only w035's, or neither.
Also the half-level width of the seed's ridge at row 0 (compare with the ridge widths in profiles.json ct_ridges).
Usage (WSL venv): python -I cal2_look_split.py      Writes out/cal2_look/split_check.json
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cal2_look as L  # noqa: E402
import cal2_look_ridge as R  # noqa: E402
import numpy as np  # noqa: E402


def half_width(prof, k0):
    fl = float(np.percentile(prof, 10))
    top = prof[k0] - fl
    a = k0
    while a > 0 and (prof[a - 1] - fl) >= 0.5 * top:
        a -= 1
    b = k0
    while b < len(prof) - 1 and (prof[b + 1] - fl) >= 0.5 * top:
        b += 1
    return float(R.S[b] - R.S[a])


def main():
    p, keys = L.load_plan()
    arr, lo = L.load_ct(keys)
    ct = L.CT(arr, lo)
    meshes, normal_at = L.load_meshes(p)
    P, kept, _ = L.load_trace("real")
    rf = json.load(open(os.path.join(L.OUT, "ridge_follow.json")))
    cut = rf["cuts"]["noise3"]
    tmax = L.S_HALF / L.UM
    out = {"conventions": {"cut": "noise3", "cut_value": cut, "tol_um": R.TOL,
                           "reach": "backward follower ends at row 0 within tol of the seed's ridge",
                           "on_sheet": "level >= 0.5 of the nearest ridge (p10 floor), as in cal2_look.py"},
           "columns": {}, "summary": {}}
    tally = {sl: {} for sl in (1.0, 0.5)}
    widths = []
    for c in range(P.shape[1]):
        rows = [r for r in range(P.shape[0]) if kept[r, c]]
        u, cum, base, n, img = L.ribbon(ct, normal_at, P[rows, c], R.DU, 0.0, R.S)
        jr = [int(np.argmin(np.abs(u - x))) for x in cum]
        cr = {w: [] for w in ("w034", "w035")}
        for j in jr:
            for w in cr:
                t = meshes[w].cross(base[j], n[j], tmax)
                cr[w].append(t[0] * L.UM if t else None)
        col = {"rows": rows}
        f1, _, _, _ = R.follow(img, cut, 0, 0.0, +1, 1.0)
        seed_ridge = f1[0][1]
        k0 = int(np.argmin(np.abs(R.S - seed_ridge)))
        col["seed_ridge_s_um"] = round(seed_ridge, 1)
        col["seed_ridge_half_width_um"] = round(half_width(img[0], k0), 1)
        col["w034_at_seed_um"] = None if cr["w034"][0] is None else round(cr["w034"][0], 1)
        widths.append(col["seed_ridge_half_width_um"])
        for sl in (1.0, 0.5):
            per = []
            for k, (r, j) in enumerate(zip(rows, jr)):
                if r == 0:
                    continue
                prof = img[j]
                idx, _ = L.ridges(prof, cut)
                rs = R.S[idx]
                fl = float(np.percentile(prof, 10))
                ent = {"row": r}
                for tag, x in (("trace", 0.0), ("w035", cr["w035"][k]), ("w034", cr["w034"][k])):
                    if x is None or rs.size == 0:
                        ent[tag] = None
                        continue
                    rk = float(rs[np.argmin(np.abs(rs - x))])
                    lev = (prof[int(np.argmin(np.abs(R.S - x)))] - fl) / max(
                        prof[int(np.argmin(np.abs(R.S - rk)))] - fl, 1e-6)
                    on = bool(lev >= L.HALF)
                    d = {"s_um": round(x, 1), "ridge_s_um": round(rk, 1), "on_sheet": on}
                    if on:
                        fb, _, _, _ = R.follow(img, cut, j, rk, -1, sl)
                        d["reaches_seed"] = bool(fb[-1][0] == 0 and abs(fb[-1][1] - seed_ridge) <= R.TOL)
                        d["back_end_s_um"] = round(fb[-1][1], 1)
                    ent[tag] = d
                t, w = ent["trace"], ent["w035"]
                ent["trace_w035_split"] = bool(t and w and t["on_sheet"] and w["on_sheet"]
                                               and abs(t["ridge_s_um"] - w["ridge_s_um"]) > 0.5)
                per.append(ent)
            split = [e for e in per if e["trace_w035_split"]]
            if not split:
                cls = "no row with trace and w035 on two different ridges"
                srow = None
            else:
                srow = split[0]["row"]
                kinds = []
                for e in split[:3]:
                    a, b = e["trace"]["reaches_seed"], e["w035"]["reaches_seed"]
                    kinds.append("touch" if a and b else "trace only" if a else "w035 only" if b else "neither")
                cls = ("touch" if "touch" in kinds else
                       "trace only" if all(x == "trace only" for x in kinds) else
                       "w035 only" if all(x == "w035 only" for x in kinds) else "mixed: " + ",".join(kinds))
            w34 = [e["w034"]["reaches_seed"] for e in per if e["w034"] and e["w034"]["on_sheet"]
                   and e["row"] <= 4]

            def shares(tag):          # rows 6..8 where the mesh is on a sheet and on the trace's ridge
                return [e["row"] for e in per if e["row"] >= 6 and e[tag] and e[tag]["on_sheet"] and e["trace"]
                        and e["trace"]["on_sheet"] and abs(e[tag]["ridge_s_um"] - e["trace"]["ridge_s_um"]) <= 0.5]
            col[f"slope{sl}"] = {"per_row": per, "split_row": srow, "class": cls,
                                 "w034_rows1to4_on_sheet_reaching_seed": f"{sum(w34)}/{len(w34)}",
                                 "rows6to8_w034_on_trace_ridge": shares("w034"),
                                 "rows6to8_w035_on_trace_ridge": shares("w035"),
                                 "rows5to8_w035_on_sheet_reaching_seed": "%d/%d" % (
                                     sum(e["w035"]["reaches_seed"] for e in per if e["row"] >= 5 and e["w035"]
                                         and e["w035"]["on_sheet"]),
                                     sum(1 for e in per if e["row"] >= 5 and e["w035"] and e["w035"]["on_sheet"])),
                                 "rows5to8_trace_on_sheet_reaching_seed": "%d/%d" % (
                                     sum(e["trace"]["reaches_seed"] for e in per if e["row"] >= 5 and e["trace"]
                                         and e["trace"]["on_sheet"]),
                                     sum(1 for e in per if e["row"] >= 5 and e["trace"] and e["trace"]["on_sheet"]))}
            grp = "cols_11_21" if 11 <= c <= 21 else "other_cols"
            tally[sl].setdefault(grp, {})
            tally[sl][grp][cls] = tally[sl][grp].get(cls, 0) + 1
        out["columns"][str(c)] = col
        m = col["slope1.0"]
        print(f"col {c:2d}: seed ridge {col['seed_ridge_s_um']:6.1f} width {col['seed_ridge_half_width_um']:5.1f} "
              f"w034@seed {col['w034_at_seed_um']} | split row {m['split_row']} | {m['class']} | "
              f"w034 r1-4 reach {m['w034_rows1to4_on_sheet_reaching_seed']} | r6-8 on trace ridge: w034 "
              f"{m['rows6to8_w034_on_trace_ridge']} w035 {m['rows6to8_w035_on_trace_ridge']} | r5-8 reach: trace "
              f"{m['rows5to8_trace_on_sheet_reaching_seed']} w035 {m['rows5to8_w035_on_sheet_reaching_seed']}")
    out["summary"] = {f"slope{sl}": v for sl, v in tally.items()}
    w = np.array(widths)
    out["summary"]["seed_ridge_half_width_um"] = {
        "cols_11_21_median": round(float(np.median(w[11:22])), 1),
        "other_cols_median": round(float(np.median(np.r_[w[:11], w[22:]])), 1)}
    print(json.dumps(out["summary"], indent=1))
    json.dump(out, open(os.path.join(L.OUT, "split_check.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
