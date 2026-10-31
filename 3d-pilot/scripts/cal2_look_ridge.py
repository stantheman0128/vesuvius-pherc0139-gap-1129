"""CT-only continuity check for cal2: in each trace column's ribbon (cal2_look.ribbon: the CT along the local w035
normal at points every level-3 voxel along the column's path through steps 0 to 8), follow brightness ridges and
ask which ridge is continuous with the sheet under the seed. The followers move on the CT alone; the meshes and the
tracer only give the start and end points.

Followers step one level-3 voxel (9.03 um) along the column and move to the nearest ridge (cal2_look.ridges, same
cut) within SLOPE x 9.03 um + JIT of their last position. JIT = 3 um allows for the crest of a broad ridge moving
between neighbouring samples. On a sample with no ridge in reach a follower takes the brightest sample in reach
and counts the sample as "lost".
  F1 forward from the seed's ridge (the ridge nearest the seed point; the seed lies on w035 by construction);
  F2 backward from the ridge nearest the trace at its last row;
  F3 backward from the ridge nearest the w035 crossing at that row.
A start more than one level-3 voxel (9.03 um) from its ridge is recorded: then that point has no ridge of its own.
F2 or F3 "reaches the seed ridge" if it gets back to the seed sample within 9.03 um of the seed ridge.
Verdicts at the last row use "on a sheet" as in cal2_look.py (inside the half-maximum of the nearest ridge):
  w035 off any sheet, trace on one   -> and whether the trace's ridge reaches the seed's ridge;
  trace off any sheet, w035 on one   -> and whether w035's ridge reaches the seed's ridge;
  both on one ridge                  -> agree (and whether that ridge reaches the seed's ridge);
  two ridges: only the trace's reaches the seed's ridge -> w035 left the seed's sheet;
  two ridges: only w035's reaches the seed's ridge      -> the trace left the seed's sheet;
  two ridges, both reach             -> they touch on the way: the CT does not decide; neither -> undecided.
These labels assume the seed's ridge is one sheet. In columns 11 to 17 the w034 mesh is only 27 to 43 um from w035
at the seed and the CT shows one ridge there, which the w034 mesh also follows later (cal2_look_split.py). There,
"reaches the seed's ridge" does not say whether a path is on w035's sheet or w034's.
Run for SLOPE 1.0 and 0.5 and cuts noise3 and noise2 to see whether the answer depends on them.
Also, per column, how dark the CT is along each surface's own path (trace, w034, w035, w036) at every sample from
the seed to the last row: a papyrus sheet is continuous and bright, so a path that runs at gap level has left one.
Usage (WSL venv): python -I cal2_look_ridge.py      Writes out/cal2_look/ridge_follow.json
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cal2_look as L  # noqa: E402
import numpy as np  # noqa: E402

DU = 8.0                     # level-0 voxels between ribbon samples = one level-3 voxel
JIT = 3.0                    # um
S = np.arange(-L.S_HALF, L.S_HALF + L.DS / 2, L.DS)
TOL = L.L3UM


def follow(img, cut, j0, s0, step, slope):
    delta = slope * DU * L.UM + JIT
    J = img.shape[0]
    cur, path, amb, lost, snap = None, [], 0, 0, None
    for j in (range(j0, J) if step > 0 else range(j0, -1, -1)):
        prof = img[j]
        idx, _ = L.ridges(prof, cut)
        rs = S[idx]
        if cur is None:                       # start: snap to the nearest ridge
            k = int(np.argmin(np.abs(rs - s0)))
            cur = float(rs[k])
            snap = round(cur - s0, 1)
        else:
            c = np.nonzero(np.abs(rs - cur) <= delta)[0]
            if c.size:
                cur = float(rs[c[np.argmin(np.abs(rs[c] - cur))]])
                amb += int(c.size > 1)
            else:
                win = np.abs(S - cur) <= delta
                cur = float(S[win][np.argmax(prof[win])])
                lost += 1
        on = np.any(np.abs(rs - cur) <= 0.5)          # neighbour gap only while on a ridge (not on "lost" samples)
        other = rs[np.abs(rs - cur) > 0.5]
        path.append((j, cur, float(np.min(np.abs(other - cur))) if (on and other.size) else None))
        if abs(cur) > L.S_HALF - 5:
            break
    return path, snap, amb, lost


def main():
    p, keys = L.load_plan()
    arr, lo = L.load_ct(keys)
    ct = L.CT(arr, lo)
    meshes, normal_at = L.load_meshes(p)
    P, kept, _ = L.load_trace("real")
    sig = {}
    for ax, nm in ((2, "x"), (0, "z")):
        dd = np.diff(arr[:, 100:200, :], axis=ax).ravel()
        sig[nm] = float(1.4826 * np.median(np.abs(dd - np.median(dd))) / np.sqrt(2))
    sigma = max(sig.values())
    cuts = {"noise3": 3 * np.sqrt(2) * sigma, "noise2": 2 * np.sqrt(2) * sigma}
    tmax = L.S_HALF / L.UM
    res = {"conventions": {"s_um": "along the local w035 outward normal, 0 = the column's path (trace points at its "
                                    "rows); + outward", "du_um": round(DU * L.UM, 3), "jit_um": JIT, "tol_um": TOL,
                           "path": "[u_um, ridge s_um] for the main setting (noise3, slope 1.0), for the figure"},
           "cuts": {k: round(float(v), 2) for k, v in cuts.items()}, "columns": {}, "summary": {}}
    settings = [(cn, sl) for cn in cuts for sl in (1.0, 0.5)]
    main_set = ("noise3", 1.0)
    tallies = {f"{cn}_slope{sl}": {"verdicts": {}, "f1_end_with": {}, "seed_on_sheet": 0, "seed_snap_le_tol": 0,
                                   "lost_share_f1_median": None, "ambiguous_share_f1_median": None}
               for cn, sl in settings}
    shares = {k: ([], []) for k in tallies}
    for c in range(P.shape[1]):
        rows = [r for r in range(P.shape[0]) if kept[r, c]]
        u, cum, base, n, img = L.ribbon(ct, normal_at, P[rows, c], DU, 0.0, S)
        jr = [int(np.argmin(np.abs(u - x))) for x in cum]
        cr = {w: [] for w in L.MESHES}
        for j in jr:
            for w in L.MESHES:
                t = meshes[w].cross(base[j], n[j], tmax)
                cr[w].append(round(t[0] * L.UM, 1) if t else None)
        col = {"rows": rows, "u_um_at_rows": [round(float(u[j] * L.UM), 1) for j in jr], "crossings_at_rows": cr,
               "n_samples": int(len(u))}
        # how dark is the CT along each surface's own path through this column (every sample, seed to last row):
        # level = (I - p10 of the line) / (p90 - p10 of the line); "dark" < 0.25 is gap level, "half" < 0.5
        dark = {}
        p10, p90 = np.percentile(img, 10, axis=1), np.percentile(img, 90, axis=1)
        for tag in ("trace",) + L.MESHES:
            lev = []
            for j in range(len(u)):
                if tag == "trace":
                    x = 0.0
                else:
                    t = meshes[tag].cross(base[j], n[j], tmax)
                    x = t[0] * L.UM if t else None
                if x is None:
                    continue
                v = img[j][int(np.argmin(np.abs(S - x)))]
                lev.append((v - p10[j]) / max(p90[j] - p10[j], 1e-6))
            lev = np.array(lev)
            if lev.size == 0:
                dark[tag] = None
                continue
            run, best = 0, 0
            for v in lev:
                run = run + 1 if v < 0.25 else 0
                best = max(best, run)
            dark[tag] = {"n": int(lev.size), "share_below_0.25": round(float(np.mean(lev < 0.25)), 3),
                         "share_below_0.5": round(float(np.mean(lev < 0.5)), 3),
                         "longest_dark_run_um": round(best * DU * L.UM, 1)}
        col["path_darkness"] = dark
        floor0 = float(np.percentile(img[0], 10))
        i0 = int(np.argmin(np.abs(S)))
        for cn, sl in settings:
            cut = cuts[cn]
            f1, snap1, a1, l1 = f1_ = follow(img, cut, 0, 0.0, +1, sl)
            seed_ridge = f1[0][1]
            k0 = int(np.argmin(np.abs(S - seed_ridge)))
            seed_level = (img[0][i0] - floor0) / max(img[0][k0] - floor0, 1e-6)
            jl = jr[-1]
            f2, snap2, a2, l2 = follow(img, cut, jl, 0.0, -1, sl)
            s35 = cr["w035"][-1]
            f3 = snap3 = a3 = l3 = None
            if s35 is not None:
                f3, snap3, a3, l3 = follow(img, cut, jl, s35, -1, sl)
            pos1 = {j: x for j, x, _ in f1}
            per_row = []
            for k, (r, j) in enumerate(zip(rows, jr)):
                x = pos1.get(j)
                if x is None:
                    per_row.append({"row": r, "f1_s_um": None, "with": []})
                    continue
                near = ["trace"] if abs(x) <= TOL else []
                near += [w for w in L.MESHES if cr[w][k] is not None and abs(x - cr[w][k]) <= TOL]
                per_row.append({"row": r, "f1_s_um": round(x, 1), "with": near})
            reach2 = f2[-1][0] == 0 and abs(f2[-1][1] - seed_ridge) <= TOL
            reach3 = None if f3 is None else (f3[-1][0] == 0 and abs(f3[-1][1] - seed_ridge) <= TOL)
            pl = img[jl]
            fl = float(np.percentile(pl, 10))

            def level(x, r_):
                return (pl[int(np.argmin(np.abs(S - x)))] - fl) / max(pl[int(np.argmin(np.abs(S - r_)))] - fl, 1e-6)

            t_on = bool(level(0.0, f2[0][1]) >= L.HALF)
            w_on = None if f3 is None else bool(level(s35, f3[0][1]) >= L.HALF)
            yn = {True: "reaches the seed's ridge", False: "does not reach the seed's ridge"}
            if f3 is None:
                verdict = "no w035 crossing at the last row"
            elif not t_on and not w_on:
                verdict = "neither trace nor w035 on a sheet at the end"
            elif t_on and not w_on:
                verdict = f"w035 off any sheet at the end; the trace's ridge {yn[reach2]}"
            elif w_on and not t_on:
                verdict = f"trace off any sheet at the end; w035's ridge {yn[reach3]}"
            elif abs(f2[0][1] - f3[0][1]) <= 0.5:
                verdict = f"trace and w035 on one ridge at the end; it {yn[reach2]}"
            elif reach2 and not reach3:
                verdict = "two ridges: the trace's reaches the seed's ridge, w035's does not"
            elif reach3 and not reach2:
                verdict = "two ridges: w035's reaches the seed's ridge, the trace's does not"
            elif reach2 and reach3:
                verdict = "two ridges, both reach the seed's ridge (they touch on the way): CT does not decide"
            else:
                verdict = "two ridges, neither reaches the seed's ridge: undecided"
            last = per_row[-1]
            f1_end = ("left the window" if last["f1_s_um"] is None else
                      "+".join(last["with"]) if last["with"] else "neither trace nor a mesh")
            key = f"{cn}_slope{sl}"
            entry = {"seed_ridge_s_um": round(seed_ridge, 1), "seed_level": round(float(seed_level), 2),
                     "seed_on_sheet": bool(seed_level >= L.HALF),
                     "f1": {"per_row": per_row, "n": len(f1), "ambiguous": a1, "lost": l1, "end_with": f1_end,
                            "min_neighbour_gap_um": round(min((g for _, _, g in f1 if g is not None),
                                                              default=float("nan")), 1)},
                     "last_row": {"trace_on_sheet": t_on, "w035_on_sheet": w_on, "s_w035_um": s35,
                                  "s_w034_um": cr["w034"][-1]},
                     "f2_from_trace": {"start_snap_um": snap2, "end_s_um": round(f2[-1][1], 1),
                                       "ended_at_seed_sample": f2[-1][0] == 0,
                                       "reaches_seed_ridge": bool(reach2), "ambiguous": a2, "lost": l2},
                     "f3_from_w035": None if f3 is None else {
                         "start_s_um": s35, "start_snap_um": snap3, "end_s_um": round(f3[-1][1], 1),
                         "ended_at_seed_sample": f3[-1][0] == 0,
                         "reaches_seed_ridge": bool(reach3), "ambiguous": a3, "lost": l3},
                     "verdict": verdict}
            col[key] = entry
            t = tallies[key]
            t["verdicts"][verdict] = t["verdicts"].get(verdict, 0) + 1
            t["f1_end_with"][f1_end] = t["f1_end_with"].get(f1_end, 0) + 1
            t["seed_on_sheet"] += int(entry["seed_on_sheet"])
            t["seed_snap_le_tol"] += int(abs(seed_ridge) <= TOL)
            shares[key][0].append(l1 / len(f1))
            shares[key][1].append(a1 / len(f1))
            if (cn, sl) == main_set:
                col["path"] = [[round(float(u[j] * L.UM), 1), x] for j, x, _ in f1]
                col["path_f2"] = [[round(float(u[j] * L.UM), 1), x] for j, x, _ in f2]
                col["path_f3"] = None if f3 is None else [[round(float(u[j] * L.UM), 1), x] for j, x, _ in f3]
        res["columns"][str(c)] = col
        m = col[f"{main_set[0]}_slope{main_set[1]}"]
        f3 = m["f3_from_w035"]
        print(f"col {c:2d}: seed ridge {m['seed_ridge_s_um']:6.1f} lvl {m['seed_level']:5.2f} | F1 end "
              f"{m['f1']['end_with']:24s} amb {m['f1']['ambiguous']:2d} lost {m['f1']['lost']:2d}/{m['f1']['n']} | "
              f"F2 snap {m['f2_from_trace']['start_snap_um']:6.1f} end {m['f2_from_trace']['end_s_um']:6.1f} | F3 "
              f"{'-' if f3 is None else (f3['start_s_um'], f3['start_snap_um'], f3['end_s_um'])} | {m['verdict']}")
    for k, t in tallies.items():
        t["lost_share_f1_median"] = round(float(np.median(shares[k][0])), 3)
        t["ambiguous_share_f1_median"] = round(float(np.median(shares[k][1])), 3)
        print(k, json.dumps(t))
    pd_ = {}
    for tag in ("trace",) + L.MESHES:
        for grp, cols in (("all", range(P.shape[1])), ("cols_11_21", range(11, 22)), ("other_cols",
                                                                                       [c for c in range(P.shape[1]) if not 11 <= c <= 21])):
            v = [res["columns"][str(c)]["path_darkness"][tag] for c in cols]
            v = [x for x in v if x is not None]
            if not v:
                continue
            pd_[f"{tag}_{grp}"] = {"columns": len(v),
                                  "share_below_0.25_median": round(float(np.median([x["share_below_0.25"] for x in v])), 3),
                                  "share_below_0.5_median": round(float(np.median([x["share_below_0.5"] for x in v])), 3),
                                  "longest_dark_run_um_median": round(float(np.median([x["longest_dark_run_um"] for x in v])), 1),
                                  "columns_with_dark_run_ge_27um": int(sum(x["longest_dark_run_um"] >= 27 for x in v))}
    res["path_darkness_summary"] = pd_
    for k, v in pd_.items():
        print(k, v)
    res["summary"] = tallies
    json.dump(res, open(os.path.join(L.OUT, "ridge_follow.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
