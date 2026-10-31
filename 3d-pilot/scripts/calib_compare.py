"""Calibration-place results: steps reached and stop reasons, per-step control-point counts by nearest mesh, the real
trace's distance to w035, and, for the +1 / -1 voxel CT shifts, the per-step mean signed displacement along w035's
outward normal next to the mean t . n (same method as ablate_compare.py: point by point at the same (step, column)
where both runs kept the point, seed row excluded, normal at the nearest w035 band vertex).
Also checks every tracer read (run_sheet_calib.py) against the downloaded chunks and, for shifted views, the
shift-fill voxels.
Usage: python -I calib_compare.py      Writes out/calib/compare.json
"""
import json
import os
import site
import sys

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pick_places as PP  # noqa: E402
from calib_prep import outward_normals, umb11_xy, CH  # noqa: E402

ROOT = PP.ROOT
C = os.path.join(ROOT, "out", "calib")
UM = 1.129
SHAPE_XYZ = np.array([4040, 4497, 2413])   # level-3 zarr shape (z, y, x) = (2413, 4497, 4040)


def summ(v):
    v = np.asarray(v, np.float64)
    if v.size == 0:
        return {"n": 0}
    return {"n": int(v.size), "mean": round(float(v.mean()), 2), "median": round(float(np.median(v)), 2),
            "p10": round(float(np.percentile(v, 10)), 2), "p90": round(float(np.percentile(v, 90)), 2)}


def load_trace(d):
    z = np.load(os.path.join(d, "trace.npz"))
    return np.stack([z[c] for c in "xyz"], -1).astype(np.float64) * 8


plans = {p["name"]: p for p in json.load(open(os.path.join(C, "plans.json")))}
views = {v["name"]: v for v in json.load(open(os.path.join(C, "views.json")))}
res = {}
for name, p in plans.items():
    Zh = p["height"]
    pd = os.path.join(C, name)
    r = {"height": Zh, "kind": p["kind"], "p10_nearer_um": p["p10_nearer_um"], "angle_range": p["angle_range"],
         "direction": p["direction"], "axis": p["axis_xyz"], "angle_normal_to_axis_deg": p["angle_normal_to_axis_deg"],
         "runs": {}}
    # meshes near the place
    lo_a, hi_a = p["angle_range"][0] - 15, p["angle_range"][1] + 15
    trees = {}
    for w in ("w034", "w035", "w036"):
        g = np.load(os.path.join(ROOT, "downloads", "mesh_rows", f"{w}_z{Zh}.npz"))["g"].astype(np.float64)
        v = np.all(g != -1, -1)
        a, _ = PP.angle24(g)
        sel = v & (a >= lo_a) & (a <= hi_a)
        rr, cc = np.nonzero(sel.any(1))[0], np.nonzero(sel.any(0))[0]
        sub = g[rr.min():rr.max() + 1, cc.min():cc.max() + 1]
        sv = np.all(sub != -1, -1)
        trees[w] = cKDTree(PP.densify(sub, sv))
        if w == "w035":
            n = outward_normals(sub, umb11_xy)
            band_tree, band_n = cKDTree(sub[sv]), n[sv]
    # chunk-presence mask over the box
    keys = np.array(p["chunks"])
    klo = keys.min(0)
    shape = tuple((keys.max(0) - klo + 1) * CH)
    have = np.zeros(shape, bool)
    for k in keys:
        o = (k - klo) * CH
        have[o[0]:o[0] + CH, o[1]:o[1] + CH, o[2]:o[2] + CH] = True
    tplus = np.array(views[name]["t_plus_um_xyz"])
    traces = {}
    for v in ("real", "p", "m"):
        d = os.path.join(pd, f"trace_{v}")
        if not os.path.exists(os.path.join(d, "trace.npz")):
            r["runs"][v] = "not run"
            continue
        P = load_trace(d)
        traces[v] = P
        rec = json.load(open(d + "_record.json"))
        raw = np.load(d + "_raw.npz")
        reads = json.load(open(d + "_reads.json"))
        steps = P.shape[0] - 1
        last = raw["rows"][-1]
        if steps >= 8:
            stop = "completed 8 steps"
        elif np.any(last < 0) or np.any(last >= SHAPE_XYZ):
            stop = "a point left the volume (inference.py:152)"
        else:
            stop = "next row had no confident point connected to the trace (inference.py:148)"
        # reads vs downloaded chunks and fill
        if v == "real":
            ok_mask = have
        else:
            ok_mask = np.unpackbits(np.load(os.path.join(pd, f"valid_{v}.npy")))[:int(np.prod(shape))].reshape(shape)
            ok_mask = ok_mask.astype(bool)
        bad = 0
        for rd in reads:
            o = np.array(rd["origin_zyx"]) - klo * CH
            s = np.array(rd["shape_zyx"])
            if np.any(o < 0) or np.any(o + s > np.array(shape)) or not ok_mask[o[0]:o[0] + s[0], o[1]:o[1] + s[1],
                                                                             o[2]:o[2] + s[2]].all():
                bad += 1
        kept = np.all(P >= 0, -1)
        per_step = []
        for i in range(1, P.shape[0]):
            cp = [c for c in range(0, P.shape[1], 6) if kept[i, c]]
            if cp:
                dd = np.stack([trees[w].query(P[i, cp])[0] for w in ("w034", "w035", "w036")])
                nearest = np.argmin(dd, 0)
                counts = [int((nearest == j).sum()) for j in range(3)]
            else:
                counts = [0, 0, 0]
            per_step.append({"step": i, "control_kept": len(cp), "nearest_w034_w035_w036": counts})
        q = P[1:][kept[1:]]
        run = {"steps_reached": int(steps), "stop": stop, "wall_s": rec["wall_s"], "status": rec["status"],
               "chunks_missing": rec["ct_chunks_missing"], "reads": len(reads), "reads_touching_fill_or_absent": bad,
               "per_step_control": per_step,
               "control_on_w035": [sum(s["nearest_w034_w035_w036"][1] for s in per_step),
                                   sum(s["control_kept"] for s in per_step)]}
        if len(q):
            dist = {w: trees[w].query(q)[0] * UM for w in ("w034", "w035", "w036")}
            run["dense_to_w035_um"] = summ(dist["w035"])
            run["dense_nearest_counts"] = [int(x) for x in np.bincount(
                np.argmin(np.stack([dist[w] for w in ("w034", "w035", "w036")]), 0), minlength=3)]
            m = dist["w035"] / np.minimum(dist["w034"], dist["w036"])
            run["margin_m"] = {"median": round(float(np.median(m)), 3), "p90": round(float(np.percentile(m, 90)), 3),
                               "max": round(float(m.max()), 3)}
        r["runs"][v] = run
    if "real" in traces:
        gr_all = traces["real"]
        for v, sgn in (("p", 1), ("m", -1)):
            if v not in traces:
                continue
            ga_all = traces[v]
            rows = min(ga_all.shape[0], gr_all.shape[0])
            ga, gr = ga_all[:rows], gr_all[:rows]
            ok = np.all(ga >= 0, -1) & np.all(gr >= 0, -1)
            ok[0] = False
            t = sgn * tplus
            st, pooled_s, pooled_t = [], [], []
            for i in range(1, rows):
                if not ok[i].any():
                    st.append({"step": i, "n": 0})
                    continue
                disp = (ga[i][ok[i]] - gr[i][ok[i]]) * UM
                nn = band_n[band_tree.query(gr[i][ok[i]])[1]]
                s_ = np.sum(disp * nn, -1)
                tn = nn @ t
                st.append({"step": i, "n": int(ok[i].sum()), "mean_signed_um": round(float(s_.mean()), 2),
                           "mean_t_dot_n_um": round(float(tn.mean()), 2)})
                if i >= 5:
                    pooled_s.append(s_)
                    pooled_t.append(tn)
            late = {}
            if pooled_s:
                s5, t5 = np.concatenate(pooled_s), np.concatenate(pooled_t)
                late = {"n": int(s5.size), "mean_signed_um": round(float(s5.mean()), 2),
                        "mean_t_dot_n_um": round(float(t5.mean()), 2),
                        "ratio": round(float(s5.mean() / t5.mean()), 3)}
            r["runs"][v]["vs_real"] = {"t_um_xyz": t.round(3).tolist(), "per_step": st, "steps_5_plus": late}
    res[name] = r
json.dump(res, open(os.path.join(C, "compare.json"), "w"), indent=1)
for name, r in res.items():
    print(name, r["kind"], "p10", r["p10_nearer_um"], "axis", r["axis"], r["angle_normal_to_axis_deg"])
    for v, run in r["runs"].items():
        if isinstance(run, str):
            print("  ", v, run)
            continue
        print(f"   {v}: steps {run['steps_reached']} stop: {run['stop']} | on w035 {run['control_on_w035']} | "
              f"to w035 {run.get('dense_to_w035_um')} | reads bad {run['reads_touching_fill_or_absent']} "
              f"missing {run['chunks_missing']}")
        print("      control per step", [(s["control_kept"], s["nearest_w034_w035_w036"]) for s in run["per_step_control"]])
        if "vs_real" in run:
            print("      per step signed vs t.n", [(s["step"], s.get("mean_signed_um"), s.get("mean_t_dot_n_um"))
                                                 for s in run["vs_real"]["per_step"]], "late", run["vs_real"]["steps_5_plus"])
