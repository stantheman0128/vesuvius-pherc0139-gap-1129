"""Descriptive, outside FINER-RULE.md: half-height width of the single seed ridge (row 0, columns 11 to 17) at
levels 1, 2 and 3, and the CT level at the w034 and w035 crossings, on the same profiles as cal2_look_finer.py.
Width: from the ridge crest outward until the profile drops below p10 + 0.5 x (crest - p10), p10 of that profile.
Level at a crossing: (I - p10) / (crest - p10). Edge: distance from the 90% to the 10% level on each flank, a
measure of how sharp the data are along the normal at each level. Writes out/cal2_look/finer_width.json.
Usage (WSL venv): python -I cal2_look_finer_width.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cal2_look as L  # noqa: E402
import cal2_look_finer as F  # noqa: E402
import numpy as np  # noqa: E402


def main():
    p, keys3 = L.load_plan()
    meshes, normal_at = L.load_meshes(p)
    P, _, _ = L.load_trace("real")
    est = json.load(open(os.path.join(L.OUT, "finer_ct_estimate.json")))
    fl = json.load(open(os.path.join(L.OUT, "finer_look.json")))
    boxes = {lv: F.Box(lv, [tuple(k) for k in est["levels"][str(lv)]["keys_zyx"]]) for lv in (1, 2)}
    boxes[3] = F.Box3(keys3)
    out = {"what": "descriptive, outside FINER-RULE.md", "columns": {}}
    for c in F.COLS:
        b = P[0, c]
        n = normal_at(b[None])[0]
        t1, t2 = F.frame(n, P[0, c + 1] - P[0, c - 1])
        ent = {}
        for lv in (1, 2, 3):
            pr, _ = F.profile(boxes[lv], b, n, t1, t2)
            col = fl["levels"][str(lv)]["columns"][str(c)]["row0"]
            s34, s35 = col["s34"], col["s35"]
            rs, _ = F.ridges(pr, fl["levels"][str(lv)]["cut"])
            w = rs[(rs >= s34 - F.W_PAD) & (rs <= s35 + F.W_PAD)]
            if w.size != 1:
                ent[f"L{lv}"] = None
                continue
            k0 = int(np.argmin(np.abs(F.S - w[0])))
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

            def rise(step):          # 10% to 90% of (crest - p10) going outward from the crest; None if not reached
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

            ent[f"L{lv}"] = {"crest_s_um": float(w[0]), "half_width_um": float(F.S[z] - F.S[a]),
                             "edge_10_90_um_inner_outer": [rise(-1), rise(+1)],
                             "half_lo_hi_um": [float(F.S[a]), float(F.S[z])],
                             "level_at_w034": lvl(s34), "level_at_w035": lvl(s35)}
        out["columns"][str(c)] = ent
        print(c, {k: (v and (v["crest_s_um"], v["half_width_um"], v["half_lo_hi_um"], v["level_at_w034"],
                              v["level_at_w035"])) for k, v in ent.items()})
    for lv in (1, 2, 3):
        ws = [e[f"L{lv}"]["half_width_um"] for e in out["columns"].values() if e[f"L{lv}"]]
        out[f"median_half_width_L{lv}_um"] = float(np.median(ws))
        eg = [x for e in out["columns"].values() if e[f"L{lv}"] for x in e[f"L{lv}"]["edge_10_90_um_inner_outer"]
              if x is not None]
        out[f"median_edge_10_90_L{lv}_um"] = float(np.median(eg)) if eg else None
        print(f"L{lv} median half width {np.median(ws):.1f} um over {len(ws)} columns; median 10-90 edge "
              f"{out[f'median_edge_10_90_L{lv}_um']} um over {len(eg)} flanks")
    json.dump(out, open(os.path.join(L.OUT, "finer_width.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
