"""Checks next to cal2_look_level0.py, outside FINER-RULE.md. (1) Pyramid: level-1 chunk (27, 49, 59) against the
2x2x2 block mean and the decimation of its 8 level-0 children, which tests the level-0 key, layout and axis order.
(2) Level-0 fine-scale texture: neighbour differences and lag correlations of a high-pass residual in one chunk.
(3) The two sub-cut maxima of level0_look.json (columns 13 and 17, row 0): where the dip is and how deep.
Writes out/cal2_look/level0_checks.json.
Usage (WSL venv): python -I cal2_look_level0_checks.py
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import cal2_look as L  # noqa: E402
import cal2_look_finer as F  # noqa: E402
import numpy as np  # noqa: E402

CH = 128
L0D = os.path.join(ROOT, "downloads", "ct_level0", "0")
L1D = os.path.join(ROOT, "downloads", "ct", "vol.zarr", "1")


def chunk(d, k):
    return np.fromfile(os.path.join(d, *map(str, k)), np.uint8).reshape(CH, CH, CH)


# 1. pyramid: level-1 chunk (27, 49, 59) against the 2x2x2 block mean of its 8 level-0 children
k1 = (27, 49, 59)
a1 = chunk(L1D, k1).astype(np.float64)
big = np.zeros((256, 256, 256), np.float64)
for dz in (0, 1):
    for dy in (0, 1):
        for dx in (0, 1):
            big[dz * CH:(dz + 1) * CH, dy * CH:(dy + 1) * CH, dx * CH:(dx + 1) * CH] = chunk(
                L0D, (2 * k1[0] + dz, 2 * k1[1] + dy, 2 * k1[2] + dx))
m = (a1 > 0)
res = {}
mean2 = big.reshape(128, 2, 128, 2, 128, 2).mean(axis=(1, 3, 5))
dec = big[::2, ::2, ::2]
for name, b in (("block_mean_2x2x2", mean2), ("decimate_even", dec), ("decimate_odd", big[1::2, 1::2, 1::2])):
    mm = m & (b > 0)
    d = a1[mm] - b[mm]
    res[name] = {"corr": round(float(np.corrcoef(a1[mm], b[mm])[0, 1]), 5), "mean_abs_diff": round(float(np.abs(d).mean()), 3),
                 "max_abs_diff": round(float(np.abs(d).max()), 2), "n": int(mm.sum())}
# shift test on the block mean: does a one-level-0-voxel shift of the block grid fit better?
for sh in (-1, 1):
    b2 = np.roll(big, sh, axis=(0, 1, 2)).reshape(128, 2, 128, 2, 128, 2).mean(axis=(1, 3, 5))
    sl = (slice(2, -2),) * 3
    mm = (m & (b2 > 0))[sl]
    d = a1[sl][mm] - b2[sl][mm]
    res[f"block_mean_rolled_{sh}"] = {"mean_abs_diff": round(float(np.abs(d).mean()), 3)}
OUTJ = {"pyramid_L1_27_49_59": res}
print("pyramid", json.dumps(res))

# fine-scale texture at level 0: fraction of zero neighbour differences along x, and lag correlations of a
# high-pass residual (voxel minus 3x3x3 mean)
c0 = chunk(L0D, (54, 98, 119)).astype(np.float64)
nz = c0 > 0
dx = np.diff(c0, axis=2)[(nz[:, :, 1:] & nz[:, :, :-1])]
OUTJ["L0_texture_chunk_54_98_119"] = {"x_diff_frac_zero": round(float((dx == 0).mean()), 3),
    "x_diff_frac_abs_le_1": round(float((np.abs(dx) <= 1).mean()), 3), "x_diff_std": round(float(dx.std()), 3),
    "highpass": "voxel minus 5x5x5 mean", "lag_corr_x": {}}
print("L0 x-diff: frac zero", round(float((dx == 0).mean()), 3), "frac |d|<=1", round(float((np.abs(dx) <= 1).mean()), 3),
      "std", round(float(dx.std()), 3))
from scipy.ndimage import uniform_filter  # noqa: E402
hp = c0 - uniform_filter(c0, 5)
inner = (slice(8, -8),) * 3
h = hp[inner]
ok = nz[inner]
for lag in (1, 2, 3, 4, 6, 8):
    a = h[:, :, :-lag][ok[:, :, :-lag] & ok[:, :, lag:]]
    b = h[:, :, lag:][ok[:, :, :-lag] & ok[:, :, lag:]]
    OUTJ["L0_texture_chunk_54_98_119"]["lag_corr_x"][str(lag)] = round(float(np.corrcoef(a, b)[0, 1]), 3)
    print(f"L0 high-pass lag {lag} vox ({lag * 1.129:.1f} um) corr along x: {np.corrcoef(a, b)[0, 1]:.3f}")

# 2. the sub-cut maxima at columns 13 and 17, row 0: where is the dip, between which maxima, how deep
p, keys3 = L.load_plan()
meshes, normal_at = L.load_meshes(p)
P, _, _ = L.load_trace("real")
est0 = json.load(open(os.path.join(L.OUT, "level0_ct_estimate.json")))
est = json.load(open(os.path.join(L.OUT, "finer_ct_estimate.json")))
F.DST = os.path.join(ROOT, "downloads", "ct_level0")
b0 = F.Box(0, [tuple(k) for k in est0["level"]["keys_zyx"]])
F.DST = os.path.join(ROOT, "downloads", "ct", "vol.zarr")
b1 = F.Box(1, [tuple(k) for k in est["levels"]["1"]["keys_zyx"]])
lk = json.load(open(os.path.join(L.OUT, "level0_look.json")))
OUTJ["subcut_dips_row0"] = {}
for c in (13, 17):
    b = P[0, c]
    n = normal_at(b[None])[0]
    t1, t2 = F.frame(n, P[0, c + 1] - P[0, c - 1])
    row = lk["levels"]["0"]["columns"][str(c)]["row0"]
    s34, s35 = row["s34"], row["s35"]
    for lv, box in ((0, b0), (1, b1)):
        pr, _ = F.profile(box, b, n, t1, t2)
        dc = lk["descriptive"]["columns"][str(c)][f"L{lv}"]["dip_check"]
        mx = sorted(x for x, _ in dc["maxima_in_W0"])
        i_a, i_b = (int(np.argmin(np.abs(F.S - x))) for x in mx)
        j = i_a + int(np.argmin(pr[i_a:i_b + 1]))
        f10 = float(np.percentile(pr, 10))
        crest = float(pr[int(np.argmin(np.abs(F.S - dc["main_s_um"])))])
        seg = pr[i_a:i_b + 1]
        OUTJ["subcut_dips_row0"][f"{c}_L{lv}"] = {"s34": s34, "s35": s35, "maxima_in_W0": dc["maxima_in_W0"],
            "dip_min_s_um": float(F.S[j]), "dip_min_value": round(float(pr[j]), 2),
            "values_at_maxima": [round(float(pr[i]), 2) for i in (i_a, i_b)], "crest_minus_p10": round(crest - f10, 1),
            "dip_depth_over_crest_minus_p10": round(dc["largest_other_prominence_in_W0"] / (crest - f10), 3),
            "dip_between_crossings": bool(min(s34, s35) <= F.S[j] <= max(s34, s35))}
        print(f"col {c} L{lv}: s34 {s34} s35 {s35}; maxima {dc['maxima_in_W0']}; dip min at s={F.S[j]:.0f} um, "
              f"value {pr[j]:.2f}; values at maxima {[round(float(pr[i]), 2) for i in (i_a, i_b)]}; crest-p10 "
              f"{crest - f10:.1f}; dip depth / (crest-p10) {dc['largest_other_prominence_in_W0'] / (crest - f10):.3f}; "
              f"dip between crossings: {min(s34, s35) <= F.S[j] <= max(s34, s35)}")
        print("   profile s=-10..+20:", [round(float(v), 1) for v in pr[(F.S >= -10) & (F.S <= 20)]])

# 4. deepest dip strictly between the w034 and w035 crossings, row 0, every level: at each sample j in the
# interval, depth = min(highest value from the interval start to j, highest value from j to the interval end) - pr[j]
F.DST = os.path.join(ROOT, "downloads", "ct", "vol.zarr")
b2 = F.Box(2, [tuple(k) for k in est["levels"]["2"]["keys_zyx"]])
b3 = F.Box3(keys3)
prof_cut = {lv: lk["levels"][str(lv)]["cut"] for lv in (0, 1, 2, 3)}
OUTJ["deepest_dip_between_crossings_row0"] = {}
for c in F.COLS:
    b = P[0, c]
    n = normal_at(b[None])[0]
    t1, t2 = F.frame(n, P[0, c + 1] - P[0, c - 1])
    row = lk["levels"]["0"]["columns"][str(c)]["row0"]
    lo, hi = sorted((row["s34"], row["s35"]))
    ent = {}
    for lv, box in ((0, b0), (1, b1), (2, b2), (3, b3)):
        pr, _ = F.profile(box, b, n, t1, t2)
        seg = pr[(F.S >= lo) & (F.S <= hi)]
        ss = F.S[(F.S >= lo) & (F.S <= hi)]
        depth = np.minimum(np.maximum.accumulate(seg), np.maximum.accumulate(seg[::-1])[::-1]) - seg
        j = int(np.argmax(depth))
        ent[f"L{lv}"] = {"depth": round(float(depth[j]), 2), "at_s_um": float(ss[j]),
                         "over_cut": round(float(depth[j]) / prof_cut[lv], 3)}
    OUTJ["deepest_dip_between_crossings_row0"][str(c)] = {"s34": row["s34"], "s35": row["s35"], **ent}
    print(f"col {c} deepest dip between crossings [{lo}, {hi}]:", ent)

json.dump(OUTJ, open(os.path.join(L.OUT, "level0_checks.json"), "w"), indent=1)
