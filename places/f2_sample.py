# Step 2: per place, fetch exactly the level-0 chunk planes needed (byte ranges) and resample:
#   cut: radial cut (dz, r), dz +-200 um, mean over +-3 um across the line, 1 um steps (as in the pilot)
#   slices at dz -150 / 0 / +150 um: (tang, r) images, mean of 5 planes (+-2.3 um), r = [w033 - 300, w037 + 300], tang +-400
#   mesh traces in each slice and mesh crossings on the ray every 10 um in dz (published 2.399 um copies)
# usage: python -I f2_sample.py E1 [--dry]
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fcommon import *

def plane_trace(g, zc):
    a, b, c, e = g[:-1, :-1], g[1:, :-1], g[:-1, 1:], g[1:, 1:]
    segs = []
    for tri in ((a, b, c), (e, c, b)):
        P = np.stack(tri, -2)
        ok = np.all(np.isfinite(P[..., 0]), -1)
        P = P[ok]
        s = P[..., 2] - zc
        cr = (s.min(1) < 0) & (s.max(1) >= 0)
        P, s = P[cr], s[cr]
        pts = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            m = (s[:, i] < 0) != (s[:, j] < 0)
            f = np.where(m, s[:, i] / np.where(m, s[:, i] - s[:, j], 1), np.nan)
            pts.append(P[:, i] + f[:, None] * (P[:, j] - P[:, i]))
        pts = np.stack(pts, 1)
        out = np.full((len(P), 2, 3), np.nan)
        for n in range(len(P)):
            fin = pts[n][np.isfinite(pts[n][:, 0])]
            if len(fin) >= 2: out[n] = fin[:2]
        segs.append(out)
    return np.concatenate(segs)

name = sys.argv[1]; dry = "--dry" in sys.argv
geo = json.load(open(os.path.join(OUT5, "geo.json")))[name]
z0, a, o, d, t = place_frame(name)
h = geo["hits9"]
lo_r = min(min(h["w033"]), min(h["w034"])); hi_r = max(max(h["w037"]), max(h["w036"]))
r0, r1 = np.floor(lo_r - 200), np.ceil(hi_r + 200)        # cut range (pilot: ~150-250 um beyond w033/w037)
s0, s1 = np.floor(lo_r - 300), np.ceil(hi_r + 300)        # slice range
TAV = 3
rc = np.arange(r0, r1 + 1, 1.0); zz = np.arange(-200, 201, 1.0)
Rc, Zc = np.meshgrid(rc, zz)
cut_pts = [to11(p24(o, d, t, Rc, np.full_like(Rc, k * UM11), Zc).reshape(-1, 3)).reshape(Rc.shape + (3,)) for k in range(-TAV, TAV + 1)]
rr = np.arange(s0, s1 + 1, 1.0); tt = np.arange(-400, 401, 1.0)
R, Tg = np.meshgrid(rr, tt)
DZS = (-150.0, 0.0, 150.0)
sl_pts = {dz: [to11(p24(o, d, t, R, Tg, np.full_like(R, dz + k)).reshape(-1, 3)).reshape(R.shape + (3,)) for k in np.arange(-2, 3) * UM11] for dz in DZS}

need_cut = merge_needs(*[needed(p[::2, ::2]) for p in cut_pts[::2]] + [needed(cut_pts[TAV])])
need_cut = merge_needs(need_cut, needed(np.concatenate([p.reshape(-1, 3) for p in cut_pts])))
need_sl = {dz: merge_needs(*[needed(p) for p in v]) for dz, v in sl_pts.items()}
allneed = merge_needs(need_cut, *need_sl.values())
print(f"{name}: angle {a:.3f}, z {z0:.0f}; cut r {r0:.0f}..{r1:.0f}, slice r {s0:.0f}..{s1:.0f}")
print(f"  cut: {len(need_cut)} chunks {plan_bytes(need_cut)/1e6:.1f} MB; slices: " +
      ", ".join(f"dz {dz:+.0f} {len(n)} chunks {plan_bytes(n)/1e6:.1f} MB" for dz, n in need_sl.items()) +
      f"; union {plan_bytes(allneed)/1e6:.1f} MB")
if dry:
    raise SystemExit
bad = fetch_need(allneed, name)
if bad:
    print("PROBLEMS", bad)

vol, org = assemble(need_cut)
cut = np.nanmean([sample(vol, org, p) for p in cut_pts], 0)
cut1 = sample(vol, org, cut_pts[TAV])
del vol
slices = {}
for dz in DZS:
    vol, org = assemble(need_sl[dz])
    slices[dz] = np.nanmean([sample(vol, org, p) for p in sl_pts[dz]], 0)
    del vol
cen = o + ((r0 + r1) / 2 / UM24) * d
hw = np.array([900, 900, 300]) / UM24
MESH = {w: mesh24_fine(w, (cen - hw, cen + hw)) for w in MESHES}
traces, cuthits = {}, {}
for w, g in MESH.items():
    if g.size == 0:
        cuthits[w] = np.zeros((0, 2)); continue
    for dz in DZS:
        segs = plane_trace(g, z0 + dz / UM24)
        rel = segs - o
        rs, ts = (rel @ d) * UM24, (rel @ t) * UM24
        keep = np.all((rs >= s0) & (rs <= s1) & (np.abs(ts) <= 400), 1)
        traces[f"{w}_{int(dz)}"] = np.stack([rs[keep], ts[keep]], -1)
    hh = []
    for dz in np.arange(-200, 201, 10.0):
        oo = o + np.array([0, 0, dz / UM24])
        for tt_, _, _ in ray_hits(g, oo, d):
            if r0 - 50 <= tt_ * UM24 <= r1 + 50:
                hh.append((tt_ * UM24, dz))
    cuthits[w] = np.array(hh).reshape(-1, 2)
hf = {w: [round(x[0] * UM24, 1) for x in ray_hits(MESH[w], o, d) if r0 - 50 <= x[0] * UM24 <= r1 + 50] if MESH[w].size else [] for w in MESHES}
print("  crossings at dz 0 (2.399 copies):", hf)
np.savez_compressed(os.path.join(OUT5, f"{name}.npz"), cut=cut, cut1=cut1, rc=rc, zz=zz, rr=rr, tt=tt,
                    **{f"slice_{int(dz)}": v for dz, v in slices.items()},
                    **{f"tr_{k}": v for k, v in traces.items()}, **{f"ch_{w}": v for w, v in cuthits.items()},
                    meta=json.dumps({"name": name, "z0": z0, "angle": a, "o": o.tolist(), "hits24": hf}))
print(f"  saved; cut nan {np.isnan(cut).mean():.4f}; slice nan " + ", ".join(f"{np.isnan(v).mean():.4f}" for v in slices.values())
      + f"; crossings per mesh in cut { {w: len(v) for w, v in cuthits.items()} }; dir {dl_bytes()/1e6:.1f} MB")
