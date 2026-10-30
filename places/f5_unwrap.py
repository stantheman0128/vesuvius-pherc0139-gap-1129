# Polar unwrap of one z slice (mean of 5 planes, +-2.3 um) over angles A0..A1 in the 1.129 um scan, with the pilot's
# run counting applied per 1-degree angle row instead of per 10 um dz row (41 rows for 55..95 deg).
# usage: f5_unwrap.py NAME Z0 [--dry]
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fcommon import *
from collections import Counter
from scipy import ndimage as ndi
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

name, z0 = sys.argv[1], float(sys.argv[2]); dry = "--dry" in sys.argv
A0, A1 = 55.0, 95.0
xu, yu = umb(z0); o = np.array([xu, yu, z0])
M9 = {w: mesh24(w) for w in MESHES}
def dvec(a): return np.array([np.cos(np.radians(a)), np.sin(np.radians(a)), 0.0])
hits9 = {a: {w: [h[0] * UM24 for h in ray_hits(g, o, dvec(a))] for w, g in M9.items()} for a in np.arange(A0, A1 + 0.1, 5)}
lo = min(min(v["w033"] + v["w034"]) for v in hits9.values() if v["w033"] + v["w034"])
hi = max(max(v["w036"] + v["w037"]) for v in hits9.values() if v["w036"] + v["w037"])
r0, r1 = np.floor(lo - 120), np.ceil(hi + 120)
rc = np.arange(r0, r1 + 1, 1.0)
rmid = (r0 + r1) / 2
da = np.degrees(1.0 / rmid)                     # ~1 um of arc at mid radius
ang = np.arange(A0 - 0.5, A1 + 0.5 + 1e-9, da)
AA, RR = np.meshgrid(ang, rc, indexing="ij")
pts = [to11((o + np.stack([RR / UM24 * np.cos(np.radians(AA)), RR / UM24 * np.sin(np.radians(AA)), np.full_like(RR, dz / UM24)], -1)).reshape(-1, 3)).reshape(AA.shape + (3,))
       for dz in np.arange(-2, 3) * UM11]
need = merge_needs(*[needed(p) for p in pts])
print(f"{name}: z {z0:.0f}, angles {A0}-{A1}, r {r0:.0f}..{r1:.0f}; {len(need)} chunks, {plan_bytes(need)/1e6:.1f} MB; hits at 5-deg steps:")
for a, v in hits9.items():
    print(f"   {a:5.1f}: " + "  ".join(f"{w} {[round(x) for x in v[w]]}" for w in MESHES))
if dry: raise SystemExit
bad = fetch_need(need, name)
vol, org = assemble(need)
img = np.nanmean([sample(vol, org, p) for p in pts], 0); del vol
print("nan", np.isnan(img).mean())

# mesh crossings every 0.25 deg (2.399 um copies, cropped)
cen = o + (rmid / UM24) * dvec((A0 + A1) / 2)
half = rmid * np.radians((A1 - A0) / 2 + 2) / UM24 + 400 / UM24
box = (cen - np.array([half, half, 300 / UM24]), cen + np.array([half, half, 300 / UM24]))
MF = {w: mesh24_fine(w, box) for w in MESHES}
cr = {w: [] for w in MESHES}
for a in np.arange(A0, A1 + 1e-9, 0.25):
    for w, g in MF.items():
        if g.size:
            for t, _, _ in ray_hits(g, o, dvec(a)):
                if r0 <= t * UM24 <= r1: cr[w].append((t * UM24, a))
cr = {w: np.array(v).reshape(-1, 2) for w, v in cr.items()}

# run counting per integer degree (rows within +-3 um of arc averaged, as in the cut)
sm = ndi.gaussian_filter(img, 2)
v = sm[np.isfinite(sm)]; h_, e_ = np.histogram(v, 256); c_ = (e_[:-1] + e_[1:]) / 2
w0 = np.cumsum(h_); w1 = w0[-1] - w0; m0 = np.cumsum(h_ * c_)
thr = c_[np.argmax(w0 * w1 * (m0 / np.maximum(w0, 1) - (m0[-1] - m0) / np.maximum(w1, 1)) ** 2)]
mask = ndi.binary_opening(sm > thr, np.ones((3, 3), bool))
lab, n = ndi.label(mask, np.ones((3, 3))); sz = ndi.sum(mask, lab, range(1, n + 1))
mask = np.isin(lab, 1 + np.where(sz >= 200)[0])
def runs(row):
    r = ndi.binary_closing(row, np.ones(4, bool)) | row
    l, _ = ndi.label(r)
    return [(rc[s[0].start], rc[s[0].stop - 1]) for s in ndi.find_objects(l) if s[0].stop - s[0].start >= 6]
def run_of(r, rr):
    for i, (a, b) in enumerate(rr):
        if a - 0.5 <= r <= b + 0.5: return i, 0.0
    ds = [min(abs(r - a), abs(r - b)) for a, b in rr]; i = int(np.argmin(ds)); return i, ds[i]
rows = []
for a in np.arange(A0, A1 + 1e-9, 1.0):
    k = int(np.argmin(np.abs(ang - a)))
    rr = runs(mask[k])
    ent = {"angle": a, "runs": [(round(x), round(y)) for x, y in rr], "mesh": {}}
    for w in MESHES:
        sel = cr[w][np.abs(cr[w][:, 1] - a) < 0.01, 0] if len(cr[w]) else []
        for j, r in enumerate(sorted(sel)):
            key = w if j == 0 else f"{w} end" if j == len(sel) - 1 else f"{w} mid{j}"
            i, d = run_of(r, rr); ent["mesh"][key] = {"r": round(float(r), 1), "run": i, "dist_um": round(float(d), 1)}
    m = ent["mesh"]
    if "w034" in m and "w036" in m:
        ent["d46"] = m["w036"]["run"] - m["w034"]["run"]
        for sh in (-30, 30):
            ent[f"d46_{sh}"] = run_of(m["w036"]["r"] + sh, rr)[0] - run_of(m["w034"]["r"] + sh, rr)[0]
    for w in MESHES:
        if f"{w} end" in m: ent[f"enddiff_{w}"] = m[f"{w} end"]["run"] - m[w]["run"]
    rows.append(ent)
rep = {"name": name, "z0": z0, "otsu": float(thr), "rows": rows,
       "d46_hist": dict(Counter(e.get("d46") for e in rows if "d46" in e)),
       "d46_hist_-30": dict(Counter(e.get("d46_-30") for e in rows if "d46" in e)),
       "d46_hist_30": dict(Counter(e.get("d46_30") for e in rows if "d46" in e)),
       **{f"enddiff_hist_{w}": dict(Counter(e[f"enddiff_{w}"] for e in rows if f"enddiff_{w}" in e)) for w in MESHES}}
json.dump(rep, open(os.path.join(OUT5, f"{name}_layers.json"), "w"), indent=1, default=float)
np.savez_compressed(os.path.join(OUT5, f"{name}.npz"), img=img, ang=ang, rc=rc, **{f"cr_{w}": v for w, v in cr.items()})

COL = {"w033": "#1f77b4", "w034": "#d62728", "w035": "#9467bd", "w036": "#2ca02c", "w037": "#ff7f0e"}
fig, axs = plt.subplots(1, 2, figsize=(22, 16))
lo_, hi_ = np.nanpercentile(img, [1, 99.5])
ext = [rc[0] - .5, rc[-1] + .5, ang[0], ang[-1]]
for i, ax in enumerate(axs):
    ax.imshow(img, cmap="gray", vmin=lo_, vmax=hi_, extent=ext, origin="lower", aspect="auto", interpolation="nearest")
    if i:
        ax.contour(rc, ang, mask.astype(float), [0.5], colors="cyan", linewidths=0.5)
        for e in rows:
            for a_, b_ in e["runs"]: ax.plot([a_, b_], [e["angle"]] * 2, color="yellow", lw=0.8)
    for w in MESHES:
        if len(cr[w]): ax.plot(cr[w][:, 0], cr[w][:, 1], ".", color=COL[w], ms=2.5, label=w)
    ax.set_xlabel("r from umbilicus (um)"); ax.set_ylabel("angle (deg)"); ax.legend(loc="upper left", markerscale=4)
axs[0].set_title(f"{name}: polar unwrap of z {z0:,.0f} in the 1.129 um scan (mean of 5 planes); mesh crossings every 0.25 deg")
axs[1].set_title(f"same, Otsu {thr:.0f} bodies (cyan) and material runs every 1 deg (yellow)")
fig.tight_layout(); fig.savefig(os.path.join(OUT5, f"{name}_unwrap.png"), dpi=90); plt.close(fig)
print(f"Otsu {thr:.1f}; d46 {rep['d46_hist']} (-30 {rep['d46_hist_-30']}, +30 {rep['d46_hist_30']}); ends " +
      str({k: v for k, v in rep.items() if k.startswith('enddiff') and v}))
for e in rows:
    print(f"  {e['angle']:5.1f}: d46 {e.get('d46')} | " + "  ".join(f"{k}:{v['r']:.0f}->{v['run']}" + (f"(+{v['dist_um']:.0f})" if v['dist_um'] else "")
                                                for k, v in sorted(e["mesh"].items(), key=lambda kv: kv[1]["r"])))
