"""Compare the copy model's _front/_back outputs at C240 with the published neighbouring meshes.

Distances are point-to-surface in um (1.129 um voxels), against the on-1.129 tifxyz copies of w034/w035/w036
upsampled 8x per grid cell (about 2.8 um spacing, so a nearest-point distance is at most about 2 um above the true
point-to-surface distance). A naive baseline moves every w035 vertex along its normal by the median w035-to-
neighbour distance; if the model is no closer than that, it is not using the CT.
Usage: python -I analyze_copy.py <run_dir> <tag>
Writes out/<tag>_metrics.json, out/png/<tag>_slice.png, out/png/<tag>_distmaps.png.
"""
import json
import os
import site
import sys

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
import tifffile  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import ROOT  # noqa: E402

UM = 1.129
run_dir, tag = sys.argv[1], sys.argv[2]
OUT = os.path.join(ROOT, "out")
PNG = os.path.join(OUT, "png")
os.makedirs(PNG, exist_ok=True)


def load_grid(d):
    g = np.stack([tifffile.imread(os.path.join(d, f"{c}.tif")).astype(np.float64) for c in "xyz"], -1)
    v = np.all(g != -1, -1) & np.all(np.isfinite(g), -1)
    return g, v


def densify(g, v, k=8):
    """Bilinear points inside every grid cell whose four corners are valid."""
    q = v[:-1, :-1] & v[1:, :-1] & v[:-1, 1:] & v[1:, 1:]
    ii, jj = np.nonzero(q)
    t = (np.arange(k) + 0.5) / k
    a, b = np.meshgrid(t, t, indexing="ij")
    a, b = a.ravel()[None, :, None], b.ravel()[None, :, None]
    p00, p10 = g[ii, jj][:, None], g[ii + 1, jj][:, None]
    p01, p11 = g[ii, jj + 1][:, None], g[ii + 1, jj + 1][:, None]
    pts = (1 - a) * (1 - b) * p00 + a * (1 - b) * p10 + (1 - a) * b * p01 + a * b * p11
    return np.concatenate([pts.reshape(-1, 3), g[v]])


def normals(g, v):
    dr = np.gradient(g, axis=0)
    dc = np.gradient(g, axis=1)
    n = np.cross(dc, dr)
    n /= np.linalg.norm(n, axis=-1, keepdims=True) + 1e-12
    return n


def stats(d):
    d = np.asarray(d) * UM
    return {"n": int(d.size), "median_um": round(float(np.median(d)), 1), "p90_um": round(float(np.percentile(d, 90)), 1),
            "max_um": round(float(d.max()), 1), "mean_um": round(float(d.mean()), 1)}


bands = {w: np.load(os.path.join(ROOT, "downloads", "mesh_bands", f"{w}_band.npz")) for w in ("w034", "w035", "w036")}
trees, dense = {}, {}
for w, b in bands.items():
    g = b["g"].astype(np.float64)
    v = np.all(g != -1, -1) & np.all(np.isfinite(g), -1)
    dense[w] = densify(g, v)
    trees[w] = cKDTree(dense[w])

patch_dir = os.path.join(OUT, "inputs", "w035_C240_patch")
pg, pv = load_grid(patch_dir)
pn = normals(pg, pv)
res = {"tag": tag, "run_dir": run_dir}
P = pg[pv]
res["input_w035_to"] = {w: stats(trees[w].query(P)[0]) for w in ("w034", "w035", "w036")}

outs = {}
for side in ("front", "back"):
    d = [x for x in os.listdir(run_dir) if x.endswith("_" + side)]
    g, v = load_grid(os.path.join(run_dir, d[0]))
    outs[side] = (g, v)
    Q = g[v]
    dd = {w: trees[w].query(Q)[0] for w in ("w034", "w035", "w036")}
    near = np.argmin(np.stack([dd["w034"], dd["w035"], dd["w036"]]), 0)
    res[side] = {"valid": int(v.sum()), "of": int(v.size),
                 "to": {w: stats(dd[w]) for w in dd},
                 "nearest_mesh_counts": {w: int((near == i).sum()) for i, w in enumerate(("w034", "w035", "w036"))}}
    # displacement length from the input vertex (same grid cell)
    m = v & pv
    res[side]["displacement_um"] = stats(np.linalg.norm(g[m] - pg[m], axis=-1))
    # signed offset along the w035 normal
    res[side]["signed_normal_offset_um_median"] = round(float(np.median(np.sum((g[m] - pg[m]) * pn[m], -1)) * UM), 1)

# naive baseline: offset w035 by the median distance to each neighbour along +/- normal (sign chosen per neighbour)
for w in ("w034", "w036"):
    dist = np.median(trees[w].query(P)[0])
    best = None
    for s in (1, -1):
        Qb = P + s * dist * pn[pv]
        dq = trees[w].query(Qb)[0]
        if best is None or np.median(dq) < np.median(best[1]):
            best = (s, dq)
    res[f"naive_offset_to_{w}"] = {"sign": best[0], "offset_um": round(float(dist * UM), 1), "to": stats(best[1])}

json.dump(res, open(os.path.join(OUT, f"{tag}_metrics.json"), "wt"), indent=1)
print(json.dumps(res, indent=1))

# ---- figure 1: CT z-slice (level 2) with mesh intersections ----
plan = json.load(open(os.path.join(OUT, "plan_ct_l2.json")))
(z0, y0, x0), (z1, y1, x1) = plan["boxes"][0]
ct_root = os.path.join(ROOT, "downloads", "ct", "vol.zarr", "2")
zc = int(round(np.median(P[:, 2]) / 4.0))
img = np.zeros((y1 - y0, x1 - x0), np.uint8)
for cy in range(y0 // 128, (y1 - 1) // 128 + 1):
    for cx in range(x0 // 128, (x1 - 1) // 128 + 1):
        fn = os.path.join(ct_root, str(zc // 128), str(cy), str(cx))
        a = np.memmap(fn, dtype=np.uint8, mode="r", shape=(128, 128, 128))[zc % 128]
        ys, xs = max(y0, cy * 128), max(x0, cx * 128)
        ye, xe = min(y1, cy * 128 + 128), min(x1, cx * 128 + 128)
        img[ys - y0:ye - y0, xs - x0:xe - x0] = a[ys - cy * 128:ye - cy * 128, xs - cx * 128:xe - cx * 128]
zl0 = zc * 4.0


def cut(g, v, z):
    """Points where grid edges cross the plane z (level-0 voxels)."""
    pts = []
    for ax in (0, 1):
        a = g[:-1] if ax == 0 else g[:, :-1]
        b = g[1:] if ax == 0 else g[:, 1:]
        va = v[:-1] if ax == 0 else v[:, :-1]
        vb = v[1:] if ax == 0 else v[:, 1:]
        ok = va & vb & ((a[..., 2] - z) * (b[..., 2] - z) <= 0) & (a[..., 2] != b[..., 2])
        t = (z - a[..., 2][ok]) / (b[..., 2][ok] - a[..., 2][ok])
        pts.append(a[ok] + t[:, None] * (b[ok] - a[ok]))
    return np.concatenate(pts) if pts else np.zeros((0, 3))


fig, axs = plt.subplots(1, 2, figsize=(15, 7.2))
ext = [x0 * 4 * UM / 1000, x1 * 4 * UM / 1000, y1 * 4 * UM / 1000, y0 * 4 * UM / 1000]
cols = {"w034": "#1f77b4", "w035": "#2ca02c", "w036": "#d62728", "front": "#ff7f0e", "back": "#9467bd"}
for k, ax in enumerate(axs):
    ax.imshow(img, cmap="gray", extent=ext, vmin=0, vmax=255)
    for w in ("w034", "w035", "w036"):
        b = bands[w]["g"].astype(np.float64)
        vb = np.all(b != -1, -1)
        c = cut(b, vb, zl0)
        ax.plot(c[:, 0] * UM / 1000, c[:, 1] * UM / 1000, ".", ms=3 if k == 0 else 5, color=cols[w],
                label=f"{w} (published)")
    for side, (g, v) in outs.items():
        c = cut(g, v, zl0)
        ax.plot(c[:, 0] * UM / 1000, c[:, 1] * UM / 1000, "x", ms=5 if k == 0 else 8, mew=1.5, color=cols[side],
                label=f"copy model {side}")
    ax.set_xlabel("x (mm, 1.129 um frame)")
    ax.set_ylabel("y (mm)")
axs[0].set_xlim(ext[0], ext[1])
axs[0].set_ylim(ext[2], ext[3])
cx_ = min(max(np.median(P[:, 0]) * UM / 1000, ext[0] + 0.35), ext[1] - 0.35)
cy_ = min(max(np.median(P[:, 1]) * UM / 1000, ext[3] + 0.3), ext[2] - 0.3)
axs[1].set_xlim(cx_ - 0.35, cx_ + 0.35)
axs[1].set_ylim(cy_ + 0.3, cy_ - 0.3)
axs[0].legend(loc="lower left", fontsize=8)
axs[0].set_title(f"C240, z = {zl0:.0f} (1.129 um frame), CT level 2 (4.5 um); input mesh = w035 patch only")
axs[1].set_title("zoom: published w034/w035/w036 (dots) vs copy-model outputs (x)")
fig.tight_layout()
fig.savefig(os.path.join(PNG, f"{tag}_slice.png"), dpi=110)

# ---- figure 2: per-vertex distance maps ----
fig, axs = plt.subplots(1, 4, figsize=(17, 4.2))
for ax, (side, w) in zip(axs[:2], (("front", "w034"), ("back", "w036"))):
    g, v = outs[side]
    dmap = np.full(v.shape, np.nan)
    dmap[v] = trees[w].query(g[v])[0] * UM
    im = ax.imshow(dmap, cmap="viridis", vmin=0, vmax=max(30, np.nanpercentile(dmap, 99)))
    ax.set_title(f"{side} to {w} (um)")
    ax.set_xlabel("grid col (22.6 um)")
    ax.set_ylabel("grid row (22.6 um, ~z)")
    plt.colorbar(im, ax=ax, fraction=0.046)
for ax, side in zip(axs[2:], ("front", "back")):
    g, v = outs[side]
    for w in ("w034", "w035", "w036"):
        ax.hist(trees[w].query(g[v])[0] * UM, bins=40, range=(0, 250), alpha=0.6, color=cols[w], label=f"to {w}")
    ax.set_title(f"{side}: distance to each published mesh")
    ax.set_xlabel("um")
    ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(PNG, f"{tag}_distmaps.png"), dpi=110)
print("figures written", os.path.join(PNG, f"{tag}_slice.png"))
