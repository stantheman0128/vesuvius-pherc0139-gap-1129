"""Figure for the C240 CT ablation: each panel shows the CT the copy model actually saw (real, constant, shifted
outward, shifted inward) at the same z-slice and zoom, the published meshes (dots), that run's outputs (x) and the
real-CT run's outputs (small +) for reference.
Usage: python -I ablate_fig.py      Writes out/png/abl_copy_panels.png
"""
import json
import os
import site
import sys

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
import tifffile  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "shim"))
from netlog import ROOT  # noqa: E402
import vc  # noqa: E402

UM = 1.129
OUT = os.path.join(ROOT, "out")
A = os.path.join(OUT, "ablate")


def grid(d):
    g = np.stack([tifffile.imread(os.path.join(d, f"{c}.tif")).astype(np.float64) for c in "xyz"], -1)
    return g, np.all(g != -1, -1)


def cut(g, v, z):
    pts = []
    for ax in (0, 1):
        a, b = (g[:-1], g[1:]) if ax == 0 else (g[:, :-1], g[:, 1:])
        va, vb = (v[:-1], v[1:]) if ax == 0 else (v[:, :-1], v[:, 1:])
        ok = va & vb & ((a[..., 2] - z) * (b[..., 2] - z) <= 0) & (a[..., 2] != b[..., 2])
        t = (z - a[..., 2][ok]) / (b[..., 2][ok] - a[..., 2][ok])
        pts.append(a[ok] + t[:, None] * (b[ok] - a[ok]))
    return np.concatenate(pts)


plan = json.load(open(os.path.join(OUT, "plan_ct_l2.json")))
(z0, y0, x0), (z1, y1, x1) = plan["boxes"][0]
pg, pv = grid(os.path.join(OUT, "inputs", "w035_C240_patch"))
P = pg[pv]
zc = int(round(np.median(P[:, 2]) / 4.0))
zl0 = zc * 4.0
bands = {w: np.load(os.path.join(ROOT, "downloads", "mesh_bands", f"{w}_band.npz"))["g"].astype(np.float64)
         for w in ("w034", "w035", "w036")}
panels = [("real CT (no TTA)", os.path.join(ROOT, "downloads", "ct", "vol.zarr"), os.path.join(OUT, "copy_notta")),
          ("(a) constant CT", os.path.join(A, "ct_const_l2", "vol.zarr"), os.path.join(A, "copy_const")),
          ("(b) CT shifted outward 31.6 um", os.path.join(A, "ct_shift_l2_p", "vol.zarr"), os.path.join(A, "copy_shift_p")),
          ("(b') CT shifted inward 31.6 um", os.path.join(A, "ct_shift_l2_m", "vol.zarr"), os.path.join(A, "copy_shift_m"))]
cols = {"w034": "#1f77b4", "w035": "#2ca02c", "w036": "#d62728", "front": "#ff7f0e", "back": "#9467bd"}
ext = [x0 * 4 * UM / 1000, x1 * 4 * UM / 1000, y1 * 4 * UM / 1000, y0 * 4 * UM / 1000]
cx_ = np.median(P[:, 0]) * UM / 1000
cy_ = np.median(P[:, 1]) * UM / 1000
ref = {s: grid(os.path.join(OUT, "copy_notta", f"w035_C240_patch_{s}")) for s in ("front", "back")}
fig, axs = plt.subplots(1, 4, figsize=(20, 5.6))
for ax, (title, vol, run) in zip(axs, panels):
    img = vc.Volume.open(vol).read_zyx([zc, y0, x0], [1, y1 - y0, x1 - x0], level=2)[0]
    ax.imshow(img, cmap="gray", extent=ext, vmin=0, vmax=255)
    for w, b in bands.items():
        c = cut(b, np.all(b != -1, -1), zl0)
        ax.plot(c[:, 0] * UM / 1000, c[:, 1] * UM / 1000, ".", ms=4, color=cols[w], label=f"{w} (published)")
    for s in ("front", "back"):
        g, v = grid(os.path.join(run, f"w035_C240_patch_{s}"))
        c = cut(g, v, zl0)
        ax.plot(c[:, 0] * UM / 1000, c[:, 1] * UM / 1000, "x", ms=7, mew=1.5, color=cols[s], label=f"output {s}")
        if run != os.path.join(OUT, "copy_notta"):
            c = cut(*ref[s], zl0)
            ax.plot(c[:, 0] * UM / 1000, c[:, 1] * UM / 1000, "+", ms=5, color="white", alpha=0.8,
                    label="real-CT output" if s == "front" else None)
    ax.set_xlim(cx_ - 0.35, cx_ + 0.35)
    ax.set_ylim(cy_ + 0.3, cy_ - 0.3)
    ax.set_title(title)
    ax.set_xlabel("x (mm, 1.129 um frame)")
axs[0].set_ylabel(f"y (mm); z = {zl0:.0f}, CT level 2; outward = up")
axs[2].legend(loc="lower left", fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "png", "abl_copy_panels.png"), dpi=100)
print("written", os.path.join(OUT, "png", "abl_copy_panels.png"))
