# Step 3b: quick-look renders of the slice and the radial cut with mesh traces (no layer picking yet).
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from repo_paths import OUT_PILOT as OUT
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
COL = {"w033": "#1f77b4", "w034": "#d62728", "w035": "#9467bd", "w036": "#2ca02c", "w037": "#ff7f0e"}

for name in sys.argv[1:]:
    D = np.load(os.path.join(OUT, f"{name}.npz"))
    s = D["slice"]; rr, tt = D["rr"], D["tt"]
    lo, hi = np.nanpercentile(s, [1, 99.5])
    fig, ax = plt.subplots(figsize=(13, 13))
    ax.imshow(s, cmap="gray", vmin=lo, vmax=hi, extent=[rr[0], rr[-1], tt[-1], tt[0]], interpolation="nearest")
    for w, c in COL.items():
        tr = D[f"tr_{w}"]
        for seg in tr:
            ax.plot(seg[:, 0], seg[:, 1], color=c, lw=0.8)
        ax.plot([], [], color=c, label=w)
    ax.axhline(0, color="yellow", lw=0.5, ls="--")
    ax.set_xlabel("r from umbilicus (um)"); ax.set_ylabel("tangential offset (um)")
    ax.legend(loc="upper left"); ax.set_title(f"{name}: z slice, 1.129 um scan, radial line = dashed")
    fig.savefig(os.path.join(OUT, f"look_{name}_slice.png"), dpi=110, bbox_inches="tight"); plt.close(fig)
    for key in ("cut", "cut1"):
        c_ = D[key]; rc, zz = D["rc"], D["zz"]
        lo, hi = np.nanpercentile(c_, [1, 99.5])
        fig, ax = plt.subplots(figsize=(16, 8))
        ax.imshow(c_, cmap="gray", vmin=lo, vmax=hi, extent=[rc[0], rc[-1], zz[-1], zz[0]], interpolation="nearest", origin="upper")
        for w, c in COL.items():
            h = D[f"ch_{w}"]
            if len(h): ax.plot(h[:, 0], h[:, 1], ".", color=c, ms=4, label=w)
        ax.invert_yaxis()
        ax.set_xlabel("r (um)"); ax.set_ylabel("dz (um)"); ax.legend(loc="upper left")
        fig.savefig(os.path.join(OUT, f"look_{name}_{key}.png"), dpi=100, bbox_inches="tight"); plt.close(fig)
