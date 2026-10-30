# Figure: E2 z slices at dz -150, 0, +150 um with mesh traces at those heights and the radial line.
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from repo_paths import OUT_PILOT as OUT
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
COL = {"w033": "#1f77b4", "w034": "#d62728", "w036": "#2ca02c", "w037": "#ff7f0e"}
L = np.load(os.path.join(OUT, "E2_levels.npz")); D = np.load(os.path.join(OUT, "E2.npz"))
rr, tt = D["rr"], D["tt"]
panels = [(-150, L["img_-150"], {w: L[f"tr_-150_{w}"] for w in COL}),
          (0, D["slice"], {w: D[f"tr_{w}"] for w in COL}),
          (150, L["img_150"], {w: L[f"tr_150_{w}"] for w in COL})]
R0, R1, T0, T1 = 9350, 10300, -450, 450
fig, axs = plt.subplots(1, 3, figsize=(24, 11.5))
for ax, (dz, img, tr) in zip(axs, panels):
    lo, hi = np.nanpercentile(img, [1, 99.5])
    ax.imshow(img, cmap="gray", vmin=lo, vmax=hi, extent=[rr[0] - .5, rr[-1] + .5, tt[0] - .5, tt[-1] + .5], origin="lower", interpolation="nearest")
    for w, c in COL.items():
        for seg in tr[w]:
            ax.plot(seg[:, 0], seg[:, 1], color=c, lw=1.0, alpha=0.9)
        ax.plot([], [], color=c, label=w)
    ax.axhline(0, color="yellow", lw=0.7, ls="--")
    ax.set_xlim(R0, R1); ax.set_ylim(T0, T1); ax.set_aspect("equal")
    ax.set_title(f"dz {dz:+d} um (z {17995 + dz/2.399:.0f} in the 2.399 um frame)")
    ax.set_xlabel("r from umbilicus (um)")
axs[0].set_ylabel("tangential offset (um)"); axs[0].legend(loc="upper left")
fig.suptitle("E2: z slices of the 1.129 um scan at three heights; mesh traces at the same heights; radial line dashed (w034 and w033 show both pieces)")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "E2_slices_3heights.png"), dpi=90); plt.close(fig)
