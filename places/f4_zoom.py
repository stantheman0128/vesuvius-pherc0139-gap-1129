# Zoomed views between w033 and w037: slices at dz -150/0/+150 (left: plain, right: with traces) and the radial cut,
# r window = [w033 - 60, w037 + 60] at dz 0, tang +-200 um. usage: f4_zoom.py NAME [R0 R1 [SUFFIX]]
# (SUFFIX is appended to the output names, e.g. f4_zoom.py X542 9640 10060 _mid)
import sys, os, json; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from repo_paths import OUT_PLACES as OUT
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
COL = {"w033": "#1f77b4", "w034": "#d62728", "w035": "#9467bd", "w036": "#2ca02c", "w037": "#ff7f0e"}
name = sys.argv[1]
SUF = sys.argv[4] if len(sys.argv) > 4 else ""
D = np.load(os.path.join(OUT, f"{name}.npz")); meta = json.loads(str(D["meta"]))
h = meta["hits24"]
allr = [r for w in COL for r in h.get(w, [])]
R0 = (float(sys.argv[2]) if len(sys.argv) > 2 else min(allr) - 60); R1 = (float(sys.argv[3]) if len(sys.argv) > 3 else max(allr) + 60)
TW = 200
rr, tt = D["rr"], D["tt"]
fig, axs = plt.subplots(3, 2, figsize=(18, 3 * 18 * (2 * TW) / (2 * (R1 - R0)) + 1.5))
for row, dz in enumerate((-150, 0, 150)):
    img = D[f"slice_{dz}"]
    m = (rr >= R0) & (rr <= R1); n = np.abs(tt) <= TW
    sub = img[np.ix_(n, m)]
    lo, hi = np.nanpercentile(sub, [1, 99.5])
    for col in (0, 1):
        ax = axs[row, col]
        ax.imshow(sub, cmap="gray", vmin=lo, vmax=hi, extent=[rr[m][0] - .5, rr[m][-1] + .5, tt[n][0] - .5, tt[n][-1] + .5], origin="lower", interpolation="nearest")
        if col:
            for w, c in COL.items():
                k = f"tr_{w}_{dz}"
                if k in D:
                    for seg in D[k]: ax.plot(seg[:, 0], seg[:, 1], color=c, lw=0.8)
                ax.plot([], [], color=c, label=w)
            ax.legend(loc="upper left", fontsize=7)
        else:
            ax.axhline(0, color="yellow", lw=0.5, ls=":")
        ax.set_xlim(R0, R1); ax.set_ylim(-TW, TW); ax.set_aspect("equal")
        ax.set_title(f"{name} dz {dz:+d} um" + (" with traces" if col else ""), fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{name}_slices_zoom{SUF}.png"), dpi=110); plt.close(fig)

cut, rc, zz = D["cut"], D["rc"], D["zz"]
m = (rc >= R0) & (rc <= R1); sub = cut[:, m]
lo, hi = np.nanpercentile(sub, [1, 99.5])
fig, axs = plt.subplots(1, 2, figsize=(18, 18 * 400 / (2 * (R1 - R0)) + 1))
for i, ax in enumerate(axs):
    ax.imshow(sub, cmap="gray", vmin=lo, vmax=hi, extent=[rc[m][0] - .5, rc[m][-1] + .5, zz[0] - .5, zz[-1] + .5], origin="lower", interpolation="nearest")
    if i:
        for w in COL:
            v = D[f"ch_{w}"]
            if len(v): ax.plot(v[:, 0], v[:, 1], ".", color=COL[w], ms=3)
    ax.set_aspect("equal"); ax.set_xlabel("r (um)"); ax.set_ylabel("dz (um)")
fig.suptitle(f"{name}: radial cut zoom (left plain, right with mesh crossings)")
fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{name}_radialcut_zoom{SUF}.png"), dpi=110); plt.close(fig)
