# Step 4: trace the layers in the CT (not from the meshes) and count sheet crossings between mesh positions.
# Method A (2D, radial cut): Gaussian smoothing (sigma 2 um), Otsu threshold over the cut, small specks removed,
#   connected components = material bodies; along each dz row the material runs between two mesh crossings are counted.
# Method B (1D, profile on the ray): mean of rows dz -5..+5 um, Gaussian sigma 2 um, peaks with prominence >= 15%
#   of the profile range and >= 10 um apart.
import sys, os, json; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from repo_paths import OUT_PILOT as OUT
import numpy as np
from scipy import ndimage as ndi
from scipy.signal import find_peaks
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
COL = {"w033": "#1f77b4", "w034": "#d62728", "w035": "#9467bd", "w036": "#2ca02c", "w037": "#ff7f0e"}
MINRUN, MINGAP = 6, 3  # um: material runs shorter than MINRUN dropped, air gaps shorter than MINGAP closed

def otsu(v):
    v = v[np.isfinite(v)]
    h, e = np.histogram(v, 256); c = (e[:-1] + e[1:]) / 2
    w0 = np.cumsum(h); w1 = w0[-1] - w0; m0 = np.cumsum(h * c)
    mu0 = m0 / np.maximum(w0, 1); mu1 = (m0[-1] - m0) / np.maximum(w1, 1)
    return c[np.argmax(w0 * w1 * (mu0 - mu1) ** 2)]

def runs(row):
    """material runs [start, end) of a boolean row, after closing short gaps and dropping short runs"""
    r = ndi.binary_closing(row, np.ones(MINGAP + 1, bool)) | row
    lab, n = ndi.label(r)
    sl = ndi.find_objects(lab)
    return [(s[0].start, s[0].stop) for s in sl if s[0].stop - s[0].start >= MINRUN]

def pieces(h):
    """split per-dz ray hits of a mesh into its inner piece and its outer (doubly traced) piece"""
    out = {}
    for dz in np.unique(h[:, 1]):
        rs = np.sort(h[h[:, 1] == dz, 0])
        out.setdefault("main", []).append((rs[0], dz))
        if len(rs) > 1:
            out.setdefault("end", []).append((rs[-1], dz))
    return {k: np.array(v) for k, v in out.items()}

summary = {}
for name in sys.argv[1:]:
    D = np.load(os.path.join(OUT, f"{name}.npz"))
    cut, rc, zz = D["cut"], D["rc"], D["zz"]
    sm = ndi.gaussian_filter(cut, 2)
    thr = otsu(sm)
    mask = sm > thr
    mask = ndi.binary_opening(mask, np.ones((3, 3), bool))
    lab, nlab = ndi.label(mask, np.ones((3, 3)))
    sizes = ndi.sum(mask, lab, range(1, nlab + 1))
    keep = np.isin(lab, 1 + np.where(sizes >= 200)[0])
    lab = np.where(keep, lab, 0)

    # mesh positions per dz
    M = {}
    for w in COL:
        h = D[f"ch_{w}"]
        if len(h) == 0: continue
        for k, v in pieces(h).items():
            if len(v) >= 5: M[f"{w}" + ("" if k == "main" else " end")] = v

    def pos(key, dz):
        v = M[key]; i = np.argmin(np.abs(v[:, 1] - dz))
        return v[i, 0] if abs(v[i, 1] - dz) < 1 else np.nan

    # per-dz material runs and their relation to mesh positions
    rows = []
    for dz in np.arange(-200, 201, 10.0):
        k = int(np.argmin(np.abs(zz - dz)))
        rr = [(rc[a], rc[b - 1]) for a, b in runs(mask[k])]
        rows.append((dz, rr))

    def run_of(r, rr):
        for i, (a, b) in enumerate(rr):
            if a - 0.5 <= r <= b + 0.5: return i, 0.0
        ds = [min(abs(r - a), abs(r - b)) for a, b in rr]
        i = int(np.argmin(ds)); return i, ds[i]

    report = {"otsu": float(thr), "per_dz": []}
    for dz, rr in rows:
        ent = {"dz": dz, "runs": [(round(a), round(b)) for a, b in rr], "mesh": {}}
        for key in M:
            r = pos(key, dz)
            if np.isfinite(r):
                i, dist = run_of(r, rr)
                ent["mesh"][key] = {"r": round(float(r), 1), "run": i, "dist_to_material_um": round(float(dist), 1)}
        report["per_dz"].append(ent)

    # profile on the ray (dz -5..+5)
    sel = np.abs(zz) <= 5
    prof = np.nanmean(cut[sel], 0)
    ps = ndi.gaussian_filter1d(prof, 2)
    rng = np.percentile(ps, 99) - np.percentile(ps, 1)
    pk, props = find_peaks(ps, prominence=0.15 * rng, distance=10)
    report["profile_peaks_um"] = [round(float(rc[i])) for i in pk]
    report["profile_peak_prominence"] = [round(float(p), 1) for p in props["prominences"]]
    k0 = int(np.argmin(np.abs(zz)))
    report["runs_dz0"] = [(round(rc[a]), round(rc[b - 1])) for a, b in runs(mask[k0])]
    summary[name] = report
    json.dump(report, open(os.path.join(OUT, f"{name}_layers.json"), "w"), indent=1)

    # ---------- figures ----------
    # radial cut with material outline, component labels and mesh positions
    fig, axs = plt.subplots(2, 1, figsize=(16, 15), gridspec_kw={"height_ratios": [1, 1]})
    lo, hi = np.nanpercentile(cut, [1, 99.5])
    ext = [rc[0] - 0.5, rc[-1] + 0.5, zz[0] - 0.5, zz[-1] + 0.5]
    for ax, overlay in zip(axs, (False, True)):
        ax.imshow(cut, cmap="gray", vmin=lo, vmax=hi, extent=ext, origin="lower", interpolation="nearest")
        if overlay:
            rgb = plt.cm.tab20(lab % 20 / 20.0); rgb[..., 3] = np.where(lab > 0, 0.35, 0)
            ax.imshow(rgb, extent=ext, origin="lower", interpolation="nearest")
            ax.contour(rc, zz, mask.astype(float), [0.5], colors="cyan", linewidths=0.6)
            for dz, rr in rows:
                for a, b in rr:
                    ax.plot([a, b], [dz, dz], color="yellow", lw=1.2)
        for key, v in M.items():
            w = key.split()[0]
            ax.plot(v[:, 0], v[:, 1], "o" if "end" not in key else "s", color=COL[w], ms=3.5, mfc="none" if "end" in key else COL[w],
                    label=key + (" (doubly traced end)" if "end" in key else ""))
        ax.axhline(0, color="white", lw=0.6, ls=":")
        ax.set_xlabel("r from umbilicus along the radial line (um)"); ax.set_ylabel("dz from z 17,995 (um)")
        ax.set_aspect("equal")
        ax.legend(loc="upper left", fontsize=8, framealpha=0.8)
    axs[0].set_title(f"{name}: radial cut through the 1.129 um scan (mean of +-3 um across the line); mesh crossings marked")
    axs[1].set_title(f"{name}: same, with automatic material bodies (Otsu {thr:.0f} on sigma-2 um smoothing, cyan outline, "
                     f"tinted by connected component) and material runs every 10 um in dz (yellow)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{name}_radialcut.png"), dpi=100); plt.close(fig)

    # profile
    fig, ax = plt.subplots(figsize=(16, 5))
    ax.plot(rc, prof, color="0.7", lw=0.7, label="mean of dz -5..+5 um, +-3 um across")
    ax.plot(rc, ps, color="k", lw=1.2, label="Gaussian sigma 2 um")
    ax.axhline(thr, color="c", lw=0.8, ls="--", label=f"Otsu threshold of the cut ({thr:.0f})")
    ax.plot(rc[pk], ps[pk], "v", color="m", ms=8, label="peaks (prominence >= 15% of range)")
    for a, b in report["runs_dz0"]:
        ax.axvspan(a, b, color="yellow", alpha=0.25)
    for key in M:
        r = pos(key, 0.0)
        w = key.split()[0]
        if np.isfinite(r):
            ax.axvline(r, color=COL[w], lw=1.5, ls="-" if "end" not in key else "--")
            ax.text(r, ax.get_ylim()[1] if False else np.nanmax(ps) * 1.02, key, color=COL[w], rotation=90, va="bottom", ha="center", fontsize=8)
    ax.set_xlabel("r from umbilicus (um)"); ax.set_ylabel("CT value (uint8)")
    ax.set_title(f"{name}: intensity profile along the radial line at z 17,995 (yellow = material runs at dz 0)")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_ylim(top=np.nanmax(ps) * 1.25)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{name}_profile.png"), dpi=100); plt.close(fig)

    # slice with traces and crossings on the line
    s = D["slice"]; rr_, tt_ = D["rr"], D["tt"]
    lo, hi = np.nanpercentile(s, [1, 99.5])
    fig, ax = plt.subplots(figsize=(13, 13))
    ax.imshow(s, cmap="gray", vmin=lo, vmax=hi, extent=[rr_[0] - .5, rr_[-1] + .5, tt_[0] - .5, tt_[-1] + .5], origin="lower", interpolation="nearest")
    for w, c in COL.items():
        for seg in D[f"tr_{w}"]:
            ax.plot(seg[:, 0], seg[:, 1], color=c, lw=0.9)
        ax.plot([], [], color=c, label=w)
    for key in M:
        r = pos(key, 0.0)
        if np.isfinite(r):
            ax.plot(r, 0, "o" if "end" not in key else "s", mfc="none", mec=COL[key.split()[0]], ms=9, mew=1.5)
    ax.axhline(0, color="yellow", lw=0.7, ls="--")
    ax.set_xlim(rr_[0], rr_[-1]); ax.set_ylim(tt_[0], tt_[-1]); ax.set_aspect("equal")
    ax.set_xlabel("r from umbilicus (um)"); ax.set_ylabel("tangential offset (um, + = counter-clockwise)")
    ax.legend(loc="upper left")
    ax.set_title(f"{name}: z slice 17,995 (2.399 um frame) in the 1.129 um scan, mean of 5 planes; radial line dashed, crossings circled (square = doubly traced end)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{name}_slice.png"), dpi=100); plt.close(fig)

for name, rep in summary.items():
    print(f"== {name}: Otsu {rep['otsu']:.1f}; profile peaks {rep['profile_peaks_um']} (prom {rep['profile_peak_prominence']})")
    print(f"   material runs at dz 0: {rep['runs_dz0']}")
    for ent in rep["per_dz"]:
        ms = "  ".join(f"{k}:{v['r']:.0f}->run{v['run']}" + (f"(+{v['dist_to_material_um']:.0f})" if v['dist_to_material_um'] else "")
                       for k, v in sorted(ent["mesh"].items(), key=lambda kv: kv[1]["r"]))
        print(f"   dz {ent['dz']:+5.0f}: {len(ent['runs'])} runs | {ms}")
