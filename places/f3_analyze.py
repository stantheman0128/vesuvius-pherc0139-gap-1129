# Step 3: the pilot's layer tracing (p5_analyze method A and B) for any place, plus
#   - w034 -> w036 run-index difference per dz row (2 = one sheet between, 1 = adjacent sheets, 0 = same band, <0 = reversed)
#   - the same with all meshes shifted by -30 / +30 um in r (registration robustness)
#   - each doubly traced end: run-index difference from its own inner piece, per dz row
#   - figures: radial cut (plain / with bodies and runs), profile at dz 0, slice at dz 0, slices at dz -150 / 0 / +150
import sys, os, json; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from repo_paths import OUT_PLACES as OUT
import numpy as np
from collections import Counter
from scipy import ndimage as ndi
from scipy.signal import find_peaks
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
COL = {"w033": "#1f77b4", "w034": "#d62728", "w035": "#9467bd", "w036": "#2ca02c", "w037": "#ff7f0e"}
MINRUN, MINGAP = 6, 3

def otsu(v):
    v = v[np.isfinite(v)]
    h, e = np.histogram(v, 256); c = (e[:-1] + e[1:]) / 2
    w0 = np.cumsum(h); w1 = w0[-1] - w0; m0 = np.cumsum(h * c)
    mu0 = m0 / np.maximum(w0, 1); mu1 = (m0[-1] - m0) / np.maximum(w1, 1)
    return c[np.argmax(w0 * w1 * (mu0 - mu1) ** 2)]

def runs(row):
    r = ndi.binary_closing(row, np.ones(MINGAP + 1, bool)) | row
    lab, n = ndi.label(r)
    return [(s[0].start, s[0].stop) for s in ndi.find_objects(lab) if s[0].stop - s[0].start >= MINRUN]

def pieces(h):
    out = {}
    for dz in np.unique(h[:, 1]):
        rs = np.sort(h[h[:, 1] == dz, 0])
        out.setdefault("main", []).append((rs[0], dz))
        if len(rs) > 1:
            out.setdefault("end", []).append((rs[-1], dz))
    return {k: np.array(v) for k, v in out.items()}

def run_of(r, rr):
    for i, (a, b) in enumerate(rr):
        if a - 0.5 <= r <= b + 0.5: return i, 0.0
    if not rr: return None, np.nan
    ds = [min(abs(r - a), abs(r - b)) for a, b in rr]
    i = int(np.argmin(ds)); return i, ds[i]

for name in sys.argv[1:]:
    D = np.load(os.path.join(OUT, f"{name}.npz"))
    meta = json.loads(str(D["meta"]))
    z0 = meta["z0"]
    cut, rc, zz = D["cut"], D["rc"], D["zz"]
    sm = ndi.gaussian_filter(cut, 2)
    thr = otsu(sm)
    mask = ndi.binary_opening(sm > thr, np.ones((3, 3), bool))
    lab, nlab = ndi.label(mask, np.ones((3, 3)))
    sizes = ndi.sum(mask, lab, range(1, nlab + 1))
    keep = np.isin(lab, 1 + np.where(sizes >= 200)[0])
    lab = np.where(keep, lab, 0); mask = keep

    M = {}
    for w in COL:
        h = D[f"ch_{w}"]
        if len(h) == 0: continue
        for k, v in pieces(h).items():
            if len(v) >= 1: M[w + ("" if k == "main" else " end")] = v

    def pos(key, dz):
        if key not in M: return np.nan
        v = M[key]; i = np.argmin(np.abs(v[:, 1] - dz))
        return v[i, 0] if abs(v[i, 1] - dz) < 1 else np.nan

    rows = []
    for dz in np.arange(-200, 201, 10.0):
        k = int(np.argmin(np.abs(zz - dz)))
        rows.append((dz, [(rc[a], rc[b - 1]) for a, b in runs(mask[k])]))

    rep = {"name": name, "z0": z0, "angle": meta["angle"], "otsu": float(thr), "per_dz": []}
    for dz, rr in rows:
        ent = {"dz": dz, "runs": [(round(a), round(b)) for a, b in rr], "mesh": {}}
        for key in M:
            r = pos(key, dz)
            if np.isfinite(r):
                i, dist = run_of(r, rr)
                ent["mesh"][key] = {"r": round(float(r), 1), "run": i, "dist_um": round(float(dist), 1)}
        # w034 -> w036, with shifts
        for sh in (-30, 0, 30):
            r4, r6 = pos("w034", dz), pos("w036", dz)
            if np.isfinite(r4) and np.isfinite(r6) and rr:
                i4, d4 = run_of(r4 + sh, rr); i6, d6 = run_of(r6 + sh, rr)
                ent[f"d46_{sh}"] = i6 - i4
                if sh == 0: ent["dist46"] = (round(d4, 1), round(d6, 1))
        # doubly traced ends vs own inner piece
        for w in COL:
            if f"{w} end" in ent["mesh"] and w in ent["mesh"]:
                ent[f"enddiff_{w}"] = ent["mesh"][f"{w} end"]["run"] - ent["mesh"][w]["run"]
                # which other mesh shares the end's run
                same = [k for k, v in ent["mesh"].items() if k != f"{w} end" and v["run"] == ent["mesh"][f"{w} end"]["run"]]
                ent[f"endshares_{w}"] = same
        # w033 -> w037 run difference
        r3, r7 = pos("w033", dz), pos("w037", dz)
        if np.isfinite(r3) and np.isfinite(r7) and rr:
            ent["d37"] = run_of(r7, rr)[0] - run_of(r3, rr)[0]
        rep["per_dz"].append(ent)

    sel = np.abs(zz) <= 5
    prof = np.nanmean(cut[sel], 0); ps = ndi.gaussian_filter1d(prof, 2)
    rng = np.percentile(ps, 99) - np.percentile(ps, 1)
    pk, props = find_peaks(ps, prominence=0.15 * rng, distance=10)
    rep["profile_peaks_um"] = [round(float(rc[i])) for i in pk]
    rep["profile_peak_prominence"] = [round(float(p), 1) for p in props["prominences"]]
    k0 = int(np.argmin(np.abs(zz)))
    rep["runs_dz0"] = [(round(rc[a]), round(rc[b - 1])) for a, b in runs(mask[k0])]
    for sh in (-30, 0, 30):
        rep[f"d46_hist_{sh}"] = dict(Counter(e.get(f"d46_{sh}") for e in rep["per_dz"] if f"d46_{sh}" in e))
    rep["d37_hist"] = dict(Counter(e.get("d37") for e in rep["per_dz"] if "d37" in e))
    for w in COL:
        c = Counter(e[f"enddiff_{w}"] for e in rep["per_dz"] if f"enddiff_{w}" in e)
        if c: rep[f"enddiff_hist_{w}"] = dict(c)
    # profile at dz 0 between w034 and w036
    r4, r6 = pos("w034", 0.0), pos("w036", 0.0)
    if np.isfinite(r4) and np.isfinite(r6):
        a_, b_ = sorted((r4, r6)); m_ = (rc >= a_) & (rc <= b_)
        rep["prof46_dz0"] = {"min": round(float(ps[m_].min()), 1), "max": round(float(ps[m_].max()), 1), "thr": round(float(thr), 1)}
    json.dump(rep, open(os.path.join(OUT, f"{name}_layers.json"), "w"), indent=1, default=lambda x: x if not isinstance(x, np.generic) else x.item())

    # ---------------- figures ----------------
    zlab = f"dz from z {z0:,.0f} (um)"
    fig, axs = plt.subplots(2, 1, figsize=(16, 15))
    lo, hi = np.nanpercentile(cut, [1, 99.5])
    ext = [rc[0] - 0.5, rc[-1] + 0.5, zz[0] - 0.5, zz[-1] + 0.5]
    for ax, overlay in zip(axs, (False, True)):
        ax.imshow(cut, cmap="gray", vmin=lo, vmax=hi, extent=ext, origin="lower", interpolation="nearest")
        if overlay:
            rgb = plt.cm.tab20(lab % 20 / 20.0); rgb[..., 3] = np.where(lab > 0, 0.35, 0)
            ax.imshow(rgb, extent=ext, origin="lower", interpolation="nearest")
            ax.contour(rc, zz, mask.astype(float), [0.5], colors="cyan", linewidths=0.6)
            for dz, rr in rows:
                for a, b in rr: ax.plot([a, b], [dz, dz], color="yellow", lw=1.2)
        for key, v in M.items():
            w = key.split()[0]
            ax.plot(v[:, 0], v[:, 1], "o" if "end" not in key else "s", color=COL[w], ms=3.5,
                    mfc="none" if "end" in key else COL[w], label=key + (" (doubly traced end)" if "end" in key else ""))
        ax.axhline(0, color="white", lw=0.6, ls=":")
        ax.set_xlabel("r from umbilicus along the radial line (um)"); ax.set_ylabel(zlab); ax.set_aspect("equal")
        ax.legend(loc="upper left", fontsize=8, framealpha=0.8)
    axs[0].set_title(f"{name}: radial cut at {meta['angle']:.1f} deg through the 1.129 um scan (mean of +-3 um across the line); mesh crossings marked")
    axs[1].set_title(f"{name}: same, automatic material bodies (Otsu {thr:.0f} on sigma-2 um smoothing, cyan) and runs every 10 um in dz (yellow)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{name}_radialcut.png"), dpi=100); plt.close(fig)

    fig, ax = plt.subplots(figsize=(16, 5))
    ax.plot(rc, prof, color="0.7", lw=0.7, label="mean of dz -5..+5 um, +-3 um across")
    ax.plot(rc, ps, color="k", lw=1.2, label="Gaussian sigma 2 um")
    ax.axhline(thr, color="c", lw=0.8, ls="--", label=f"Otsu threshold of the cut ({thr:.0f})")
    ax.plot(rc[pk], ps[pk], "v", color="m", ms=8, label="peaks (prominence >= 15% of range)")
    for a, b in rep["runs_dz0"]: ax.axvspan(a, b, color="yellow", alpha=0.25)
    for key in M:
        r = pos(key, 0.0); w = key.split()[0]
        if np.isfinite(r):
            ax.axvline(r, color=COL[w], lw=1.5, ls="-" if "end" not in key else "--")
            ax.text(r, np.nanmax(ps) * 1.02, key, color=COL[w], rotation=90, va="bottom", ha="center", fontsize=8)
    ax.set_xlabel("r from umbilicus (um)"); ax.set_ylabel("CT value (uint8)")
    ax.set_title(f"{name}: profile along the radial line at z {z0:,.0f}, angle {meta['angle']:.1f} deg (yellow = material runs at dz 0)")
    ax.legend(loc="lower right", fontsize=8); ax.set_ylim(top=np.nanmax(ps) * 1.25)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{name}_profile.png"), dpi=100); plt.close(fig)

    rr_, tt_ = D["rr"], D["tt"]
    def slice_ax(ax, dz, mark):
        img = D[f"slice_{dz}"]
        lo, hi = np.nanpercentile(img, [1, 99.5])
        ax.imshow(img, cmap="gray", vmin=lo, vmax=hi, extent=[rr_[0] - .5, rr_[-1] + .5, tt_[0] - .5, tt_[-1] + .5], origin="lower", interpolation="nearest")
        for w, c in COL.items():
            k = f"tr_{w}_{dz}"
            if k in D:
                for seg in D[k]: ax.plot(seg[:, 0], seg[:, 1], color=c, lw=0.9)
            ax.plot([], [], color=c, label=w)
        if mark:
            for key in M:
                r = pos(key, float(dz))
                if np.isfinite(r):
                    ax.plot(r, 0, "o" if "end" not in key else "s", mfc="none", mec=COL[key.split()[0]], ms=9, mew=1.5)
        ax.axhline(0, color="yellow", lw=0.7, ls="--")
        ax.set_xlim(rr_[0], rr_[-1]); ax.set_ylim(tt_[0], tt_[-1]); ax.set_aspect("equal")
        ax.set_title(f"dz {dz:+d} um (z {z0 + dz/2.399:,.0f} in the 2.399 um frame)")
        ax.set_xlabel("r from umbilicus (um)"); ax.set_ylabel("tangential (um)")
    fig, ax = plt.subplots(figsize=(14, 8))
    slice_ax(ax, 0, True); ax.legend(loc="upper left")
    ax.set_title(f"{name}: z slice {z0:,.0f} in the 1.129 um scan, mean of 5 planes; radial line dashed; crossings circled (square = doubly traced end)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{name}_slice.png"), dpi=100); plt.close(fig)
    fig, axs = plt.subplots(3, 1, figsize=(14, 20))
    for ax, dz in zip(axs, (-150, 0, 150)): slice_ax(ax, dz, True)
    axs[0].legend(loc="upper left")
    fig.suptitle(f"{name}: z slices of the 1.129 um scan at three heights with mesh traces at those heights; radial line dashed")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{name}_slices_3heights.png"), dpi=90); plt.close(fig)

    print(f"== {name}: angle {meta['angle']:.2f}, z {z0:.0f}; Otsu {thr:.1f}; peaks {rep['profile_peaks_um']}")
    print(f"   runs dz0: {rep['runs_dz0']}")
    print(f"   d46 hist (w036 run - w034 run; 2 = one sheet between): shift-30 {rep['d46_hist_-30']}  0 {rep['d46_hist_0']}  +30 {rep['d46_hist_30']}; d37 {rep['d37_hist']}")
    print(f"   prof46 {rep.get('prof46_dz0')}; ends: " + str({k: v for k, v in rep.items() if k.startswith('enddiff_hist')}))
    for e in rep["per_dz"]:
        ms = "  ".join(f"{k}:{v['r']:.0f}->{v['run']}" + (f"(+{v['dist_um']:.0f})" if v['dist_um'] else "")
                       for k, v in sorted(e["mesh"].items(), key=lambda kv: kv[1]["r"]))
        print(f"   dz {e['dz']:+5.0f}: {len(e['runs'])} runs d46 {e.get('d46_0')} | {ms}")
