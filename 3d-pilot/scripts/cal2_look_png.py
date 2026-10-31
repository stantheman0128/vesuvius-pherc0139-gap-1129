"""Cross-sections of the cal2 CT (level 3, trilinear) on ribbons that contain the local w035 outward normal at every
image column (cal2_look.ribbon), with the w034 / w035 / w036 crossings of each column's line and the real-CT trace.
  step ribbon:   along the trace front at that step (37 columns, about z), split into two strips.
  column ribbon: along one trace column through steps 0 to 8 (along the arc). If out/cal2_look/ridge_follow.json
                 exists, the CT-only ridge path for that column (cal2_look_ridge.py) is drawn in green.
Each image column is the CT along that point's own normal, so the trace sits exactly on the centre line at its
points and the mesh marks are exact line crossings. Outward (toward w036) is up. Ribbons bend with the sheet, so
distances along the normal are true but distances along the ribbon are only true on the centre line.
Only cal2's 30 chunks are read (cal2_look.load_plan / load_ct).
Usage (WSL venv): python -I cal2_look_png.py [--steps 2,5,8] [--cols 12] [--out DIR]
Writes DIR/cal2_step<i>.png, DIR/cal2_col<c>.png (default DIR out/cal2_look/png) and DIR/figures.json.
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cal2_look as L  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402

DU = 2.0                    # level-0 voxels between image columns (2.26 um)
UP = 2                      # nearest upscaling for display
S_FIG = 250.0               # um each side of the ribbon
COL = {"w034": (0, 200, 255), "w035": (255, 220, 0), "w036": (255, 60, 255), "trace": (255, 40, 40),
       "ridge": (40, 220, 40), "f2": (255, 150, 150), "f3": (255, 140, 0)}
LO, HI = 20.0, 200.0        # grey window (CT counts)


def font(sz):
    try:
        return ImageFont.load_default(size=sz)
    except TypeError:
        return ImageFont.load_default()


def strip(img, cross, marks, ridge=None):
    """One ribbon image (rows = s from +S_FIG down, cols = u) with overlays; returns a PIL image."""
    g = np.clip((img.T[::-1] - LO) / (HI - LO), 0, 1) * 255
    g = np.where(np.isnan(img.T[::-1]), 90, g).astype(np.uint8)
    H, W = g.shape
    im = Image.fromarray(np.repeat(g[..., None], 3, -1)).resize((W * UP, H * UP), Image.NEAREST)
    dr = ImageDraw.Draw(im)
    ds = 2 * S_FIG / (H - 1)

    def yy(s):
        return (S_FIG - s) / ds * UP + UP / 2

    for j, cr in enumerate(cross):
        for w in L.MESHES:
            for x in cr[w]:
                if abs(x) <= S_FIG:
                    dr.rectangle((j * UP, yy(x) - 1, j * UP + UP - 1, yy(x)), fill=COL[w])
    if ridge is not None:
        for j, x, tag in ridge:
            if x is not None and abs(x) <= S_FIG:
                off = {"ridge": 0, "f2": 3, "f3": -3}[tag]       # small offsets so coincident paths stay visible
                dr.rectangle((j * UP, yy(x) - 1 + off, j * UP + UP - 1, yy(x) + off), fill=COL[tag])
    for j, lab in marks:
        x, y, r = j * UP + UP / 2, yy(0.0), 6
        dr.line((x - r, y, x + r, y), fill=COL["trace"], width=2)
        dr.line((x, y - r, x, y + r), fill=COL["trace"], width=2)
        if lab:
            dr.text((x + 4, y - 18), lab, fill=COL["trace"], font=font(13))
    return im


def frame(strips, title, legend_extra, path):
    pad_top, gap, pad_bot = 46, 10, 34
    W = max(s.width for s in strips)
    Hs = sum(s.height for s in strips) + gap * (len(strips) - 1)
    canvas = Image.new("RGB", (W, pad_top + Hs + pad_bot), (255, 255, 255))
    y = pad_top
    for s in strips:
        canvas.paste(s, (0, y))
        y += s.height + gap
    dr = ImageDraw.Draw(canvas)
    dr.text((6, 4), title, fill=(0, 0, 0), font=font(15))
    xl = 6
    items = [("w034", "w034"), ("w035", "w035"), ("w036", "w036"), ("trace", "real-CT trace (+)")] + legend_extra
    for k, lab in items:
        dr.rectangle((xl, 28, xl + 10, 38), fill=COL[k])
        dr.text((xl + 14, 25), lab, fill=(0, 0, 0), font=font(13))
        xl += 24 + 8 * len(lab)
    dr.text((xl + 10, 25), "up = outward (toward w036); each column = CT along that point's normal",
            fill=(0, 0, 0), font=font(13))
    y0 = pad_top + Hs + 16
    bar = 100.0 / L.UM / DU * UP
    dr.line((10, y0, 10 + bar, y0), fill=(0, 0, 0), width=3)
    dr.text((16 + bar, y0 - 8), "100 um (both axes)", fill=(0, 0, 0), font=font(13))
    canvas.save(path)
    return [canvas.width, canvas.height]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--steps", default="2,5,8")
    ap.add_argument("--cols", default="12")
    ap.add_argument("--out", default=os.path.join(L.OUT, "png"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    p, keys = L.load_plan()
    arr, lo = L.load_ct(keys)
    ct = L.CT(arr, lo)
    meshes, normal_at = L.load_meshes(p)
    P, kept, _ = L.load_trace("real")
    s = np.linspace(S_FIG, -S_FIG, int(round(2 * S_FIG / (DU * L.UM))) + 1)[::-1]
    tmax = S_FIG / L.UM
    rf_path = os.path.join(L.OUT, "ridge_follow.json")
    rf = json.load(open(rf_path)) if os.path.exists(rf_path) else None
    info = {}

    def crossings(base, n):
        return [{w: [t * L.UM for t in meshes[w].cross(b, nn, tmax)] for w in L.MESHES} for b, nn in zip(base, n)]

    for i in [int(x) for x in a.steps.split(",") if x]:
        cols = np.nonzero(kept[i])[0]
        u, cum, base, n, img = L.ribbon(ct, normal_at, P[i][cols], DU, 20 * 8, s, masked=True)
        cr = crossings(base, n)
        marks = [(int(np.argmin(np.abs(u - cu))), (str(c) if c % 6 == 0 else "")) for cu, c in zip(cum, cols)]
        half = len(u) // 2
        strips = []
        for j0, j1 in ((0, half), (half, len(u))):
            strips.append(strip(img[j0:j1], cr[j0:j1], [(j - j0, lab) for j, lab in marks if j0 <= j < j1]))
        path = os.path.join(a.out, f"cal2_step{i}.png")
        size = frame(strips, f"cal2, real CT, step {i}: ribbon along the trace front (labels = trace column; "
                             f"top strip columns 0 to about 18, bottom strip the rest)", [], path)
        info[os.path.basename(path)] = {"kind": "step ribbon", "step": i, "columns": cols.tolist(),
                                         "size_px": size, "um_per_px": round(DU * L.UM / UP, 3)}
        print(path, size)
    for c in [int(x) for x in a.cols.split(",") if x]:
        rows = [r for r in range(P.shape[0]) if kept[r, c]]
        u, cum, base, n, img = L.ribbon(ct, normal_at, P[rows, c], DU, 30 * 8, s, masked=True)
        cr = crossings(base, n)
        marks = [(int(np.argmin(np.abs(u - cu))), f"s{r}") for cu, r in zip(cum, rows)]
        ridge, extra = None, []
        if rf is not None and str(c) in rf["columns"]:
            # paths are [u_um, s_um] on the ridge script's grid (same DU, u = 0 at the seed); map by u value
            ridge = []
            for key, tag in (("path", "ridge"), ("path_f2", "f2"), ("path_f3", "f3")):
                for uu, x in (rf["columns"][str(c)].get(key) or []):
                    ridge.append((int(np.argmin(np.abs(u - uu / L.UM))), x, tag))
            extra = [("ridge", "F1 from seed"), ("f2", "F2 back from trace"), ("f3", "F3 back from w035")]
        st = strip(img, cr, marks, ridge)
        path = os.path.join(a.out, f"cal2_col{c}.png")
        size = frame([st], f"cal2, real CT, trace column {c}: ribbon along the arc through steps 0 to 8 "
                           f"(s0 = seed on w035)", extra, path)
        info[os.path.basename(path)] = {"kind": "column ribbon", "col": c, "rows": rows, "size_px": size,
                                         "um_per_px": round(DU * L.UM / UP, 3), "ridge_path_drawn": ridge is not None}
        print(path, size)
    fj = os.path.join(a.out, "figures.json")
    old = json.load(open(fj)) if os.path.exists(fj) else {}
    old.update(info)
    json.dump(old, open(fj, "w"), indent=1)


if __name__ == "__main__":
    main()
