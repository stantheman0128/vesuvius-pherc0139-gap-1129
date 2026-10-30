# Follow-up download sizing: count 1.129 um chunks (L0 and L1) the follow-up would need.
# Run as: python -I cov2167c.py [<data-dir>]
import sys, json, gzip
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from repo_paths import DATA
import numpy as np
import tifffile

DL = sys.argv[1] if len(sys.argv) > 1 else DATA
def A34(m):
    a = np.eye(4); a[:3, :] = np.array(m, float); return a
raw = open(f"{DL}/catalog-metadata.json", "rb").read()
raw = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
cat = json.loads(raw.decode("utf-8"))["samples"]["PHerc0139"]
vt = {(e["from_volume_id"], t["to_volume_id"]): A34(t["matrix"])
      for e in cat["sample"]["properties"]["volume_transforms"] for t in e["transforms"]}
T924, T2411 = vt[("20250728140407", "20260102150214")], vt[("20260102150214", "20260413113053")]
u = json.load(open(f"{DL}/umbilicus.json")); cp = sorted(u["control_points"], key=lambda c: c["z"])
uz, ux, uy = (np.array([c[k] for c in cp], float) for k in "zxy")
CH = 2097152  # bytes per chunk (128^3 uint8, compressor null)

def dense(w, up=4):
    g = [tifffile.imread(f"{DL}/seg/{w}-9.362/{c}.tif").astype(float) for c in "xyz"]
    v = (g[0] > 0) & (g[1] > 0) & (g[2] > 0)
    H, W = v.shape
    out = []
    for fy in np.linspace(0, 1, up, endpoint=False):
        for fx in np.linspace(0, 1, up, endpoint=False):
            ok = v[:-1, :-1] & v[1:, :-1] & v[:-1, 1:] & v[1:, 1:]
            c = [(1-fy)*(1-fx)*a[:-1, :-1] + fy*(1-fx)*a[1:, :-1] + (1-fy)*fx*a[:-1, 1:] + fy*fx*a[1:, 1:] for a in g]
            out.append(np.stack([c[0][ok], c[1][ok], c[2][ok]], 1))
    p = np.concatenate(out)
    p24 = (np.c_[p, np.ones(len(p))] @ T924.T)[:, :3]
    return p24, (np.c_[p24, np.ones(len(p24))] @ T2411.T)[:, :3]

allc0, allc1 = set(), set()
for w in ("w033", "w034", "w035", "w036", "w037"):
    p24, p11 = dense(w)
    a = np.degrees(np.arctan2(p24[:, 1] - np.interp(p24[:, 2], uz, uy), p24[:, 0] - np.interp(p24[:, 2], uz, ux))) % 360
    sel = ((a >= 327) | (a <= 85)) & (p11[:, 2] >= 16) & (p11[:, 2] <= 19296)
    q = p11[sel]
    c0 = {tuple(t) for t in np.floor(q / 128).astype(int)}
    c1 = {tuple(t) for t in np.floor(q / 256).astype(int)}
    allc0 |= c0; allc1 |= c1
    print(f"{w}: {sel.sum()} dense pts in covered gap envelope; L0 chunks {len(c0)}, L1 chunks {len(c1)}")

def dil(s):
    out = set()
    for x, y, z in s:
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                out.add((x+dx, y+dy, z))
    return out
for name, s in (("L0", allc0), ("L1", allc1)):
    d = dil(s)
    print(f"union w033-w037 {name}: {len(s)} chunks = {len(s)*CH/1e9:.1f} GB; dilated 1 chunk in x,y: {len(d)} = {len(d)*CH/1e9:.1f} GB")

# per CT place: box +-0.75 mm in x,y and +-0.25 mm in z around the place (as in the note's look)
places = {"E1": (16642, 14927, 14856), "E2": (18358, 12754, 17995), "E5": (17427, 13862, 13741),
          "E6": (17029, 14405, 14109), "ex w034 434,248": (17115, 14541, 12837), "ex w036 542,1016": (18213, 11594, 18432)}
tot = 0
for n, p in places.items():
    c = (np.r_[p, 1] @ T2411.T)[:3]
    h = np.array([750, 750, 250]) / 1.129
    lo, hi = np.floor((c - h) / 128).astype(int), np.floor((c + h) / 128).astype(int)
    k = int(np.prod(hi - lo + 1)); tot += k
    print(f"{n}: 1.129 L0 xyz ({c[0]:.0f},{c[1]:.0f},{c[2]:.0f}); chunks x{lo[0]}-{hi[0]} y{lo[1]}-{hi[1]} z{lo[2]}-{hi[2]}: {k} = {k*CH/1e9:.2f} GB")
print(f"six places total {tot} chunks = {tot*CH/1e9:.1f} GB")
