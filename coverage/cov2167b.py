# Second check: (a) margin of covered gap vertices to the nearest empty L3 chunk, (b) voxel content at mesh
# positions inside two downloaded level-5 chunks.  Run as: python -I cov2167b.py [<data-dir> [<keys-L3.txt>]]
import sys, json, gzip
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from repo_paths import DATA
import numpy as np
import tifffile

DL = sys.argv[1] if len(sys.argv) > 1 else DATA
KEYS = sys.argv[2] if len(sys.argv) > 2 else f"{DL}/keys-L3.txt"
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

occ = np.zeros((19, 36, 32), bool)
for line in open(KEYS):
    p = line.split()[0].split(".zarr/")[1].split("/")
    if len(p) == 4 and p[0] == "3":
        occ[int(p[1]), int(p[2]), int(p[3])] = True
print("L3 occupancy per z-chunk (of 36x32):", occ.sum(axis=(1, 2)).tolist())

def load(w):
    g = [tifffile.imread(f"{DL}/seg/{w}-9.362/{c}.tif").astype(float) for c in "xyz"]
    v = (g[0] > 0) & (g[1] > 0) & (g[2] > 0)
    p = np.stack([g[0][v], g[1][v], g[2][v]], 1)
    p24 = (np.c_[p, np.ones(len(p))] @ T924.T)[:, :3]
    return p24, (np.c_[p24, np.ones(len(p24))] @ T2411.T)[:, :3]

# (a) margin: Chebyshev distance (in L3 chunks, xy only, same z chunk) to the nearest empty chunk
from itertools import product
def margin(ci):
    zc, yc, xc = ci
    for d in range(0, 40):
        for dy, dx in product(range(-d, d + 1), repeat=2):
            if max(abs(dy), abs(dx)) != d: continue
            y, x = yc + dy, xc + dx
            if y < 0 or x < 0 or y >= 36 or x >= 32 or not occ[zc, y, x]:
                return d
    return 99
for w in ("w034", "w036"):
    p24, p11 = load(w)
    a = np.degrees(np.arctan2(p24[:, 1] - np.interp(p24[:, 2], uz, uy), p24[:, 0] - np.interp(p24[:, 2], uz, ux))) % 360
    sel = ((a >= 327) | (a <= 85)) & (p24[:, 2] >= 9000) & (p24[:, 2] <= 27800) & (p11[:, 2] >= 16) & (p11[:, 2] <= 19296)
    ci = np.floor(p11[sel] / 1024).astype(int)[:, ::-1]
    ms = np.array([margin(tuple(c)) for c in {tuple(c) for c in ci}])
    print(f"{w}: covered gap vertices sit in {len(ms)} distinct L3 chunks; min xy distance to an empty chunk = {ms.min()} chunks"
          f" ({ms.min()*1024*1.129/1000:.1f} mm), median {np.median(ms):.0f}")

# (b) voxel content at mesh positions inside the two L5 chunks
for name in ("1_5_5", "3_4_6"):
    zc, yc, xc = map(int, name.split("_"))
    blk = np.frombuffer(open(f"{DL}/L5/{name}", "rb").read(), np.uint8).reshape(128, 128, 128)
    print(f"\nL5 chunk {name}: nonzero voxels {100*(blk>0).mean():.1f}%, mean {blk.mean():.1f}")
    for w in ("w033", "w034", "w036", "w037"):
        p24, p11 = load(w)
        q = np.floor(p11 / 32).astype(int)  # L5 voxel (x,y,z)
        lo = np.array([xc, yc, zc]) * 128
        m = np.all((q >= lo) & (q < lo + 128), axis=1)
        if not m.any():
            print(f"  {w}: no vertices in this chunk"); continue
        loc = q[m] - lo
        vals = blk[loc[:, 2], loc[:, 1], loc[:, 0]]
        print(f"  {w}: {m.sum()} vertices in chunk (2.399 z {p24[m,2].min():.0f}..{p24[m,2].max():.0f});"
              f" voxel value >0 at {100*(vals>0).mean():.1f}%, median value {np.median(vals):.0f}")
