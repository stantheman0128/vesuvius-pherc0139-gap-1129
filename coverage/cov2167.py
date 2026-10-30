# Coverage check: does the 1.129 um PHerc0139 scan (20260413113053) cover the #2167 gap sector?
# Run as: python -I cov2167.py [<data-dir> [<keys-L3.txt>]]
import sys, json, gzip, math
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from repo_paths import DATA
import numpy as np
import tifffile

DL = sys.argv[1] if len(sys.argv) > 1 else DATA
KEYS = sys.argv[2] if len(sys.argv) > 2 else f"{DL}/keys-L3.txt"

def A34(m):
    a = np.eye(4); a[:3, :] = np.array(m, dtype=float); return a

raw = open(f"{DL}/catalog-metadata.json", "rb").read()
if raw[:2] == b"\x1f\x8b":
    raw = gzip.decompress(raw)
cat = json.loads(raw.decode("utf-8"))["samples"]["PHerc0139"]
vt = {}
for e in cat["sample"]["properties"]["volume_transforms"]:
    for t in e["transforms"]:
        vt[(e["from_volume_id"], t["to_volume_id"])] = A34(t["matrix"])
T_9_to_24 = vt[("20250728140407", "20260102150214")]
T_24_to_11 = vt[("20260102150214", "20260413113053")]
T_11_to_24 = vt[("20260413113053", "20260102150214")]
tj = json.load(open(f"{DL}/20260413113053-1.129um-0.2m-59keV-masked/transform.json"))
assert np.allclose(A34(tj["transformation_matrix"]), T_11_to_24), "transform.json != catalog"
print("catalog 1.129->2.399 matrix == volume transform.json: True")
print("max |inv(T_11_to_24) - catalog T_24_to_11| =", np.abs(np.linalg.inv(T_11_to_24) - T_24_to_11).max())
print("landmark residuals (catalog matrix, 2.399 vox):",
      np.round(np.linalg.norm((np.c_[np.array(tj["moving_landmarks"]), np.ones(12)] @ T_11_to_24.T)[:, :3]
                              - np.array(tj["fixed_landmarks"]), axis=1), 1))

SHAPE11 = json.load(open(f"{DL}/20260413113053-1.129um-0.2m-59keV-masked/0_.zarray"))["shape"]  # z,y,x
ZC0, ZC1 = 16, 19296  # z_crop_start/end from metadata.json zarr_export
print("1.129 shape (z,y,x):", SHAPE11)

# extent of the 1.129 box in 2.399 frame (8 corners)
cs = np.array([[x, y, z, 1] for x in (0, SHAPE11[2]) for y in (0, SHAPE11[1]) for z in (0, SHAPE11[0])], float)
cf = cs @ T_11_to_24.T
print("1.129 volume corners in 2.399 frame: x %.0f..%.0f  y %.0f..%.0f  z %.0f..%.0f" % (
    cf[:, 0].min(), cf[:, 0].max(), cf[:, 1].min(), cf[:, 1].max(), cf[:, 2].min(), cf[:, 2].max()))
for zz in (0, ZC0, ZC1, SHAPE11[0]):
    zs = (np.array([[x, y, zz, 1] for x in (0, SHAPE11[2]) for y in (0, SHAPE11[1])], float) @ T_11_to_24.T)[:, 2]
    print(f"  z_1.129={zz}: z_2.399 {zs.min():.0f}..{zs.max():.0f}")

# L3 chunk occupancy
present = set()
for line in open(KEYS):
    k = line.split()[0].split(".zarr/")[1]
    p = k.split("/")
    if len(p) == 4 and p[0] == "3":
        present.add((int(p[1]), int(p[2]), int(p[3])))
print("L3 chunks present:", len(present))

def in_scan(p11):
    """p11: (N,3) xyz in 1.129 L0 voxels -> (inbox, inchunk)"""
    x, y, z = p11[:, 0], p11[:, 1], p11[:, 2]
    inbox = (x >= 0) & (x < SHAPE11[2]) & (y >= 0) & (y < SHAPE11[1]) & (z >= ZC0) & (z <= ZC1)
    ci = np.floor(p11 / 1024.0).astype(int)  # L3 chunk = 128 * 8 L0 voxels
    inchunk = np.array([inbox[i] and (ci[i, 2], ci[i, 1], ci[i, 0]) in present for i in range(len(p11))])
    return inbox, inchunk

# umbilicus (2.399 frame)
u = json.load(open(f"{DL}/umbilicus.json"))
cp = sorted(u["control_points"], key=lambda c: c["z"])
uz = np.array([c["z"] for c in cp], float); ux = np.array([c["x"] for c in cp], float); uy = np.array([c["y"] for c in cp], float)
def umb(z):
    return np.interp(z, uz, ux), np.interp(z, uz, uy)
for zz in (9000, 11400, 16000, 20490, 27800):
    print("umbilicus at z %d: (%.0f, %.0f)" % ((zz,) + tuple(umb(zz))))

def angle(p):
    xu, yu = umb(p[:, 2])
    return np.degrees(np.arctan2(p[:, 1] - yu, p[:, 0] - xu)) % 360.0, np.hypot(p[:, 0] - xu, p[:, 1] - yu)

# named points from ChaseHendrick note (x,y,z in 2.399 frame)
pts = {
    "ex w034 434,248 gap": (17115, 14541, 12837), "ex w036 542,1016 gap": (18213, 11594, 18432),
    "ex w036 854,1110 gap(#2167 body)": (17834, 12866, 24543), "ex w034 472,322 skip": (16038, 15479, 13682),
    "ex w036 530,1102 skip": (18188, 13254, 18335), "E1": (16642, 14927, 14856), "E2": (18358, 12754, 17995),
    "E4": (17611, 13846, 10437), "E5": (17427, 13862, 13741), "E6": (17029, 14405, 14109),
}
P = np.array(list(pts.values()), float)
P11 = (np.c_[P, np.ones(len(P))] @ T_24_to_11.T)[:, :3]
ib, ic = in_scan(P11)
ang, rad = angle(P)
print("\nNamed places (2.399 xyz -> 1.129 xyz), angle recomputed, in box, in nonempty L3 chunk:")
for (n, p), q, a, r, b, c in zip(pts.items(), P11, ang, rad, ib, ic):
    print(f"  {n:34s} {p} -> ({q[0]:.0f},{q[1]:.0f},{q[2]:.0f})  ang {a:5.1f}  r {r*2.399/1000:.2f} mm  box={b} chunk={c}")
# G1-G4, C1-C2 are z-only places; report z in/out
for n, z in (("G1", 11000), ("G2", 19500), ("G3", 23000), ("G4", 26800), ("C1", 13000), ("C2", 22000),
             ("unwrap", 14856), ("unwrap", 19500)):
    print(f"  {n:6s} z {z}: inside 1.129 z band? {11395 <= z <= 20490}")

# meshes
def load9(w):
    g = [tifffile.imread(f"{DL}/seg/{w}-9.362/{c}.tif").astype(np.float64) for c in "xyz"]
    valid = (g[0] > 0) & (g[1] > 0) & (g[2] > 0) & np.isfinite(g[0])
    p = np.stack([g[0][valid], g[1][valid], g[2][valid]], 1)
    return (np.c_[p, np.ones(len(p))] @ T_9_to_24.T)[:, :3]

print("\nPer-mesh check (9.362 um copies mapped to 20260102150214 with the catalog matrix):")
rows = []
for w in ("w033", "w034", "w035", "w036", "w037"):
    p = load9(w)
    mb = json.load(open(f"{DL}/seg/{w}-on-20260102150214-2.399um-meta.json"))["bbox"]
    print(f" {w}: {len(p)} vertices; mapped bbox x {p[:,0].min():.0f}..{p[:,0].max():.0f} y {p[:,1].min():.0f}..{p[:,1].max():.0f} z {p[:,2].min():.0f}..{p[:,2].max():.0f}"
          f" | 2.399 meta bbox x {mb[0][0]:.0f}..{mb[1][0]:.0f} y {mb[0][1]:.0f}..{mb[1][1]:.0f} z {mb[0][2]:.0f}..{mb[1][2]:.0f}")
    a, r = angle(p)
    p11 = (np.c_[p, np.ones(len(p))] @ T_24_to_11.T)[:, :3]
    ibox, ichk = in_scan(p11)
    for label, sel_ang in (("envelope 327-85 deg", (a >= 327) | (a <= 85)), ("core 5-71 deg", (a >= 5) & (a <= 71))):
        sel = sel_ang & (p[:, 2] >= 9000) & (p[:, 2] <= 27800)
        n = sel.sum()
        if n == 0:
            print(f"   {label}: no vertices"); continue
        inb = (sel & ibox).sum(); inc = (sel & ichk).sum()
        zin = p[sel & ichk, 2]
        ain = a[sel & ichk]
        # angle span among covered vertices (wrap at 327)
        aw = np.where(ain >= 300, ain - 360, ain)
        if inc == 0:
            print(f"   {label}: {n} verts in z9000-27800 (z {p[sel,2].min():.0f}..{p[sel,2].max():.0f}); none in 1.129 scan"); continue
        print(f"   {label}: {n} verts in z9000-27800; in 1.129 box {inb} ({100*inb/n:.1f}%), in nonempty L3 chunk {inc} ({100*inc/n:.1f}%);"
              f" covered z {zin.min():.0f}..{zin.max():.0f}, covered angle {aw.min():.1f}..{aw.max():.1f}; radius {r[sel].min()*2.399/1000:.2f}-{r[sel].max()*2.399/1000:.2f} mm")
        # in-z-band but not covered?
        band = sel & (p[:, 2] >= 11420) & (p[:, 2] <= 20460)
        nb = band.sum(); nbc = (band & ichk).sum()
        print(f"      within z 11420-20460: {nb} verts, covered {nbc} ({100*nbc/max(nb,1):.2f}%)")
        if label.startswith("envelope"):
            hist = []
            for z0 in range(9000, 27800, 2000):
                s2 = sel & (p[:, 2] >= z0) & (p[:, 2] < z0 + 2000)
                if s2.sum():
                    hist.append(f"{z0}-{z0+2000}:{100*(s2&ichk).sum()/s2.sum():.0f}%")
            print("      by z band:", " ".join(hist))
