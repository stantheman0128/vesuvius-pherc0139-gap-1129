# Cross-check: catalog 1.129->9.362 chain vs the independent #1381 transform (axis_order zyx).
import sys, json, gzip
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from repo_paths import DATA
import numpy as np
DL = sys.argv[1] if len(sys.argv) > 1 else DATA
def A34(m):
    a = np.eye(4); a[:3, :] = np.array(m, float); return a
raw = open(f"{DL}/catalog-metadata.json", "rb").read()
raw = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
cat = json.loads(raw.decode("utf-8"))["samples"]["PHerc0139"]
vt = {(e["from_volume_id"], t["to_volume_id"]): A34(t["matrix"]) for e in cat["sample"]["properties"]["volume_transforms"] for t in e["transforms"]}
print("catalog has 20260413113053->20250728140407:", ("20260413113053", "20250728140407") in vt)
C11_24 = vt[("20260413113053", "20260102150214")]; C24_9 = np.linalg.inv(vt[("20250728140407", "20260102150214")])
C11_9 = C24_9 @ C11_24
C9_24 = vt[("20250728140407", "20260102150214")]
R = json.load(open(f"{DL}/pr1381-20260413113053-to-20250728140407.json"))
Mz = A34(R["transformation_matrix"])  # zyx -> zyx
P = np.array([[0,0,0,0,0],[0,0,0,0,0]])
def r1381(xyz):  # xyz in 1.129 -> xyz in 9.362
    zyx = xyz[:, ::-1]
    out = (np.c_[zyx, np.ones(len(zyx))] @ Mz.T)[:, :3]
    return out[:, ::-1]
# sample points: grid over the 1.129 volume
sh = [19297, 35971, 32318]
g = np.array([[x, y, z] for x in np.linspace(0, sh[2], 5) for y in np.linspace(0, sh[1], 5) for z in np.linspace(16, 19296, 5)], float)
a = (np.c_[g, np.ones(len(g))] @ C11_9.T)[:, :3]
b = r1381(g)
d = (a - b) * 9.362
print("catalog chain vs #1381, 1.129 grid points mapped to 9.362 frame, diff in um: median %.0f, max %.0f; z-diff median %.0f max %.0f" % (
    np.median(np.linalg.norm(d, axis=1)), np.linalg.norm(d, axis=1).max(), np.median(np.abs(d[:, 2])), np.abs(d[:, 2]).max()))
# 1.129 z-extent mapped into 2.399 via #1381 then catalog 9.362->2.399
for zz in (16, 19296):
    gg = np.array([[x, y, zz] for x in np.linspace(0, sh[2], 9) for y in np.linspace(0, sh[1], 9)], float)
    z24 = (np.c_[r1381(gg), np.ones(len(gg))] @ C9_24.T)[:, 2]
    z24c = (np.c_[gg, np.ones(len(gg))] @ C11_24.T)[:, 2]
    print(f"z_1.129={zz}: 2.399 z via #1381+catalog9->2.4 {z24.min():.0f}..{z24.max():.0f}; via catalog direct {z24c.min():.0f}..{z24c.max():.0f}")
# decomposition of catalog chain 9.362->1.129
M = np.linalg.inv(C11_9)[:3, :3]
s = np.cbrt(abs(np.linalg.det(M)))
print("catalog 9.362->1.129 isotropic scale %.5f, nominal 9.362/1.129 = %.5f, ratio %.5f" % (s, 9.362/1.129, s/(9.362/1.129)))
print("rotation about z (deg): %.2f" % np.degrees(np.arctan2(M[1,0], M[0,0])))
