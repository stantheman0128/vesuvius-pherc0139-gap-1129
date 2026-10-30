# Step 1: radial line through E2 and mesh crossings; candidate control angles at the same z.
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pcommon import *

E2 = np.array([18358., 12754., 17995.])
z0 = E2[2]
xu, yu = umb(z0)
ang = np.degrees(np.arctan2(E2[1] - yu, E2[0] - xu)) % 360
rE2 = np.hypot(E2[0] - xu, E2[1] - yu)
print(f"umbilicus at z {z0:.0f}: ({xu:.2f}, {yu:.2f}); E2 angle {ang:.3f} deg; E2 radius {rE2:.1f} vox = {rE2*UM24:.0f} um")
print("umbilicus in 1.129 frame:", np.round(to11([xu, yu, z0])[0], 1), " E2 in 1.129:", np.round(to11(E2)[0], 1))

meshes = {w: mesh24(w) for w in ("w033", "w034", "w035", "w036", "w037")}
o = np.array([xu, yu, z0])

def crossings(a_deg):
    d = np.array([np.cos(np.radians(a_deg)), np.sin(np.radians(a_deg)), 0.0])
    res = {}
    for w, g in meshes.items():
        res[w] = [(t, r, c) for t, r, c in ray_hits(g, o, d)]
    return d, res

d, res = crossings(ang)
print(f"\nE2 radial line: origin = umbilicus at z {z0:.0f}, direction {ang:.3f} deg; crossings (r in 2.399 vox, um; 1.129 xyz; 9.362 grid row,col):")
for w, hs in res.items():
    for t, r, c in hs:
        p = o + t * d
        print(f"  {w}: r {t:8.1f} vox = {t*UM24:7.0f} um; 1.129 xyz {np.round(to11(p)[0], 1)}; grid ({r},{c})")

print("\nControl candidates at the same z (all five meshes, r in um):")
for a in list(range(90, 330, 10)):
    d, res = crossings(a)
    row = []
    for w in ("w033", "w034", "w035", "w036", "w037"):
        row.append(",".join(f"{t*UM24:.0f}" for t, _, _ in res[w]) or "-")
    p = o + (np.mean([h[0] for h in res["w034"]]) if res["w034"] else 0) * d
    q = to11(p)[0]
    print(f"  {a:3d} deg: " + " | ".join(f"{w} {s}" for w, s in zip(meshes, row)) + f"   w034 at 1.129 xyz {np.round(q)}")
