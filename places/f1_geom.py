# Step 1: radial line per place and all mesh crossings on it (9.362 um copies over the whole ray). No downloads.
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fcommon import *

M9 = {w: mesh24(w) for w in MESHES}
geo = {}
for name in sys.argv[1:] or PLACES:
    z0, a, o, d, t = place_frame(name)
    xy = PLACES[name][2]
    rp = np.hypot(xy[0] - o[0], xy[1] - o[1]) * UM24 if xy else None
    print(f"== {name}: z {z0:.0f}, angle {a:.3f} deg, umbilicus ({o[0]:.1f}, {o[1]:.1f})" + (f", place r {rp:.0f} um" if rp else ""))
    hits = {}
    for w, g in M9.items():
        hs = [round(h[0] * UM24, 1) for h in ray_hits(g, o, d)]
        hits[w] = hs
        print(f"   {w}: {hs}")
    geo[name] = {"z0": z0, "angle": a, "r_place": rp, "hits9": hits}
json.dump(geo, open(os.path.join(OUT5, "geo.json"), "w"), indent=1)
