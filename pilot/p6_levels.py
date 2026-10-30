# Step 4b: two more z slices at E2 (dz -150 and +150 um), to see in xy whether the bands that merge in the radial cut
# are separate sheets touching or one sheet. Byte-range GETs of planes 47..63 of z chunks 108 and 110 (11x11 chunks).
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pcommon import *
from concurrent.futures import ThreadPoolExecutor
from scipy.ndimage import map_coordinates
import p3_sample as P3  # reuses assemble/sample/plane_trace/frame (module runs its main loop on import: cheap, cached data)

plan = json.load(open(os.path.join(OUT, "plan.json")))
ZP0, ZP1 = plan["ZP"]
full = [tuple(k) for k in plan["E2"]["full"]]
sl109 = [tuple(k) for k in plan["E2"]["slice"]]
yx = sorted({(y, x) for _, y, x in sl109} | {(y, x) for _, y, x in full})
out = {}
for zc, dz in ((108, -150.0), (110, 150.0)):
    jobs = [(f"{zc}/{y}/{x}", (ZP0 * 16384, (ZP1 + 1) * 16384 - 1)) for y, x in yx if (zc, y, x) not in set(full)]
    assert dl_bytes() + len(jobs) * (ZP1 - ZP0 + 1) * 16384 < 500_000_000
    with ThreadPoolExecutor(8) as ex:
        res = list(ex.map(lambda j: fetch(*j), jobs))
    print(f"z chunk {zc}: {len(jobs)} range GETs, 404 {sum(r is None for r in res)}")
    keys = [(zc, y, x) for y, x in yx]
    svol, sorg = P3.assemble(keys, planes=(ZP0, ZP1))
    a, r0, r1 = P3.PLACES["E2"]
    d, t = P3.frame(a)
    rmid = (r0 + r1) / 2
    rr = np.arange(rmid - 650, rmid + 651, 1.0); tt = np.arange(-650, 651, 1.0)
    R, Tg = np.meshgrid(rr, tt)
    img = np.nanmean([P3.sample(svol, sorg, P3.p24(d, t, R, Tg, np.full_like(R, dz + k))) for k in np.arange(-2, 3) * UM11], 0)
    zq = to11(P3.p24(d, t, np.array([rmid]), np.array([0.0]), np.array([dz])))[0, 2]
    print(f"  dz {dz:+.0f}: z11 at line centre {zq:.1f} (planes held {128*zc+ZP0}..{128*zc+ZP1}); nan {np.isnan(img).mean():.3f}")
    cen = P3.o + (rmid / UM24) * d
    hw = np.array([800, 800, 260]) / UM24
    tr = {}
    for w in ("w033", "w034", "w036", "w037"):
        g = mesh24_fine(w, (cen - hw, cen + hw))
        segs = P3.plane_trace(g, P3.z0 + dz / UM24)
        rel = segs - P3.o
        rs, ts = (rel @ d) * UM24, (rel @ t) * UM24
        keep = np.all((np.abs(rs - rmid) < 700) & (np.abs(ts) < 700), 1)
        tr[w] = np.stack([rs[keep], ts[keep]], -1)
    out[dz] = (img, tr)
D = np.load(os.path.join(OUT, "E2.npz"))
np.savez_compressed(os.path.join(OUT, "E2_levels.npz"), rr=D["rr"], tt=D["tt"],
                    **{f"img_{int(k)}": v[0] for k, v in out.items()},
                    **{f"tr_{int(k)}_{w}": s for k, v in out.items() for w, s in v[1].items()})
print("dir total", dl_bytes() / 1e6, "MB")
