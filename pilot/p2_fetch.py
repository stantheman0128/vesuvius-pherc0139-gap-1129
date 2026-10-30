# Step 2: fetch level-0 chunks of the 1.129 um scan for E2 and the control.
#  (a) z slice: planes ZP0..ZP1 of z chunk 109 (byte-range GETs, chunks are uncompressed C-order z,y,x), 11x11 chunks
#  (b) radial cut: full chunks covering a +-SLAB voxel slab along the radial line, z chunks 108-110
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pcommon import *
from concurrent.futures import ThreadPoolExecutor

CAP = 500_000_000
ZP0, ZP1 = 47, 63           # planes of z chunk 109 for the slice (z11 13999..14015)
SLAB = 16                   # tangential half-width (1.129 vox) of the radial slab
PLACES = {  # name: (angle deg, r range um on the 2.399 frame), z = 17995
    "E2": (359.805, 9380, 10250),
    "C240": (240.0, 7980, 8830),
}
z0 = 17995.0
xu, yu = umb(z0)
o = np.array([xu, yu, z0])

def line11(a, r0, r1, off=0.0, dz=0.0, n=2000):
    d = np.array([np.cos(np.radians(a)), np.sin(np.radians(a)), 0.0])
    t = np.array([-d[1], d[0], 0.0])
    r = np.linspace(r0, r1, n) / UM24
    p = o + r[:, None] * d + (off / UM24) * t + np.array([0, 0, dz / UM24])
    return to11(p)

plan = {}
for name, (a, r0, r1) in PLACES.items():
    pts = np.concatenate([line11(a, r0, r1, off) for off in np.linspace(-SLAB * UM11, SLAB * UM11, 9)])
    cxy = {(int(y // 128), int(x // 128)) for x, y, _ in pts}
    zc = int(np.median(pts[:, 2]) // 128)
    full = [(z, y, x) for z in (zc - 1, zc, zc + 1) for (y, x) in sorted(cxy)]
    mid = line11(a, (r0 + r1) / 2, (r0 + r1) / 2, n=1)[0]
    cy, cx = int(mid[1] // 128), int(mid[0] // 128)
    sl = [(zc, y, x) for y in range(cy - 5, cy + 6) for x in range(cx - 5, cx + 6) if (y, x) not in cxy]
    plan[name] = (full, sl)
    print(f"{name}: line z11 {pts[:,2].min():.1f}..{pts[:,2].max():.1f} (z chunk {zc}, local {np.median(pts[:,2])-128*zc:.1f});"
          f" radial-cut full chunks {len(full)}; slice range-GETs {len(sl)} (+{len(cxy)} slice chunks inside the full set)")

est = sum(len(f) * 2097152 + len(s) * (ZP1 - ZP0 + 1) * 16384 for f, s in plan.values())
print(f"planned download {est/1e6:.1f} MB (cap {CAP/1e6:.0f} MB); already in dir {dl_bytes()/1e6:.1f} MB")
assert est + dl_bytes() < CAP

jobs = []
for name, (full, sl) in plan.items():
    jobs += [(f"{z}/{y}/{x}", None) for z, y, x in full]
    jobs += [(f"{z}/{y}/{x}", (ZP0 * 16384, (ZP1 + 1) * 16384 - 1)) for z, y, x in sl]
with ThreadPoolExecutor(8) as ex:
    res = list(ex.map(lambda j: fetch(*j), jobs))
miss = [j[0] for j, r in zip(jobs, res) if r is None]
bad = [j[0] for j, r in zip(jobs, res) if r is not None and len(r) != (2097152 if j[1] is None else (ZP1 - ZP0 + 1) * 16384)]
print(f"fetched {len(jobs)} requests; 404: {miss}; wrong size: {bad}; bytes this run {NBYTES[0]/1e6:.1f} MB; dir total {dl_bytes()/1e6:.1f} MB")
json.dump({k: {"full": v[0], "slice": v[1]} for k, v in plan.items()} | {"ZP": [ZP0, ZP1]},
          open(os.path.join(OUT, "plan.json"), "w"))
