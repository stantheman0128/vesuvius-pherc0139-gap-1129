# Step 3a: resample the fetched 1.129 um data onto the radial frame of each place, and compute mesh traces.
# Output: $PHERC0139_OUT/pilot/<place>.npz with
#   slice   (tang, r) image at z0, mean of 5 planes (+-2 vox), r/tang in um on the 2.399 frame (1 um steps)
#   cut     (dz, r)  radial cut, mean over tangential +-TAV vox, 1 um steps
#   meshes: slice traces (r, tang) per mesh, cut traces (r, dz) per mesh
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pcommon import *
from scipy.ndimage import map_coordinates

plan = json.load(open(os.path.join(OUT, "plan.json")))
ZP0, ZP1 = plan["ZP"]
PLACES = {"E2": (359.805, 9380, 10250), "C240": (240.0, 7980, 8830)}
z0 = 17995.0
xu, yu = umb(z0)
o = np.array([xu, yu, z0])

def frame(a):
    d = np.array([np.cos(np.radians(a)), np.sin(np.radians(a)), 0.0])
    return d, np.array([-d[1], d[0], 0.0])

def p24(d, t, r_um, tang_um, dz_um):
    return o + np.multiply.outer(r_um / UM24, d) + np.multiply.outer(tang_um / UM24, t) + np.multiply.outer(dz_um / UM24, [0, 0, 1.0])

def assemble(keys, planes=None):
    ks = [tuple(k) for k in keys]
    z0c, y0c, x0c = (min(k[i] for k in ks) for i in range(3))
    z1c, y1c, x1c = (max(k[i] for k in ks) for i in range(3))
    nz = (planes[1] - planes[0] + 1) if planes else 128 * (z1c - z0c + 1)
    vol = np.zeros((nz, 128 * (y1c - y0c + 1), 128 * (x1c - x0c + 1)), np.uint8)
    have = np.zeros_like(vol, dtype=bool) if False else None
    for (z, y, x) in ks:
        fn = os.path.join(DLP, f"{z}_{y}_{x}")
        if os.path.exists(fn):
            blk = np.frombuffer(open(fn, "rb").read(), np.uint8).reshape(128, 128, 128)
            if planes:
                blk = blk[planes[0]:planes[1] + 1]
        else:
            fn = os.path.join(DLP, f"{z}_{y}_{x}_r{planes[0]*16384}-{(planes[1]+1)*16384-1}")
            blk = np.frombuffer(open(fn, "rb").read(), np.uint8).reshape(-1, 128, 128)
        zs = 0 if planes else 128 * (z - z0c)
        vol[zs:zs + blk.shape[0], 128 * (y - y0c):128 * (y - y0c + 1), 128 * (x - x0c):128 * (x - x0c + 1)] = blk
    org = np.array([128 * x0c, 128 * y0c, (128 * z0c + planes[0]) if planes else 128 * z0c], float)  # xyz
    return vol, org

def sample(vol, org, pts24):
    q = to11(pts24.reshape(-1, 3)) - org
    v = map_coordinates(vol.astype(np.float32), [q[:, 2], q[:, 1], q[:, 0]], order=1, mode="constant", cval=np.nan)
    return v.reshape(pts24.shape[:-1])

def plane_trace(g, zc):
    """Intersection of the triangulated grid with plane z = zc (2.399 frame): list of segment endpoints (N,2,3)."""
    a, b, c, e = g[:-1, :-1], g[1:, :-1], g[:-1, 1:], g[1:, 1:]
    segs = []
    for tri in ((a, b, c), (e, c, b)):
        P = np.stack(tri, -2)  # (...,3,3)
        ok = np.all(np.isfinite(P[..., 0]), -1)
        P = P[ok]
        s = P[..., 2] - zc
        cross = (s.min(1) < 0) & (s.max(1) >= 0)
        P, s = P[cross], s[cross]
        pts = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            m = (s[:, i] < 0) != (s[:, j] < 0)
            f = np.where(m, s[:, i] / np.where(m, s[:, i] - s[:, j], 1), np.nan)
            pts.append(P[:, i] + f[:, None] * (P[:, j] - P[:, i]))
        pts = np.stack(pts, 1)  # (N,3,3), two finite
        out = np.full((len(P), 2, 3), np.nan)
        for n in range(len(P)):
            fin = pts[n][np.isfinite(pts[n][:, 0])]
            if len(fin) >= 2: out[n] = fin[:2]
        segs.append(out)
    return np.concatenate(segs)

TAV = 3  # tangential half-width (1.129 vox ~ um) averaged in the radial cut
for name, (a, r0, r1) in PLACES.items():
    d, t = frame(a)
    cen = o + ((r0 + r1) / 2 / UM24) * d
    hw = np.array([800, 800, 260]) / UM24
    MESH = {w: mesh24_fine(w, (cen - hw, cen + hw)) for w in ("w033", "w034", "w035", "w036", "w037")}
    for w in MESH:
        h9 = [round(x[0] * UM24) for x in ray_hits(mesh24(w), o, d) if r0 - 50 <= x[0] * UM24 <= r1 + 50]
        hf = [round(x[0] * UM24, 1) for x in ray_hits(MESH[w], o, d) if r0 - 50 <= x[0] * UM24 <= r1 + 50] if MESH[w].size else []
        print(f"  {name} {w}: crossings at z0 (um from umbilicus): 2.399 copy {hf}; 9.362 copy {h9}; grid crop {MESH[w].shape[:2]}")
    full = plan[name]["full"]; sl = plan[name]["slice"]
    # slice volume: planes ZP0..ZP1 of all slice chunks
    svol, sorg = assemble(full_zc := [k for k in full if k[0] == sl[0][0]] + sl, planes=(ZP0, ZP1))
    rmid = (r0 + r1) / 2
    rr = np.arange(rmid - 650, rmid + 651, 1.0)
    tt = np.arange(-650, 651, 1.0)
    R, Tg = np.meshgrid(rr, tt)
    sl_img = np.nanmean([sample(svol, sorg, p24(d, t, R, Tg, np.full_like(R, dz))) for dz in np.arange(-2, 3) * UM11], 0)
    # radial cut
    cvol, corg = assemble(full)
    rc = np.arange(r0, r1 + 1, 1.0)
    zz = np.arange(-200, 201, 1.0)
    Rc, Zc = np.meshgrid(rc, zz)
    cut = np.nanmean([sample(cvol, corg, p24(d, t, Rc, np.full_like(Rc, k * UM11), Zc)) for k in range(-TAV, TAV + 1)], 0)
    cut1 = sample(cvol, corg, p24(d, t, Rc, np.zeros_like(Rc), Zc))
    # mesh traces in the slice and in the cut
    traces, cuthits = {}, {}
    for w, g in MESH.items():
        segs = plane_trace(g, z0)
        rel = segs - o
        rs, ts = (rel @ d) * UM24, (rel @ t) * UM24
        keep = np.all((np.abs(rs - rmid) < 700) & (np.abs(ts) < 700), 1)
        traces[w] = np.stack([rs[keep], ts[keep]], -1)
        h = []
        for dz in np.arange(-200, 201, 10.0):
            oo = o + np.array([0, 0, dz / UM24])
            for tt_, _, _ in ray_hits(g, oo, d):
                if r0 - 50 <= tt_ * UM24 <= r1 + 50:
                    h.append((tt_ * UM24, dz))
        cuthits[w] = np.array(h).reshape(-1, 2)
    np.savez_compressed(os.path.join(OUT, f"{name}.npz"), slice=sl_img, rr=rr, tt=tt, cut=cut, cut1=cut1, rc=rc, zz=zz,
                        **{f"tr_{w}": v for w, v in traces.items()}, **{f"ch_{w}": v for w, v in cuthits.items()})
    print(name, "slice", sl_img.shape, "nan%", np.isnan(sl_img).mean().round(3), "cut", cut.shape, "nan%", np.isnan(cut).mean().round(3),
          {w: len(v) for w, v in cuthits.items()})
