# Shared helpers for the #2167 E2 pilot on the 1.129 um scan. Imported by the pilot and places scripts (run with python -I).
import sys, json, gzip, os, struct, urllib.request
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from repo_paths import DATA, OUT_PILOT, mesh_url
import numpy as np
import tifffile

DL = DATA                                   # catalog, umbilicus, 9.362 um meshes (fetch_inputs.py)
DLP = os.path.join(DATA, "chunks-pilot")    # this pilot's CT chunk downloads
MESHDIR = os.path.join(DATA, "mesh24")      # row bands of the 2.399 um mesh tifs
OUT = OUT_PILOT
os.makedirs(DLP, exist_ok=True)
BASE = "https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc0139/volumes/20260413113053-1.129um-0.2m-59keV-masked.zarr/0/"
UM11, UM24 = 1.129, 2.399

def A34(m):
    a = np.eye(4); a[:3, :] = np.array(m, float); return a

raw = open(f"{DL}/catalog-metadata.json", "rb").read()
raw = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
_cat = json.loads(raw.decode("utf-8"))["samples"]["PHerc0139"]
VT = {(e["from_volume_id"], t["to_volume_id"]): A34(t["matrix"])
      for e in _cat["sample"]["properties"]["volume_transforms"] for t in e["transforms"]}
T924 = VT[("20250728140407", "20260102150214")]
T2411 = VT[("20260102150214", "20260413113053")]
T1124 = VT[("20260413113053", "20260102150214")]

_u = json.load(open(f"{DL}/umbilicus.json")); _cp = sorted(_u["control_points"], key=lambda c: c["z"])
UZ, UX, UY = (np.array([c[k] for c in _cp], float) for k in "zxy")

def umb(z):
    return np.interp(z, UZ, UX), np.interp(z, UZ, UY)

def to11(p24):
    p24 = np.atleast_2d(p24)
    return (np.c_[p24, np.ones(len(p24))] @ T2411.T)[:, :3]

def to24(p11):
    p11 = np.atleast_2d(p11)
    return (np.c_[p11, np.ones(len(p11))] @ T1124.T)[:, :3]

def mesh24(w):
    """9.362 um tifxyz grid of mesh w mapped to 20260102150214 voxels: (H,W,3) array, NaN where invalid."""
    g = [tifffile.imread(f"{DL}/seg/{w}-9.362/{c}.tif").astype(np.float64) for c in "xyz"]
    v = (g[0] > 0) & (g[1] > 0) & (g[2] > 0) & np.isfinite(g[0])
    p = np.stack(g, -1)
    q = (np.c_[p.reshape(-1, 3), np.ones(p.size // 3)] @ T924.T)[:, :3].reshape(p.shape)
    q[~v] = np.nan
    return q

def _get(url, rng):
    req = urllib.request.Request(url, headers={"Range": f"bytes={rng[0]}-{rng[1]}"})
    data = urllib.request.urlopen(req, timeout=120).read()
    assert len(data) == rng[1] - rng[0] + 1, (url, rng, len(data))
    return data

def _tif_layout(w, c):
    """(data offset, height, width) of the 2.399 um tif c of mesh w, read from its header and IFD by range GETs.
    The published tifxyz tifs are little-endian, uncompressed float32 with a single strip."""
    fn = os.path.join(MESHDIR, w, f"{c}.layout.json")
    if os.path.exists(fn):
        return tuple(json.load(open(fn)))
    url = f"{mesh_url(w, '2.399')}/{c}.tif"
    head = _get(url, (0, 7))
    assert head[:4] == b"II*\x00", head
    ifd = struct.unpack("<I", head[4:])[0]
    n = struct.unpack("<H", _get(url, (ifd, ifd + 1)))[0]
    raw = _get(url, (ifd + 2, ifd + 1 + 12 * n))
    tags = {}
    for i in range(n):
        tag, typ, cnt = struct.unpack("<HHI", raw[12 * i:12 * i + 8])
        tags[tag] = (cnt, struct.unpack("<H" if typ == 3 else "<I", raw[12 * i + 8:12 * i + (10 if typ == 3 else 12)])[0])
    W, H = tags[256][1], tags[257][1]
    assert tags[258][1] == 32 and tags[259][1] == 1 and tags[339][1] == 3 and tags[273][0] == 1, tags
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    json.dump([tags[273][1], H, W], open(fn, "w"))
    return tags[273][1], H, W

def _tif_rows(w, c, r0, r1):
    """rows r0..r1 (inclusive) of the 2.399 um tif c of mesh w as float64 (r1 - r0 + 1, W). Cached in MESHDIR."""
    off, H, W = _tif_layout(w, c)
    d = os.path.join(MESHDIR, w)
    for f in os.listdir(d):
        m = f.startswith(f"{c}_rows") and f.endswith(".f32") and tuple(map(int, f[len(c) + 5:-4].split("-")))
        if m and m[0] <= r0 and m[1] >= r1:
            a = np.fromfile(os.path.join(d, f), "<f4").reshape(-1, W)
            return a[r0 - m[0]:r1 - m[0] + 1].astype(np.float64)
    data = _get(f"{mesh_url(w, '2.399')}/{c}.tif", (off + r0 * W * 4, off + (r1 + 1) * W * 4 - 1))
    open(os.path.join(d, f"{c}_rows{r0}-{r1}.f32"), "wb").write(data)
    return np.frombuffer(data, "<f4").reshape(-1, W).astype(np.float64)

def mesh24_fine(w, box=None):
    """Published 2.399 um tifxyz of mesh w (already in 20260102150214 voxels, 48 um grid), optionally cropped to
    the grid rows/cols that have a vertex inside box = (lo_xyz, hi_xyz). NaN where invalid.
    Only a band of grid rows is downloaded. The band is chosen from the 9.362 um copy of the same mesh: its rows with
    a vertex within 500 um of the box, scaled to the 2.399 grid, plus 16 rows on each side (the two copies agree to
    about 320 um, about 7 rows). The crop is then the same as on the whole tif; an assert fails if the band was too
    narrow for it."""
    _, H, _ = _tif_layout(w, "x")
    if box is None:
        r0, r1 = 0, H - 1
    else:
        lo, hi = np.asarray(box[0]), np.asarray(box[1])
        g9, pad = mesh24(w), 500 / UM24
        with np.errstate(invalid="ignore"):
            near = np.all((g9 >= lo - pad) & (g9 <= hi + pad), -1)
        if not near.any():
            return np.zeros((0, 0, 3))
        rows, s = np.where(near.any(1))[0], H / g9.shape[0]
        r0, r1 = max(int(np.floor(rows[0] * s)) - 16, 0), min(int(np.ceil(rows[-1] * s)) + 16, H - 1)
    g = [_tif_rows(w, c, r0, r1) for c in "xyz"]
    v = (g[0] > 0) & (g[1] > 0) & (g[2] > 0) & np.isfinite(g[0])
    p = np.stack(g, -1); p[~v] = np.nan
    if box is not None:
        with np.errstate(invalid="ignore"):
            m = np.all((p >= lo) & (p <= hi), -1)
        if not m.any():
            return p[:0, :0]
        rs, cs = np.where(m.any(1))[0] + r0, np.where(m.any(0))[0]
        a, b = max(rs[0] - 2, 0), min(rs[-1] + 3, H)
        assert a >= r0 and b <= r1 + 1, f"{w}: row band {r0}..{r1} too narrow for crop {a}..{b - 1}"
        p = p[a - r0:b - r0, max(cs[0] - 2, 0):cs[-1] + 3]
    return p

def ray_hits(grid, o, d):
    """All crossings of the ray o + t d (t>0) with the triangulated grid. Returns sorted list of (t, row, col)."""
    a, b, c, e = grid[:-1, :-1], grid[1:, :-1], grid[:-1, 1:], grid[1:, 1:]
    rr, cc = np.meshgrid(np.arange(grid.shape[0] - 1), np.arange(grid.shape[1] - 1), indexing="ij")
    hits = []
    for (p0, p1, p2) in ((a, b, c), (e, c, b)):
        ok = np.isfinite(p0[..., 0]) & np.isfinite(p1[..., 0]) & np.isfinite(p2[..., 0])
        P0, P1, P2 = p0[ok], p1[ok], p2[ok]
        R, C = rr[ok], cc[ok]
        e1, e2 = P1 - P0, P2 - P0
        h = np.cross(d, e2)
        det = (e1 * h).sum(1)
        good = np.abs(det) > 1e-9
        f = np.where(good, 1.0 / np.where(good, det, 1), 0)
        s = o - P0
        uu = f * (s * h).sum(1)
        q = np.cross(s, e1)
        vv = f * (q * d).sum(1)
        t = f * (q * e2).sum(1)
        m = good & (uu >= 0) & (vv >= 0) & (uu + vv <= 1) & (t > 0)
        for ti, r, cl in zip(t[m], R[m], C[m]):
            hits.append((float(ti), int(r), int(cl)))
    hits.sort()
    # merge duplicates on shared edges
    out = []
    for h in hits:
        if not out or abs(h[0] - out[-1][0]) > 0.5:
            out.append(h)
    return out

NBYTES = [0]
def fetch(key, rng=None):
    """GET one level-0 chunk (or a byte range of it). Cached in DLP. Returns bytes or None on 404."""
    fn = os.path.join(DLP, key.replace("/", "_") + (f"_r{rng[0]}-{rng[1]}" if rng else ""))
    if os.path.exists(fn):
        return open(fn, "rb").read()
    req = urllib.request.Request(BASE + key)
    if rng:
        req.add_header("Range", f"bytes={rng[0]}-{rng[1]}")
    try:
        data = urllib.request.urlopen(req, timeout=120).read()
    except urllib.error.HTTPError as ex:
        if ex.code == 404:
            open(fn + ".404", "wb").close(); return None
        raise
    open(fn, "wb").write(data)
    NBYTES[0] += len(data)
    return data

def dl_bytes():
    return sum(os.path.getsize(os.path.join(DLP, f)) for f in os.listdir(DLP))
