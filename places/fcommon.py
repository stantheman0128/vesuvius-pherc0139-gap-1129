# Shared helpers for the #2167 five-place run on the 1.129 um scan (generalises the E2 pilot to any place).
# Reuses the pilot's frames, mesh readers and ray-triangle code (pilot/pcommon.py); chunk downloads go to
# $PHERC0139_DATA/chunks-places. Run with python -I.
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pilot"))
from pcommon import np, umb, to11, to24, mesh24, mesh24_fine, ray_hits, UM11, UM24, BASE
import urllib.request
from concurrent.futures import ThreadPoolExecutor

from repo_paths import DATA, OUT_PLACES as OUT5
DL5 = os.path.join(DATA, "chunks-places")
os.makedirs(DL5, exist_ok=True)
CAP = 600_000_000
CS = 128
MESHES = ("w033", "w034", "w035", "w036", "w037")

# name: (z in 2.399 frame, angle deg or None, xy of the place or None)
PLACES = {
    "E1": (14856.0, None, (16642.0, 14927.0)),
    "E5": (13741.0, None, (17427.0, 13862.0)),
    "E6": (14109.0, None, (17029.0, 14405.0)),
    "X434": (12837.0, None, (17115.0, 14541.0)),    # w034 row 434 col 248, note example "in the gap"
    "X542": (18432.0, None, (18213.0, 11594.0)),    # w036 row 542 col 1016, note example "in the gap"
    "G2": (19500.0, 15.0, None),
}

def place_frame(name):
    z0, a, xy = PLACES[name]
    xu, yu = umb(z0)
    if a is None:
        a = float(np.degrees(np.arctan2(xy[1] - yu, xy[0] - xu)) % 360)
    d = np.array([np.cos(np.radians(a)), np.sin(np.radians(a)), 0.0])
    t = np.array([-d[1], d[0], 0.0])
    return z0, a, np.array([xu, yu, z0]), d, t

def p24(o, d, t, r_um, tang_um, dz_um):
    return (o + np.multiply.outer(r_um / UM24, d) + np.multiply.outer(tang_um / UM24, t)
            + np.multiply.outer(dz_um / UM24, [0, 0, 1.0]))

def dl_bytes():
    return sum(os.path.getsize(os.path.join(DL5, f)) for f in os.listdir(DL5))

def needed(pts11):
    """{(zc,yc,xc): [p0,p1]} z-plane range per chunk needed for trilinear sampling at the 1.129 points (N,3 xyz)."""
    p = pts11.reshape(-1, 3)
    lo = np.floor(p).astype(np.int64); hi = lo + 1
    keys = []
    for zz in (lo[:, 2], hi[:, 2]):
        for yy in (lo[:, 1], hi[:, 1]):
            for xx in (lo[:, 0], hi[:, 0]):
                keys.append(np.unique(((zz // CS * 1000 + yy // CS) * 1000 + xx // CS) * CS + zz % CS))
    u = np.unique(np.concatenate(keys))
    ck, pl = u // CS, u % CS
    uc, inv = np.unique(ck, return_inverse=True)
    mn = np.full(len(uc), 999); mx = np.full(len(uc), -1)
    np.minimum.at(mn, inv.ravel(), pl); np.maximum.at(mx, inv.ravel(), pl)
    return {(int(c // 1000000), int(c // 1000 % 1000), int(c % 1000)): [int(a), int(b)] for c, a, b in zip(uc, mn, mx)}

def merge_needs(*ns):
    out = {}
    for n in ns:
        for k, (a, b) in n.items():
            out[k] = [min(out[k][0], a), max(out[k][1], b)] if k in out else [a, b]
    return out

def plan_bytes(need):
    return sum((b - a + 1) * CS * CS for a, b in need.values())

def _fn(k, a, b):
    return os.path.join(DL5, f"{k[0]}_{k[1]}_{k[2]}_p{a}-{b}")

def have_planes(k, a, b):
    """cached file covering planes a..b of chunk k, or None"""
    pre = f"{k[0]}_{k[1]}_{k[2]}_p"
    for f in os.listdir(DL5):
        if f.startswith(pre) and not f.endswith(".404"):
            pa, pb = map(int, f[len(pre):].split("-"))
            if pa <= a and pb >= b:
                return os.path.join(DL5, f), pa
    return None

NB = [0]
def fetch_one(k, a, b):
    h = have_planes(k, a, b)
    if h:
        return "cached"
    req = urllib.request.Request(BASE + f"{k[0]}/{k[1]}/{k[2]}")
    if not (a == 0 and b == CS - 1):
        req.add_header("Range", f"bytes={a*CS*CS}-{(b+1)*CS*CS-1}")
    try:
        data = urllib.request.urlopen(req, timeout=120).read()
    except urllib.error.HTTPError as ex:
        if ex.code == 404:
            open(_fn(k, a, b) + ".404", "wb").close(); return "404"
        raise
    if len(data) != (b - a + 1) * CS * CS:
        return f"badsize {len(data)}"
    open(_fn(k, a, b), "wb").write(data)
    NB[0] += len(data)
    return "ok"

def fetch_need(need, label=""):
    todo = {k: v for k, v in need.items() if not have_planes(k, *v)}
    est = plan_bytes(todo)
    cur = dl_bytes()
    print(f"[{label}] chunks {len(need)}, to fetch {len(todo)} = {est/1e6:.1f} MB; dir now {cur/1e6:.1f} MB; cap {CAP/1e6:.0f} MB", flush=True)
    if cur + est > CAP:
        raise SystemExit(f"CAP: would reach {(cur+est)/1e6:.1f} MB")
    with ThreadPoolExecutor(8) as ex:
        res = list(ex.map(lambda kv: fetch_one(kv[0], *kv[1]), todo.items()))
    bad = [(k, r) for k, r in zip(todo, res) if r not in ("ok", "cached")]
    print(f"[{label}] fetched {sum(r == 'ok' for r in res)}; problems {bad[:10]}; dir now {dl_bytes()/1e6:.1f} MB", flush=True)
    return bad

def assemble(need):
    """dense uint8 subvolume over the bbox of the needed chunks/planes; returns vol (z,y,x), origin xyz"""
    ks = list(need)
    z0 = min(k[0] * CS + need[k][0] for k in ks); z1 = max(k[0] * CS + need[k][1] for k in ks)
    y0 = min(k[1] for k in ks) * CS; y1 = (max(k[1] for k in ks) + 1) * CS
    x0 = min(k[2] for k in ks) * CS; x1 = (max(k[2] for k in ks) + 1) * CS
    vol = np.zeros((z1 - z0 + 1, y1 - y0, x1 - x0), np.uint8)
    for k, (a, b) in need.items():
        h = have_planes(k, a, b)
        if h is None:
            continue
        fn, pa = h
        blk = np.frombuffer(open(fn, "rb").read(), np.uint8).reshape(-1, CS, CS)
        blk = blk[a - pa:b - pa + 1]
        zs = k[0] * CS + a - z0
        vol[zs:zs + blk.shape[0], k[1] * CS - y0:(k[1] + 1) * CS - y0, k[2] * CS - x0:(k[2] + 1) * CS - x0] = blk
    return vol, np.array([x0, y0, z0], float)

def sample(vol, org, pts11):
    from scipy.ndimage import map_coordinates
    q = pts11.reshape(-1, 3) - org
    v = map_coordinates(vol, [q[:, 2], q[:, 1], q[:, 0]], order=1, mode="constant", cval=np.nan, output=np.float32)
    return v.reshape(pts11.shape[:-1])
