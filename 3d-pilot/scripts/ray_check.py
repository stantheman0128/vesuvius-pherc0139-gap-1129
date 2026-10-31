"""Second check, in the frame of the gap-1129 README: crossings of the C240 ray (umbilicus at z 17,995 of the
2.399 um frame, 240 deg) with w034/w035/w036 and with the copy model's outputs, after mapping the 1.129 um
vertices to the 2.399 um frame with the zarr's transform.json. The README gives w034 8257.8, w035 8382.6,
w036 8550.9 um (from the 2.399 um mesh copies) and CT material runs at dz 0 of 8237-8292, 8319-8360 and
8530-8688 um. Rays at dz = -100..+100 um (2.399 voxels: about +/-42) give a small height profile.
Usage: python -I ray_check.py <run_dir> [<run_dir> ...]
"""
import json
import os
import site
import sys

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
import tifffile  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import ROOT  # noqa: E402

UM24 = 2.399
loc = json.load(open(os.path.join(ROOT, "out", "locate.json")))
T = np.array(json.load(open(os.path.join(ROOT, "downloads", "meta", "zarr", "transform.json")))["transformation_matrix"])
ux, uy = loc["umb24_xy"]
Z0, ANG = 17995.0, 240.0
d = np.array([np.cos(np.radians(ANG)), np.sin(np.radians(ANG)), 0.0])


def to24(g):
    return np.c_[g.reshape(-1, 3), np.ones(g.size // 3)] @ T.T


def hits(g, v, o):
    """Ray o + t d against the grid's triangles (Moller-Trumbore); returns sorted t > 0."""
    p = to24(g).reshape(g.shape)
    q = v[:-1, :-1] & v[1:, :-1] & v[:-1, 1:] & v[1:, 1:]
    out = []
    for a, b, c in ((p[:-1, :-1], p[1:, :-1], p[1:, 1:]), (p[:-1, :-1], p[1:, 1:], p[:-1, 1:])):
        a, b, c = a[q], b[q], c[q]
        e1, e2 = b - a, c - a
        h = np.cross(d, e2)
        det = np.einsum("ij,ij->i", e1, h)
        ok = np.abs(det) > 1e-9
        f = np.where(ok, 1.0 / np.where(ok, det, 1), 0)
        s = o - a
        u = f * np.einsum("ij,ij->i", s, h)
        qq = np.cross(s, e1)
        w = f * (qq @ d)
        t = f * np.einsum("ij,ij->i", e2, qq)
        m = ok & (u >= 0) & (w >= 0) & (u + w <= 1) & (t > 0)
        out += list(t[m])
    return sorted(out)


def load(dname):
    g = np.stack([tifffile.imread(os.path.join(dname, f"{c}.tif")).astype(np.float64) for c in "xyz"], -1)
    return g, np.all(g != -1, -1)


surfs = {}
for w in ("w034", "w035", "w036"):
    b = np.load(os.path.join(ROOT, "downloads", "mesh_bands", f"{w}_band.npz"))["g"].astype(np.float64)
    surfs[w] = (b, np.all(b != -1, -1))
for rd in sys.argv[1:]:
    tag = os.path.basename(rd.rstrip("\\/"))
    for side in ("front", "back"):
        dn = [x for x in os.listdir(rd) if x.endswith("_" + side)][0]
        surfs[f"{tag}:{side}"] = load(os.path.join(rd, dn))

res = {}
dzs = np.arange(-100, 101, 25) / UM24
for name, (g, v) in surfs.items():
    res[name] = {}
    for dz in dzs:
        o = np.array([ux, uy, Z0 + dz])
        ts = hits(g, v, o)
        res[name][f"{dz * UM24:+.0f}"] = [round(t * UM24, 1) for t in ts]
print("crossings r (um) along the C240 ray, by dz (um):")
for name, row in res.items():
    print(f"  {name:22s} " + "  ".join(f"{k}:{','.join(map(str, v)) or '-'}" for k, v in row.items()))
json.dump(res, open(os.path.join(ROOT, "out", "ray_check.json"), "wt"), indent=1)
