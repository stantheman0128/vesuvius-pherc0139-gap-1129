"""CT for the four calibration places picked by pick_places.py (out/places/places.json "calibration"), and the
real / +1 / -1 voxel views for the sheet tracer. Never touches a control: every chunk is checked against the
control regions, and every chunk corner must be at least 30 deg outside the 327..85 deg envelope (per-height
umbilicus, 2.399 um frame).

Seed: the place's w035 column line, trimmed to 361..375 level-3 voxels (C240's was 364.9; 7 front points at the
checkpoint's point spacing 60). CT chunks: the union over the 8 tracer reads of the read box the tracer will ask for
(inference.py read_block: front +- view_radius + step + 4 = 60.5 level-3 voxels), with the front predicted from the
w035 mesh advanced 20 level-3 voxels per step along the arc, plus an 8-voxel margin.
Shift: shifted[p] = real[p - t], t = +-1 level-3 voxel (9.03 um) along the zyx axis nearest the place's mean outward
w035 normal. Voxels whose source is outside the downloaded chunks take the median; analysis checks reads avoid them.
Usage: python -I calib_prep.py plan|fetch|views      Writes out/calib/<place>/...
"""
import json
import math
import os
import shutil
import site
import sys
import threading
import urllib.error
from concurrent.futures import ThreadPoolExecutor

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
import tifffile  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import netlog  # noqa: E402
import pick_places as PP  # noqa: E402

ROOT = netlog.ROOT
OUT = os.path.join(ROOT, "out", "calib")
VOL = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc0139/volumes/"
       "20260413113053-1.129um-0.2m-59keV-masked.zarr")
DST = os.path.join(ROOT, "downloads", "ct", "vol.zarr")
CH = 128
READ_PAD = math.sqrt(11.5 ** 2 + 23.5 ** 2 + (15.5 + 10.0) ** 2) + 20.0 + 4.0   # view_radius + step + 4 = 60.5
MARGIN = 8.0
UM = 1.129


def outward_normals(g, umb_xy_of):
    n = np.cross(np.gradient(g, axis=1), np.gradient(g, axis=0))
    n /= np.linalg.norm(n, axis=-1, keepdims=True) + 1e-12
    flip = np.sum(n[..., :2] * (g[..., :2] - umb_xy_of(g)), -1) < 0
    n[flip] *= -1
    return n


def umb11_xy(g):
    """Umbilicus (x, y) in the 1.129 frame at each point's height (2.399-frame control points mapped back)."""
    q = PP.to24(g)
    ux, uy = np.interp(q[..., 2], PP.UZ, PP.UX), np.interp(q[..., 2], PP.UZ, PP.UY)
    A = np.eye(4)
    A[:3] = PP.T
    Ai = np.linalg.inv(A)
    u24 = np.stack([ux, uy, q[..., 2]], -1)
    u11 = u24 @ Ai[:3, :3].T + Ai[:3, 3]
    return u11[..., :2]


def places():
    P = json.load(open(os.path.join(ROOT, "out", "places", "places.json")))
    res = []
    for k, pick in enumerate(P["calibration"]):
        Zh = pick["height"]
        reg = P["heights"][str(Zh)]["calib_candidates"][pick["candidate_index"]]
        name = f"cal{k + 1}_z{Zh}_a{reg['centre_grid_deg']}"
        res.append((name, int(Zh), reg, pick))
    return P, res


def plan():
    P, pl = places()
    allctl = []
    for Zh, h in P["heights"].items():
        d = np.load(os.path.join(ROOT, "downloads", "mesh_rows", f"w035_z{Zh}.npz"))
        for r in h["controls"].values():
            if r:
                allctl.append(d["g"][r["band_rows"][0]:r["band_rows"][1], r["cols"][0]:r["cols"][1]].reshape(-1, 3))
    allctl = np.concatenate(allctl).astype(np.float64) / 8.0
    plans = []
    for name, Zh, reg, pick in pl:
        d = np.load(os.path.join(ROOT, "downloads", "mesh_rows", f"w035_z{Zh}.npz"))
        G, rows = d["g"].astype(np.float64), d["rows"]
        cs, sgn = reg["seed_col"], (1 if reg["direction"] == "+" else -1)
        r0, r1 = reg["band_rows"]
        while True:
            L = np.linalg.norm(np.diff(G[r0:r1, cs] / 8.0, axis=0), axis=1).sum()
            if L <= 375.0:
                break
            if (r1 - r0) % 2:
                r1 -= 1
            else:
                r0 += 1
        assert L >= 361.0, (name, L)
        cols_per_step = 20.0 * 8 * UM / reg["col_step_um"]
        chunks, boxes = set(), []
        for t in range(8):      # reads happen at the front of steps 1..8, i.e. after t = 0..7 advances
            c = cs + sgn * t * cols_per_step
            c0 = int(math.floor(c))
            f = c - c0
            line = (1 - f) * G[r0:r1, c0] + f * G[r0:r1, min(c0 + 1, G.shape[1] - 1)]
            assert np.all(line > -1), (name, t)
            p = line / 8.0
            lo, hi = p.min(0) - READ_PAD - MARGIN, p.max(0) + READ_PAD + MARGIN
            boxes.append([lo.round(1).tolist(), hi.round(1).tolist()])
            for z in range(int(lo[2] // CH), int(hi[2] // CH) + 1):
                for y in range(int(lo[1] // CH), int(hi[1] // CH) + 1):
                    for x in range(int(lo[0] // CH), int(hi[0] // CH) + 1):
                        chunks.add((z, y, x))
        chunks = sorted(chunks)
        guard = PP.chunk_guard(chunks)
        ck = np.array(chunks)
        ctl_in = 0
        for (z, y, x) in chunks:
            lo = np.array([x, y, z]) * CH
            ctl_in += int(np.all((allctl >= lo) & (allctl < lo + CH), -1).sum())
        # shift axis: nearest the mean outward normal over the region
        g = G[r0:r1, min(cs, reg["end_col"]):max(cs, reg["end_col"]) + 1]
        n = outward_normals(g, umb11_xy)
        nm = n.reshape(-1, 3).mean(0)
        nm /= np.linalg.norm(nm)
        ax_xyz = int(np.argmax(np.abs(nm)))
        pl_ = {"name": name, "height": Zh, "kind": pick["kind"], "p10_nearer_um": pick["p10_nearer_um"],
               "seed_col": cs, "direction": reg["direction"], "seed_band_rows": [int(r0), int(r1)],
               "seed_global_rows": [int(rows[r0]), int(rows[r1 - 1])], "seed_len_l3": round(float(L), 1),
               "line": f"col:{cs}", "range": f"{r0}:{r1}", "angle_range": [reg["angle_min"], reg["angle_max"]],
               "chunks": [list(c) for c in chunks], "n_chunks": len(chunks), "n_chunks_global_box": reg["n_chunks"],
               "guard_min_deg_from_envelope": round(guard, 1), "control_vertices_in_chunks": ctl_in,
               "read_boxes_xyz_l3": boxes, "n_mean_xyz_outward": nm.round(4).tolist(), "axis_xyz": "xyz"[ax_xyz],
               "outward_sign_on_axis": int(np.sign(nm[ax_xyz])),
               "angle_normal_to_axis_deg": round(float(np.degrees(np.arccos(abs(nm[ax_xyz])))), 1)}
        assert guard >= 30 and ctl_in == 0, pl_
        plans.append(pl_)
        print(f"{name}: seed {L:.1f} L3, chunks {len(chunks)} (global box {reg['n_chunks']}), guard {guard:.1f} deg, "
              f"axis {'xyz'[ax_xyz]} ({pl_['angle_normal_to_axis_deg']} deg from normal)")
    total = sum(p["n_chunks"] for p in plans)
    have = sum(1 for p in plans for c in p["chunks"]
               if os.path.exists(os.path.join(DST, "3", *map(str, c))))
    print(f"total chunks {total} ({total * 2097152} bytes), already on disk {have}; ledger {netlog.total_fetched()}, "
          f"room {netlog.CAP_BYTES - netlog.total_fetched()}")
    os.makedirs(OUT, exist_ok=True)
    json.dump(plans, open(os.path.join(OUT, "plans.json"), "w"), indent=1)


def fetch():
    plans = json.load(open(os.path.join(OUT, "plans.json")))
    base = netlog.total_fetched()
    acc, lk = [0], threading.Lock()
    netlog.total_fetched = lambda: base + acc[0]
    os.makedirs(os.path.join(DST, "3"), exist_ok=True)
    if not os.path.exists(os.path.join(DST, "3", ".zarray")):
        shutil.copy(os.path.join(ROOT, "downloads", "meta", "zarr", "3.zarray"), os.path.join(DST, "3", ".zarray"))
    todo = sorted({tuple(c) for p in plans for c in p["chunks"]
                   if not (os.path.exists(os.path.join(DST, "3", *map(str, c))) and
                           os.path.getsize(os.path.join(DST, "3", *map(str, c))) == 2097152)})
    need = len(todo) * 2097152
    if base + need > netlog.CAP_BYTES:
        raise SystemExit(f"cap: {base} + {need} > {netlog.CAP_BYTES}; not fetching")
    print(f"fetching {len(todo)} chunks ({need} bytes); ledger before {base}", flush=True)

    def one(key):
        fn = os.path.join(DST, "3", *map(str, key))
        try:
            data = netlog.fetch(f"{VOL}/3/{key[0]}/{key[1]}/{key[2]}", expect_len=2097152, tag="ct-L3-calib")
        except urllib.error.HTTPError as exc:
            if exc.code in (403, 404):
                return key, "absent"
            raise
        with lk:
            acc[0] += len(data)
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        with open(fn, "wb") as f:
            f.write(data)
        return key, "ok"

    with ThreadPoolExecutor(6) as ex:
        res = list(ex.map(one, todo))
    absent = [k for k, s in res if s == "absent"]
    print(f"fetched {len(res) - len(absent)}, absent {absent}; ledger now {netlog.total_fetched()}")


def views():
    plans = json.load(open(os.path.join(OUT, "plans.json")))
    meta = json.load(open(os.path.join(ROOT, "downloads", "meta", "zarr", "3.zarray")))
    info = []
    for p in plans:
        keys = [tuple(c) for c in p["chunks"]]
        k = np.array(keys)
        lo = k.min(0)
        shape = tuple((k.max(0) - lo + 1) * CH)
        arr = np.zeros(shape, np.uint8)
        have = np.zeros(shape, bool)
        for key in keys:
            b = np.fromfile(os.path.join(DST, "3", *map(str, key)), np.uint8).reshape(CH, CH, CH)
            o = (np.array(key) - lo) * CH
            arr[o[0]:o[0] + CH, o[1]:o[1] + CH, o[2]:o[2] + CH] = b
            have[o[0]:o[0] + CH, o[1]:o[1] + CH, o[2]:o[2] + CH] = True
        med = int(np.median(arr[have]))
        ax_zyx = 2 - "xyz".index(p["axis_xyz"])
        sgn = p["outward_sign_on_axis"]
        pd = os.path.join(OUT, p["name"])
        for tag, d in (("real", 0), ("p", sgn), ("m", -sgn)):
            if d == 0:
                out, valid = arr, have
            else:
                out = np.full_like(arr, med)
                valid = np.zeros_like(have)
                src, dst = [slice(None)] * 3, [slice(None)] * 3
                if d > 0:
                    dst[ax_zyx], src[ax_zyx] = slice(d, None), slice(None, -d)
                else:
                    dst[ax_zyx], src[ax_zyx] = slice(None, d), slice(-d, None)
                out[tuple(dst)] = arr[tuple(src)]
                valid[tuple(dst)] = have[tuple(src)]
                out = np.where(valid & have, out, med).astype(np.uint8)
            vd = os.path.join(pd, f"view_{tag}", "0")
            if os.path.exists(vd):
                shutil.rmtree(vd)
            os.makedirs(vd)
            json.dump(meta, open(os.path.join(vd, ".zarray"), "w"))
            for key in keys:
                o = (np.array(key) - lo) * CH
                dd = os.path.join(vd, str(key[0]), str(key[1]))
                os.makedirs(dd, exist_ok=True)
                np.ascontiguousarray(out[o[0]:o[0] + CH, o[1]:o[1] + CH, o[2]:o[2] + CH]).tofile(
                    os.path.join(dd, str(key[2])))
            if d != 0:
                np.save(os.path.join(pd, f"valid_{tag}.npy"), np.packbits(valid & have))
        t_zyx = [0, 0, 0]
        t_zyx[ax_zyx] = sgn
        t_um_xyz = (np.array(t_zyx[::-1], float) * 8 * UM).tolist()
        # seed tifxyz: the band in level-3 voxels (-1 kept)
        d = np.load(os.path.join(ROOT, "downloads", "mesh_rows", f"w035_z{p['height']}.npz"))
        g = d["g"]
        g3 = np.where(np.all(g != -1, -1, keepdims=True), g / 8.0, -1.0).astype(np.float32)
        sd = os.path.join(pd, "seed_l3")
        os.makedirs(sd, exist_ok=True)
        for kk, c in enumerate("xyz"):
            tifffile.imwrite(os.path.join(sd, f"{c}.tif"), np.ascontiguousarray(g3[..., kk]))
        # read-back check: on-disk shifted chunk at q equals the real CT at q - t (interior of the middle chunk)
        rb = {}
        key = keys[len(keys) // 2]
        o = (np.array(key) - lo) * CH
        for tag, dsign in (("p", sgn), ("m", -sgn)):
            disk = np.fromfile(os.path.join(pd, f"view_{tag}", "0", *map(str, key)), np.uint8).reshape(CH, CH, CH)
            tt = np.zeros(3, int)
            tt[ax_zyx] = dsign
            s0 = o + 1 - tt
            ref = arr[s0[0]:s0[0] + CH - 2, s0[1]:s0[1] + CH - 2, s0[2]:s0[2] + CH - 2]
            rb[tag] = {"equal_to_real_at_q_minus_t": bool(np.array_equal(disk[1:-1, 1:-1, 1:-1], ref)),
                       "differs_from_real": bool(not np.array_equal(disk[1:-1, 1:-1, 1:-1],
                                                                    arr[o[0] + 1:o[0] + CH - 1, o[1] + 1:o[1] + CH - 1,
                                                                        o[2] + 1:o[2] + CH - 1]))}
        info.append({"name": p["name"], "median": med, "box_lo_chunk_zyx": lo.tolist(), "shape_zyx": [int(s) for s in shape],
                     "shift_axis_zyx": ax_zyx, "t_plus_um_xyz": t_um_xyz,
                     "t_plus_dot_n_mean_um": round(float(np.array(t_um_xyz) @ np.array(p["n_mean_xyz_outward"])), 2),
                     "readback": rb})
        print(info[-1], flush=True)
    json.dump(info, open(os.path.join(OUT, "views.json"), "w"), indent=1)


if __name__ == "__main__":
    {"plan": plan, "fetch": fetch, "views": views}[sys.argv[1]]()
