"""Size estimate for finer CT (levels 0, 1, 2) over the cal2 seed region. HTTP HEAD only: no chunk body is read.

Region: trace columns 11 to 17, rows 0 to 2 of the real-CT trace (a bilinear patch of the trace grid, level-0
coordinates), and at every patch point its own local w035 outward normal (cal2_look.load_meshes, as in the ribbons)
from -130 to +60 um. Sampled every 2 um in all three directions.
Margin per point at level L: trilinear reads voxel floor(v) and floor(v) + 1; add 0.5 voxel for a possible
centre/corner offset between levels (not checked in CAL2-LOOK.md) and half the 2 um sampling step.
Chunk key (z, y, x) = floor(level-0 xyz / 2**L / 128), the convention of cal2_look.CT and fetch_ct.py.
Chunk shape, dtype and compressor come from the cached downloads/meta/zarr/<L>.zarray; HEAD Content-Length checks
them. Each level-L key is also mapped to its level-3 parent: if every parent is one of cal2's 30 chunks, the finer
chunks inherit cal2's checks (no control vertices, at least 30 deg from the gap envelope). fetch_ct.outside_gap is
run on every key too, since fetch_ct.py would refuse a key that fails it.
Nothing is appended to logs/fetch_log.jsonl (HEAD reads no body bytes); requests are listed in the output instead.
Usage (WSL venv): python -I cal2_look_finer_estimate.py      Writes out/cal2_look/finer_ct_estimate.json
"""
import json
import math
import os
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cal2_look as L  # noqa: E402
import numpy as np  # noqa: E402
from fetch_ct import DST, VOL, outside_gap  # noqa: E402
from netlog import CAP_BYTES, total_fetched  # noqa: E402

COLS = (11, 17)
ROWS = (0, 2)
S_LO, S_HI = -130.0, 60.0
STEP_UM = 2.0
LEVELS = (0, 1, 2, 3)
SEP_UM = (27.0, 43.0)        # w034 to w035 at the seed, columns 11 to 17 (CAL2-LOOK.md)


def patch(P):
    """Bilinear patch of the trace grid over ROWS x COLS, sampled at most STEP_UM apart."""
    g = P[ROWS[0]:ROWS[1] + 1, COLS[0]:COLS[1] + 1]
    dc = np.linalg.norm(np.diff(g, axis=1), axis=-1).max() * L.UM
    dr = np.linalg.norm(np.diff(g, axis=0), axis=-1).max() * L.UM
    nc = int(math.ceil((COLS[1] - COLS[0]) * dc / STEP_UM)) + 1
    nr = int(math.ceil((ROWS[1] - ROWS[0]) * dr / STEP_UM)) + 1
    cf = np.linspace(0, COLS[1] - COLS[0], nc)
    rf = np.linspace(0, ROWS[1] - ROWS[0], nr)
    c0 = np.minimum(np.floor(cf).astype(int), g.shape[1] - 2)
    r0 = np.minimum(np.floor(rf).astype(int), g.shape[0] - 2)
    fc, fr = (cf - c0)[None, :, None], (rf - r0)[:, None, None]
    a = g[r0][:, c0]
    b = g[r0][:, c0 + 1]
    c = g[r0 + 1][:, c0]
    d = g[r0 + 1][:, c0 + 1]
    base = (1 - fr) * ((1 - fc) * a + fc * b) + fr * ((1 - fc) * c + fc * d)
    return base.reshape(-1, 3), {"n_rows": nr, "n_cols": nc, "max_col_gap_um": round(float(dc), 1),
                                 "max_row_gap_um": round(float(dr), 1)}


def keys_for(pts, lv):
    v = pts / (2 ** lv)                                   # level-L voxel index, xyz
    pad = 0.5 + 0.5 * STEP_UM / (L.UM * 2 ** lv)
    lo = np.floor((v - pad) / L.CH).astype(np.int64)
    hi = np.floor((v + 1 + pad) / L.CH).astype(np.int64)
    out = set()
    for ix in (lo[:, 0], hi[:, 0]):
        for iy in (lo[:, 1], hi[:, 1]):
            for iz in (lo[:, 2], hi[:, 2]):
                u = np.unique(np.stack([iz, iy, ix], 1), axis=0)
                out.update(map(tuple, u.tolist()))
    vmin, vmax = (v - pad).min(0), (v + 1 + pad).max(0)
    return out, vmin, vmax


def head(url):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "pherc0139-smoke/0.1"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.status, int(r.headers.get("Content-Length", -1))
        except urllib.error.HTTPError as exc:
            return exc.code, None
        except Exception as exc:  # transient network errors
            if attempt == 3:
                return f"error: {type(exc).__name__}", None
            time.sleep(2 + 3 * attempt)


def main():
    p, keys3 = L.load_plan()
    keys3 = set(keys3)
    meshes, normal_at = L.load_meshes(p)
    P, kept, _ = L.load_trace("real")
    assert kept[ROWS[0]:ROWS[1] + 1, COLS[0]:COLS[1] + 1].all()
    base, pinfo = patch(P)
    n = normal_at(base)
    s = np.arange(S_LO, S_HI + 1e-9, STEP_UM)
    umb11 = json.load(open(os.path.join(L.ROOT, "out", "locate.json")))["umb11_xyz"]
    acc = {lv: [set(), np.full(3, np.inf), np.full(3, -np.inf)] for lv in LEVELS}
    for i in range(0, len(base), 2000):
        pts = (base[i:i + 2000, None, :] + (s / L.UM)[None, :, None] * n[i:i + 2000, None, :]).reshape(-1, 3)
        for lv in LEVELS:
            k, vmin, vmax = keys_for(pts, lv)
            acc[lv][0] |= k
            acc[lv][1] = np.minimum(acc[lv][1], vmin)
            acc[lv][2] = np.maximum(acc[lv][2], vmax)
    already = total_fetched()
    res = {"region": {"place": L.NAME, "trace_columns": list(COLS), "trace_rows": list(ROWS),
                      "normal_s_um": [S_LO, S_HI], "sample_step_um": STEP_UM, "patch": pinfo,
                      "n_patch_points": int(len(base)),
                      "margin": "per point: voxel floor and floor + 1 (trilinear), plus 0.5 voxel (possible "
                                "centre/corner offset between levels) and half the sampling step",
                      "key_convention": "(z, y, x) = floor(level-0 xyz / 2**L / 128), as cal2_look.CT and "
                                        "fetch_ct.py"},
           "source": {"volume": VOL, "zarray": "cached downloads/meta/zarr/<L>.zarray (not re-fetched)",
                      "method": "HTTP HEAD per chunk (Content-Length); no chunk body read; nothing added to "
                                "logs/fetch_log.jsonl"},
           "ledger": {"cap_bytes": CAP_BYTES, "fetched_bytes": already, "left_bytes": CAP_BYTES - already},
           "levels": {}}
    n_head = 0
    for lv in LEVELS:
        za = json.load(open(os.path.join(L.ROOT, "downloads", "meta", "zarr", f"{lv}.zarray")))
        chunk_bytes = int(np.prod(za["chunks"])) * np.dtype(za["dtype"]).itemsize
        keys = sorted(acc[lv][0])
        nmax = [int(math.ceil(sh / ch)) for sh, ch in zip(za["shape"], za["chunks"])]
        assert all(0 <= k[a] < nmax[a] for k in keys for a in range(3)), "key outside the array"
        parents = {(z >> (3 - lv), y >> (3 - lv), x >> (3 - lv)) for z, y, x in keys}
        vox_um = L.UM * 2 ** lv
        ent = {"voxel_um": round(vox_um, 3), "chunk_shape": za["chunks"], "dtype": za["dtype"],
               "compressor": za["compressor"], "bytes_per_chunk_from_zarray": chunk_bytes,
               "n_chunks": len(keys),
               "region_bbox_voxels_xyz": [np.floor(acc[lv][1]).astype(int).tolist(),
                                          np.ceil(acc[lv][2]).astype(int).tolist()],
               "level3_parents_all_in_cal2_box": parents <= keys3,
               "fetch_ct_outside_gap_all": all(outside_gap(lv, k, umb11) for k in keys),
               "sheet_separation_voxels": [round(x / vox_um, 1) for x in SEP_UM]}
        if lv == 3:
            ent["note"] = "reference only: cal2's level-3 chunks are already on disk"
            ent["already_on_disk"] = sum(os.path.exists(os.path.join(DST, "3", *map(str, k))) for k in keys)
            res["levels"][str(lv)] = ent
            continue
        urls = [f"{VOL}/{lv}/{k[0]}/{k[1]}/{k[2]}" for k in keys]
        with ThreadPoolExecutor(6) as ex:
            hs = list(ex.map(head, urls))
        n_head += len(urls)
        present = [(k, ln) for k, (st, ln) in zip(keys, hs) if st == 200]
        missing = [[list(k), st] for k, (st, _) in zip(keys, hs) if st in (403, 404)]
        errors = [[list(k), st] for k, (st, _) in zip(keys, hs) if st not in (200, 403, 404)]
        on_disk = [k for k, _ in present
                   if os.path.exists(os.path.join(DST, str(lv), *map(str, k)))
                   and os.path.getsize(os.path.join(DST, str(lv), *map(str, k))) == chunk_bytes]
        total = sum(ln for _, ln in present)
        ent.update({"present_200": len(present), "missing_403_404": missing, "errors": errors,
                    "content_length_values": sorted({ln for _, ln in present}),
                    "total_bytes": total,
                    "already_on_disk": len(on_disk),
                    "bytes_to_fetch": total - chunk_bytes * len(on_disk),
                    "fits_in_ledger_left": total - chunk_bytes * len(on_disk) <= CAP_BYTES - already,
                    "keys_zyx": [list(k) for k in keys]})
        res["levels"][str(lv)] = ent
        print(f"level {lv}: voxel {vox_um:.3f} um | {len(keys)} chunks | 200: {len(present)} | missing "
              f"{len(missing)} | errors {len(errors)} | total {total:,} B | on disk {len(on_disk)} | to fetch "
              f"{ent['bytes_to_fetch']:,} B | sep {ent['sheet_separation_voxels']} vox | parents in cal2 "
              f"{ent['level3_parents_all_in_cal2_box']} | outside_gap {ent['fetch_ct_outside_gap_all']}")
    l3 = res["levels"]["3"]
    print(f"level 3 (reference): {l3['n_chunks']} chunks, all in cal2 box: {l3['level3_parents_all_in_cal2_box']}, "
          f"on disk {l3['already_on_disk']}")
    res["source"]["head_requests"] = n_head
    print(f"ledger left {CAP_BYTES - already:,} B; HEAD requests {n_head}")
    json.dump(res, open(os.path.join(L.OUT, "finer_ct_estimate.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
