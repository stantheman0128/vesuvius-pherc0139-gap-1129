"""Level-0 size estimate for the cal2 seed, made the same way as cal2_look_finer_estimate.py (its patch, keys_for and
head functions are imported, not copied), for level 0 only. HTTP HEAD only: no chunk body is read and nothing is
added to logs/fetch_log.jsonl.

Budget (LEVEL0-ADDENDUM.md): 160,000,000 bytes for this step; the ledger may reach 1,633,529,433 bytes.
Usage (WSL venv): python -I cal2_look_level0_estimate.py      Writes out/cal2_look/level0_ct_estimate.json
"""
import json
import math
import os
import sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cal2_look as L  # noqa: E402
import cal2_look_finer_estimate as E  # noqa: E402
import numpy as np  # noqa: E402
from fetch_ct import VOL, outside_gap  # noqa: E402
from netlog import LOG  # noqa: E402

LV = 0
STEP_BUDGET = 160_000_000
LEDGER_LIMIT = 1_633_529_433
NEW_DST = os.path.join(L.ROOT, "downloads", "ct_level0")


def ledger_total():
    """Fresh read of logs/fetch_log.jsonl (not netlog.total_fetched, see FINER-LOOK.md on the stale read)."""
    n = lines = 0
    with open(LOG, "rt") as f:
        for line in f:
            line = line.strip()
            if line:
                n += int(json.loads(line)["bytes"])
                lines += 1
    return n, lines


def main():
    p, keys3 = L.load_plan()
    keys3 = set(keys3)
    meshes, normal_at = L.load_meshes(p)
    P, kept, _ = L.load_trace("real")
    assert kept[E.ROWS[0]:E.ROWS[1] + 1, E.COLS[0]:E.COLS[1] + 1].all()
    base, pinfo = E.patch(P)
    n = normal_at(base)
    s = np.arange(E.S_LO, E.S_HI + 1e-9, E.STEP_UM)
    umb11 = json.load(open(os.path.join(L.ROOT, "out", "locate.json")))["umb11_xyz"]
    keys, vmin, vmax = set(), np.full(3, np.inf), np.full(3, -np.inf)
    for i in range(0, len(base), 2000):
        pts = (base[i:i + 2000, None, :] + (s / L.UM)[None, :, None] * n[i:i + 2000, None, :]).reshape(-1, 3)
        k, a, b = E.keys_for(pts, LV)
        keys |= k
        vmin, vmax = np.minimum(vmin, a), np.maximum(vmax, b)
    keys = sorted(keys)
    za = json.load(open(os.path.join(L.ROOT, "downloads", "meta", "zarr", f"{LV}.zarray")))
    chunk_bytes = int(np.prod(za["chunks"])) * np.dtype(za["dtype"]).itemsize
    nmax = [int(math.ceil(sh / ch)) for sh, ch in zip(za["shape"], za["chunks"])]
    assert all(0 <= k[a] < nmax[a] for k in keys for a in range(3)), "key outside the array"
    parents = {(z >> 3, y >> 3, x >> 3) for z, y, x in keys}
    urls = [f"{VOL}/{LV}/{k[0]}/{k[1]}/{k[2]}" for k in keys]
    with ThreadPoolExecutor(6) as ex:
        hs = list(ex.map(E.head, urls))
    present = [(k, ln) for k, (st, ln) in zip(keys, hs) if st == 200]
    missing = [[list(k), st] for k, (st, _) in zip(keys, hs) if st in (403, 404)]
    errors = [[list(k), st] for k, (st, _) in zip(keys, hs) if st not in (200, 403, 404)]
    total = sum(ln for _, ln in present)
    ledger, lines = ledger_total()
    old = json.load(open(os.path.join(L.OUT, "finer_ct_estimate.json")))["levels"]["0"]["keys_zyx"]
    res = {"what": "level-0 HEAD estimate for LEVEL0-ADDENDUM.md, same method as cal2_look_finer_estimate.py",
           "region": {"place": L.NAME, "trace_columns": list(E.COLS), "trace_rows": list(E.ROWS),
                      "normal_s_um": [E.S_LO, E.S_HI], "sample_step_um": E.STEP_UM, "patch": pinfo,
                      "n_patch_points": int(len(base))},
           "source": {"volume": VOL, "method": "HTTP HEAD per chunk (Content-Length); no chunk body read",
                      "head_requests": len(urls)},
           "ledger": {"total_bytes_now": ledger, "lines": lines, "limit_bytes": LEDGER_LIMIT,
                      "step_budget_bytes": STEP_BUDGET},
           "level": {"level": LV, "voxel_um": L.UM, "chunk_shape": za["chunks"], "dtype": za["dtype"],
                     "compressor": za["compressor"], "bytes_per_chunk_from_zarray": chunk_bytes,
                     "n_chunks": len(keys),
                     "region_bbox_voxels_xyz": [np.floor(vmin).astype(int).tolist(),
                                                np.ceil(vmax).astype(int).tolist()],
                     "level3_parents_all_in_cal2_box": parents <= keys3,
                     "fetch_ct_outside_gap_all": all(outside_gap(LV, k, umb11) for k in keys),
                     "present_200": len(present), "missing_403_404": missing, "errors": errors,
                     "content_length_values": sorted({ln for _, ln in present}),
                     "total_bytes": total,
                     "same_keys_as_finer_ct_estimate": sorted(map(tuple, old)) == keys,
                     "keys_zyx": [list(k) for k in keys]},
           "fetch_dir": os.path.relpath(NEW_DST, L.ROOT),
           "fetch_dir_exists_now": os.path.exists(NEW_DST)}
    fits = (not errors) and total <= STEP_BUDGET and ledger + total <= LEDGER_LIMIT
    res["fits_budget"] = fits
    res["ledger_after_if_fetched"] = ledger + total
    json.dump(res, open(os.path.join(L.OUT, "level0_ct_estimate.json"), "w"), indent=1)
    print(f"level 0: {len(keys)} chunks | 200: {len(present)} | missing {len(missing)} | errors {len(errors)} | "
          f"total {total:,} B | CL {res['level']['content_length_values']} | parents in cal2 "
          f"{res['level']['level3_parents_all_in_cal2_box']} | outside_gap {res['level']['fetch_ct_outside_gap_all']}"
          f" | same keys as finer estimate {res['level']['same_keys_as_finer_ct_estimate']}")
    print(f"ledger now {ledger:,} B ({lines} lines); after {ledger + total:,} B; limit {LEDGER_LIMIT:,}; "
          f"step budget {STEP_BUDGET:,}; fits {fits}")


if __name__ == "__main__":
    main()
