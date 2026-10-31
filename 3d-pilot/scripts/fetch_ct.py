"""Fetch whole CT chunks of the 1.129 um zarr at one level into downloads/ct/vol.zarr/<level>/<z>/<y>/<x>.

Usage: python -I fetch_ct.py <plan_json>   (plan from plan_ct.py: {"level": L, "chunks": [[z, y, x], ...]})
Chunks that do not exist on S3 (404, masked-out space) are listed, not created.
Refuses any chunk whose region touches the gap sector (327 to 85 deg about the umbilicus) at this z.
"""
import json
import math
import os
import shutil
import sys
import urllib.error
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import fetch, ROOT  # noqa: E402

VOL = ("https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc0139/volumes/"
       "20260413113053-1.129um-0.2m-59keV-masked.zarr")
DST = os.path.join(ROOT, "downloads", "ct", "vol.zarr")


def outside_gap(level, key, umb11):
    """True if every corner of the chunk is at least 30 deg outside the 327..85 deg gap envelope."""
    s = 128 * (2 ** level)
    z, y, x = key
    for dy in (0, 1):
        for dx in (0, 1):
            a = math.degrees(math.atan2((y + dy) * s - umb11[1], (x + dx) * s - umb11[0])) % 360
            # angles here are in the 1.129 frame; the frames differ by a near-identity rotation (< 0.1 deg)
            if a >= 327 - 30 or a <= 85 + 30:
                return False
    return True


def main():
    plan = json.load(open(sys.argv[1]))
    L = int(plan["level"])
    loc = json.load(open(os.path.join(ROOT, "out", "locate.json")))
    umb11 = loc["umb11_xyz"]
    os.makedirs(os.path.join(DST, str(L)), exist_ok=True)
    shutil.copy(os.path.join(ROOT, "downloads", "meta", "zarr", f"{L}.zarray"), os.path.join(DST, str(L), ".zarray"))
    todo = []
    for key in plan["chunks"]:
        key = tuple(int(v) for v in key)
        if not outside_gap(L, key, umb11):
            raise SystemExit(f"chunk {key} at level {L} is within 30 deg of the gap sector; refusing")
        fn = os.path.join(DST, str(L), *map(str, key))
        if os.path.exists(fn) and os.path.getsize(fn) == 2097152:
            continue
        todo.append((key, fn))

    def one(item):
        key, fn = item
        try:
            data = fetch(f"{VOL}/{L}/{key[0]}/{key[1]}/{key[2]}", expect_len=2097152, tag=f"ct-L{L}")
        except urllib.error.HTTPError as exc:
            if exc.code in (403, 404):
                return key, "absent"
            raise
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        with open(fn, "wb") as f:
            f.write(data)
        return key, "ok"

    with ThreadPoolExecutor(6) as ex:
        res = list(ex.map(one, todo))
    absent = [k for k, s in res if s == "absent"]
    print(f"level {L}: {len(plan['chunks'])} chunks planned, fetched {len(res) - len(absent)}, absent {len(absent)}, "
          f"already present {len(plan['chunks']) - len(todo)}")
    if absent:
        print("absent:", absent)


if __name__ == "__main__":
    main()
