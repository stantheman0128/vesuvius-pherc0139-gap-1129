"""Fetch the 14 finer-CT chunks of finer_ct_estimate.json (levels 1 and 2 only) through the existing fetch path.

Same request, tag and file layout as fetch_ct.py (netlog.fetch -> logs/fetch_log.jsonl; downloads/ct/vol.zarr/<L>/z/y/x),
but sequential, and the ledger is checked before every chunk: stop if the total would pass netlog.CAP_BYTES.
fetch_ct.outside_gap is checked for every key, and each key's level-3 parent must be one of cal2's chunks.
Usage (WSL venv): python -I cal2_look_finer_fetch.py
"""
import json
import os
import shutil
import sys
import urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cal2_look as L  # noqa: E402
from fetch_ct import DST, VOL, outside_gap  # noqa: E402
from netlog import CAP_BYTES, ROOT, fetch, total_fetched  # noqa: E402

CHUNK = 2097152
LEVELS = (1, 2)


def main():
    est = json.load(open(os.path.join(L.OUT, "finer_ct_estimate.json")))
    _, keys3 = L.load_plan()
    keys3 = set(keys3)
    umb11 = json.load(open(os.path.join(ROOT, "out", "locate.json")))["umb11_xyz"]
    start = total_fetched()
    print(f"ledger before: {start:,} B")
    got = 0
    for lv in LEVELS:
        keys = [tuple(k) for k in est["levels"][str(lv)]["keys_zyx"]]
        os.makedirs(os.path.join(DST, str(lv)), exist_ok=True)
        za = os.path.join(DST, str(lv), ".zarray")
        if not os.path.exists(za):
            shutil.copy(os.path.join(ROOT, "downloads", "meta", "zarr", f"{lv}.zarray"), za)
        for key in keys:
            assert (key[0] >> (3 - lv), key[1] >> (3 - lv), key[2] >> (3 - lv)) in keys3, key
            if not outside_gap(lv, key, umb11):
                raise SystemExit(f"chunk {key} at level {lv} is within 30 deg of the gap sector; refusing")
            fn = os.path.join(DST, str(lv), *map(str, key))
            if os.path.exists(fn) and os.path.getsize(fn) == CHUNK:
                print(f"L{lv} {key} already on disk")
                continue
            now = total_fetched()
            if now + CHUNK > CAP_BYTES:
                raise SystemExit(f"stop: ledger {now:,} + {CHUNK:,} would pass the cap {CAP_BYTES:,}")
            try:
                data = fetch(f"{VOL}/{lv}/{key[0]}/{key[1]}/{key[2]}", expect_len=CHUNK, tag=f"ct-L{lv}")
            except urllib.error.HTTPError as exc:
                print(f"L{lv} {key} HTTP {exc.code}: absent")
                continue
            os.makedirs(os.path.dirname(fn), exist_ok=True)
            with open(fn, "wb") as f:
                f.write(data)
            got += len(data)
            print(f"L{lv} {key} ok; ledger {total_fetched():,} B")
    end = total_fetched()
    print(f"fetched this run: {got:,} B; ledger after: {end:,} B (was {start:,}); left {CAP_BYTES - end:,} B")


if __name__ == "__main__":
    main()
