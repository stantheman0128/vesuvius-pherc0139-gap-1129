"""Fetch the level-0 chunks of out/cal2_look/level0_ct_estimate.json (LEVEL0-ADDENDUM.md) into a new, empty
directory, downloads/ct_level0/0/z/y/x. Standard library only, for Windows Python (py -I), so the ledger is read
from C:\\ and not through /mnt/c (the stale netlog.total_fetched() read noted in FINER-LOOK.md).

Same request and ledger as fetch_ct.py: netlog.fetch appends one line per chunk to logs/fetch_log.jsonl.
Sequential. Before every chunk the ledger total is recomputed from fetch_log.jsonl, and the fetch stops if it would
pass 1,633,529,433 bytes, or if this run would pass 160,000,000 bytes. After every chunk the file is re-read and the
bytes of this run's own lines are checked against the in-process count; a mismatch (a stale read) stops the run.
fetch_ct.outside_gap is checked for every key, and each key's level-3 parent must be one of cal2's chunks.
Usage (Windows): py -I cal2_look_level0_fetch.py
"""
import ctypes
import json
import os
import sys
import time
import urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fetch_ct import VOL, outside_gap  # noqa: E402
from netlog import LOG, ROOT, fetch  # noqa: E402

LV = 0
CHUNK = 2097152
STEP_BUDGET = 160_000_000
LEDGER_LIMIT = 1_633_529_433
TAG = "ct-L0"
NEW_DST = os.path.join(ROOT, "downloads", "ct_level0")
OUT = os.path.join(ROOT, "out", "cal2_look")


def ledger():
    """(total bytes, bytes of this run's lines, number of lines), from a fresh read of fetch_log.jsonl."""
    tot = own = lines = 0
    with open(LOG, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            e = json.loads(line)
            tot += int(e["bytes"])
            lines += 1
            if e.get("tag") == TAG:
                own += int(e["bytes"])
    return tot, own, lines


def main():
    if os.name == "nt":   # the "nice" of this run: below-normal priority class
        k32 = ctypes.windll.kernel32
        k32.SetPriorityClass(k32.GetCurrentProcess(), 0x4000)
    est = json.load(open(os.path.join(OUT, "level0_ct_estimate.json")))
    assert est["fits_budget"], "estimate is over budget"
    keys = [tuple(k) for k in est["level"]["keys_zyx"]]
    plan = [x for x in json.load(open(os.path.join(ROOT, "out", "calib", "plans.json")))
            if x["name"] == "cal2_z14856_a245"][0]
    keys3 = {tuple(int(v) for v in c) for c in plan["chunks"]}
    umb11 = json.load(open(os.path.join(ROOT, "out", "locate.json")))["umb11_xyz"]
    for key in keys:
        assert (key[0] >> 3, key[1] >> 3, key[2] >> 3) in keys3, key
        if not outside_gap(LV, key, umb11):
            raise SystemExit(f"chunk {key} at level {LV} is within 30 deg of the gap sector; refusing")
    os.makedirs(NEW_DST, exist_ok=False)   # must be new and empty
    start, own0, lines0 = ledger()
    assert own0 == 0, f"ledger already has {own0} bytes tagged {TAG}"
    print(f"{time.strftime('%Y-%m-%d %H:%M:%S %z')} ledger before: {start:,} B, {lines0} lines; "
          f"{len(keys)} chunks planned; dir {NEW_DST}")
    got, first_t, absent = 0, None, []
    for key in keys:
        now, own, _ = ledger()
        if own != got:
            raise SystemExit(f"stop: ledger shows {own:,} B tagged {TAG}, this run counted {got:,} B (stale read?)")
        if now + CHUNK > LEDGER_LIMIT:
            raise SystemExit(f"stop: ledger {now:,} + {CHUNK:,} would pass {LEDGER_LIMIT:,}")
        if got + CHUNK > STEP_BUDGET:
            raise SystemExit(f"stop: this run {got:,} + {CHUNK:,} would pass the step budget {STEP_BUDGET:,}")
        fn = os.path.join(NEW_DST, str(LV), *map(str, key))
        if first_t is None:
            first_t = time.strftime("%Y-%m-%d %H:%M:%S %z")
        try:
            data = fetch(f"{VOL}/{LV}/{key[0]}/{key[1]}/{key[2]}", expect_len=CHUNK, tag=TAG)
        except urllib.error.HTTPError as exc:
            print(f"L{LV} {key} HTTP {exc.code}: absent")
            absent.append(list(key))
            continue
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        with open(fn, "wb") as f:
            f.write(data)
        got += len(data)
        after, own, _ = ledger()
        if own != got:
            raise SystemExit(f"stop: after {key}, ledger shows {own:,} B tagged {TAG}, counted {got:,} B")
        print(f"L{LV} {key} ok; this run {got:,} B; ledger {after:,} B")
    end, own, lines = ledger()
    rec = {"first_fetch_started": first_t, "finished": time.strftime("%Y-%m-%d %H:%M:%S %z"),
           "ledger_before": start, "ledger_after": end, "lines_before": lines0, "lines_after": lines,
           "fetched_this_run": got, "tag": TAG, "absent": absent, "n_chunks": len(keys) - len(absent),
           "dir": os.path.relpath(NEW_DST, ROOT)}
    json.dump(rec, open(os.path.join(ROOT, "logs", "level0_fetch_record.json"), "w"), indent=1)
    print(f"fetched this run: {got:,} B; ledger after: {end:,} B (was {start:,}); "
          f"left to the limit {LEDGER_LIMIT - end:,} B; absent {absent}")


if __name__ == "__main__":
    main()
