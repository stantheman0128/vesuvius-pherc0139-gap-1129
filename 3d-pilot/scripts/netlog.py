"""HTTP range fetches with a byte ledger (logs/fetch_log.jsonl) and a hard cap.

Every byte read from the network by this task's scripts goes through fetch() so the
total can be checked against the 1.5 GB cap (CAP_BYTES) before each request.
"""
import json
import os
import threading
import time
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOG = os.path.join(ROOT, "logs", "fetch_log.jsonl")
CAP_BYTES = 1_500_000_000
_LOCK = threading.Lock()


def total_fetched():
    if not os.path.exists(LOG):
        return 0
    n = 0
    with open(LOG, "rt") as f:
        for line in f:
            line = line.strip()
            if line:
                n += int(json.loads(line)["bytes"])
    return n


def fetch(url, start=None, end_inclusive=None, expect_len=None, tag="", allow_short=False):
    """GET url (optionally a byte range). Returns bytes. Refuses if the cap would be exceeded."""
    want = None if start is None else (end_inclusive - start + 1)
    already = total_fetched()
    if want is not None and already + want > CAP_BYTES:
        raise RuntimeError(f"cap: {already} + {want} > {CAP_BYTES}")
    req = urllib.request.Request(url, headers={"User-Agent": "pherc0139-smoke/0.1"})
    if start is not None:
        req.add_header("Range", f"bytes={start}-{end_inclusive}")
    t0 = time.time()
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
                status = r.status
            break
        except Exception as exc:  # retry transient errors
            if attempt == 3:
                raise
            time.sleep(2 + 3 * attempt)
    if start is not None and status != 206:
        raise RuntimeError(f"expected 206 for range request, got {status} ({url})")
    if want is not None and len(data) != want and not (allow_short and 0 < len(data) < want):
        raise RuntimeError(f"short read {len(data)} != {want} ({url})")
    if expect_len is not None and len(data) != expect_len:
        raise RuntimeError(f"length {len(data)} != expected {expect_len} ({url})")
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with _LOCK, open(LOG, "at") as f:
        f.write(json.dumps({"t": time.strftime("%Y-%m-%d %H:%M:%S"), "url": url, "start": start,
                            "end": end_inclusive, "bytes": len(data), "tag": tag,
                            "secs": round(time.time() - t0, 2)}) + "\n")
    return data
