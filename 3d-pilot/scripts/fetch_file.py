"""Download one whole file through the byte ledger and check its sha256.
Usage: python -I fetch_file.py <url> <dest> <expected_size> <expected_sha256>
"""
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import fetch  # noqa: E402

url, dest, size, sha = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
if not (os.path.exists(dest) and os.path.getsize(dest) == size):
    # ranged in 64 MB pieces so a dropped connection costs one piece
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest + ".part", "wb") as f:
        for s in range(0, size, 64 << 20):
            f.write(fetch(url, s, min(size, s + (64 << 20)) - 1, tag="file"))
    os.replace(dest + ".part", dest)
h = hashlib.sha256()
with open(dest, "rb") as f:
    for blk in iter(lambda: f.read(1 << 22), b""):
        h.update(blk)
print(os.path.getsize(dest), h.hexdigest(), "OK" if h.hexdigest() == sha else "SHA256 MISMATCH")
