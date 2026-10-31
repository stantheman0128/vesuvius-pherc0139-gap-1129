"""Fetch named members of a remote zip (stored, uncompressed) by Range requests and check CRC32.

Usage: python -I hf_zip_member.py <index_json> <out_dir> <member> [<member> ...]
       python -I hf_zip_member.py <index_json> <out_dir> @<list_file>
Each member is written to <out_dir>/<member path>. Members whose file already exists with the
right size and CRC are skipped. Adjacent members are fetched in merged ranges (<= 64 MB each).
"""
import json
import os
import struct
import sys
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import fetch  # noqa: E402

MAX_SPAN = 64 * 1024 * 1024


def main():
    index_json, out_dir = sys.argv[1], sys.argv[2]
    names = sys.argv[3:]
    if len(names) == 1 and names[0].startswith("@"):
        with open(names[0][1:], "rt") as f:
            names = [ln.strip() for ln in f if ln.strip()]
    with open(index_json, "rt") as f:
        idx = json.load(f)
    by_name = {e["name"]: e for e in idx["entries"]}
    todo = []
    for n in names:
        e = by_name[n]
        if e["comp"] != 0:
            raise SystemExit(f"{n}: compressed member (method {e['comp']}) not supported")
        dst = os.path.join(out_dir, *n.split("/"))
        if os.path.exists(dst) and os.path.getsize(dst) == e["usize"]:
            with open(dst, "rb") as f:
                if zlib.crc32(f.read()) & 0xFFFFFFFF == e["crc"]:
                    continue
        todo.append(e)
    todo.sort(key=lambda e: e["lho"])
    # group into spans; local header is 30 + name + extra (<= 30 + 256 + 128 assumed, checked below)
    groups, cur = [], []
    for e in todo:
        end = e["lho"] + 30 + 512 + e["csize"]
        if cur and (end - cur[0]["lho"] > MAX_SPAN or e["lho"] - (cur[-1]["lho"] + 30 + 512 + cur[-1]["csize"]) > 4096):
            groups.append(cur)
            cur = []
        cur.append(e)
    if cur:
        groups.append(cur)
    nbytes = 0
    for g in groups:
        s = g[0]["lho"]
        last = g[-1]
        e_end = min(idx["size"] - 1, last["lho"] + 30 + 512 + last["csize"])
        buf = fetch(idx["url"], s, e_end, tag="zip-member")
        nbytes += len(buf)
        for e in g:
            o = e["lho"] - s
            sig, _, _, _, _, _, _, _, _, nlen, xlen = struct.unpack("<IHHHHHIIIHH", buf[o:o + 30])
            if sig != 0x04034B50:
                raise SystemExit(f"bad local header for {e['name']}")
            d0 = o + 30 + nlen + xlen
            data = buf[d0:d0 + e["csize"]]
            if len(data) != e["csize"]:
                raise SystemExit(f"{e['name']}: span too short (header {nlen}+{xlen})")
            crc = zlib.crc32(data) & 0xFFFFFFFF
            if crc != e["crc"]:
                raise SystemExit(f"{e['name']}: CRC mismatch {crc:08x} != {e['crc']:08x}")
            dst = os.path.join(out_dir, *e["name"].split("/"))
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            with open(dst, "wb") as f:
                f.write(data)
    print(f"fetched {len(todo)} members in {len(groups)} ranges, {nbytes} bytes; skipped {len(names) - len(todo)}")


if __name__ == "__main__":
    main()
