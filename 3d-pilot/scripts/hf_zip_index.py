"""List the members of a torch zip checkpoint on Hugging Face without downloading it.

Usage: python -I hf_zip_index.py <url> <total_size> <out_json> [tail_bytes]
Reads only the tail of the file (EOCD + central directory) with a Range request.
"""
import json
import os
import struct
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import fetch  # noqa: E402


def main():
    url, size, out_json = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    tail = int(sys.argv[4]) if len(sys.argv) > 4 else 262144
    start = max(0, size - tail)
    buf = fetch(url, start, size - 1, tag="zip-tail")
    eocd = buf.rfind(b"PK\x05\x06")
    if eocd < 0:
        raise SystemExit("no EOCD in tail")
    (_, _, _, n_this, n_total, cd_size, cd_off, _) = struct.unpack("<IHHHHIIH", buf[eocd:eocd + 22])
    if cd_off == 0xFFFFFFFF or n_total == 0xFFFF or cd_size == 0xFFFFFFFF:
        loc = buf.rfind(b"PK\x06\x07", 0, eocd)
        (_, _, z64_off, _) = struct.unpack("<IIQI", buf[loc:loc + 20])
        z = z64_off - start
        (_, _, _, _, _, _, n_this, n_total, cd_size, cd_off) = struct.unpack("<IQHHIIQQQQ", buf[z:z + 56])
    if cd_off < start:
        raise SystemExit(f"central directory starts before tail: need tail >= {size - cd_off}")
    p = cd_off - start
    entries = []
    for _ in range(n_total):
        sig, _, _, flags, comp, _, _, crc, csize, usize, nlen, xlen, clen, _, _, _, lho = struct.unpack(
            "<IHHHHHHIIIHHHHHII", buf[p:p + 46])
        if sig != 0x02014B50:
            raise SystemExit(f"bad central dir signature at {p}")
        name = buf[p + 46:p + 46 + nlen].decode("utf-8")
        extra = buf[p + 46 + nlen:p + 46 + nlen + xlen]
        # zip64 extra field
        q = 0
        while q + 4 <= len(extra):
            hid, hlen = struct.unpack("<HH", extra[q:q + 4])
            if hid == 0x0001:
                vals = list(struct.unpack("<" + "Q" * (hlen // 8), extra[q + 4:q + 4 + (hlen // 8) * 8]))
                if usize == 0xFFFFFFFF:
                    usize = vals.pop(0)
                if csize == 0xFFFFFFFF:
                    csize = vals.pop(0)
                if lho == 0xFFFFFFFF:
                    lho = vals.pop(0)
            q += 4 + hlen
        entries.append({"name": name, "comp": comp, "csize": csize, "usize": usize, "lho": lho, "crc": crc})
        p += 46 + nlen + xlen + clen
    with open(out_json, "wt") as f:
        json.dump({"url": url, "size": size, "cd_off": cd_off, "n": len(entries), "entries": entries}, f, indent=1)
    tot = sum(e["usize"] for e in entries)
    print(f"entries={len(entries)} sum_usize={tot} cd_off={cd_off}")
    big = sorted(entries, key=lambda e: -e["usize"])[:8]
    for e in big:
        print(e["name"], e["usize"])
    for e in entries:
        if not e["name"].split("/")[-2:-1] == ["data"]:
            print("non-data:", e["name"], e["usize"])


if __name__ == "__main__":
    main()
