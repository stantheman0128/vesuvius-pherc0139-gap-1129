"""Read the layout of a remote uncompressed float32 TIFF and fetch row bands by Range requests.

layout(url) -> dict(width, height, offsets per strip, rows_per_strip, ...)
rows(url, layout, r0, r1) -> float32 array (r1 - r0, width)
"""
import json
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import fetch  # noqa: E402


class _Spliced:
    """Bytes view of [0, len(head)) plus [off, off + len(tail)) at their file offsets."""

    def __init__(self, head, off, tail):
        self.head, self.off, self.tail = head, off, tail

    def __len__(self):
        return self.off + len(self.tail)

    def __getitem__(self, s):
        if s.start >= self.off:
            return self.tail[s.start - self.off:s.stop - self.off]
        if s.stop <= len(self.head):
            return self.head[s]
        raise IndexError("slice spans the unfetched middle of the file")


def layout(url, cache_json=None):
    if cache_json and os.path.exists(cache_json):
        return json.load(open(cache_json))
    head = fetch(url, 0, 65535, tag="tif-header")
    if head[:2] != b"II":
        raise RuntimeError("not little-endian TIFF")
    magic = struct.unpack("<H", head[2:4])[0]
    if magic != 42:
        raise RuntimeError(f"unsupported TIFF magic {magic}")
    ifd = struct.unpack("<I", head[4:8])[0]
    if ifd + 4096 > len(head):
        # IFD stored after the pixel data: fetch it and address it as if it were at its file offset
        tail = fetch(url, ifd, ifd + 4095, tag="tif-ifd", allow_short=True)
        head = _Spliced(head, ifd, tail)
    n = struct.unpack("<H", head[ifd:ifd + 2])[0]
    tags = {}
    for i in range(n):
        tag, typ, cnt, val = struct.unpack("<HHII", head[ifd + 2 + 12 * i: ifd + 14 + 12 * i])
        tags[tag] = (typ, cnt, val, head[ifd + 2 + 12 * i + 8: ifd + 14 + 12 * i])
    def scalar(tag):
        typ, cnt, val, raw = tags[tag]
        return struct.unpack("<H", raw[:2])[0] if typ == 3 else val
    def array(tag):
        typ, cnt, val, raw = tags[tag]
        fmt = {3: "H", 4: "I"}[typ]
        size = struct.calcsize(fmt) * cnt
        if size <= 4:
            return list(struct.unpack("<" + fmt * cnt, raw[:size]))
        try:
            if val + size > len(head):
                raise IndexError
            return list(struct.unpack("<" + fmt * cnt, head[val:val + size]))
        except IndexError:
            more = fetch(url, val, val + size - 1, tag="tif-ifd-array")
            return list(struct.unpack("<" + fmt * cnt, more))
    W, H = scalar(256), scalar(257)
    comp = scalar(259) if 259 in tags else 1
    bps = scalar(258)
    fmt = scalar(339) if 339 in tags else 1
    if comp != 1 or bps != 32 or fmt != 3:
        raise RuntimeError(f"need uncompressed float32, got compression={comp} bps={bps} format={fmt}")
    if 322 in tags:
        raise RuntimeError("tiled TIFF not supported")
    rps = scalar(278) if 278 in tags else H
    offs = array(273)
    lay = {"url": url, "width": W, "height": H, "rows_per_strip": rps, "offsets": offs}
    if cache_json:
        json.dump(lay, open(cache_json, "wt"))
    return lay


def rows(lay, r0, r1):
    W, rps = lay["width"], lay["rows_per_strip"]
    out = np.empty((r1 - r0, W), np.float32)
    r = r0
    while r < r1:
        s = r // rps
        within = r - s * rps
        n = min(r1 - r, rps - within)
        start = lay["offsets"][s] + within * W * 4
        buf = fetch(lay["url"], start, start + n * W * 4 - 1, tag="tif-rows")
        out[r - r0:r - r0 + n] = np.frombuffer(buf, "<f4").reshape(n, W)
        r += n
    return out
