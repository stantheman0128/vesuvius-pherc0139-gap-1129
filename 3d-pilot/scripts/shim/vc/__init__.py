"""Stand-in for the volume-cartographer `vc` Python bindings (a compiled C++ module that is not built on this
machine). It implements only what infer_rowcol_triplet_wraps.py calls on vc.Volume: open, shape_at, chunk_shape,
dtype, has_scale_level, present_scale_levels, read_zyx(..., missing_policy="all_fill"), and
set_chunk_cache_budget. It reads a local zarr v2 pyramid of uncompressed chunks (<root>/<level>/.zarray and
<root>/<level>/<z>/<y>/<x>) with numpy.

Chunks that the tool asks for but that are absent on disk are filled with 0 (like missing_policy="all_fill") and
recorded in MISSING; the launcher reports them, so a CT box that is too small cannot pass silently.
"""
import json
import os

import numpy as np

MISSING = set()
READ = set()
_CACHE_BUDGET = [0]


def set_chunk_cache_budget(nbytes):
    _CACHE_BUDGET[0] = int(nbytes)


class Volume:
    def __init__(self, root):
        self.root = str(root)
        self._levels = {}
        for name in sorted(os.listdir(self.root)):
            p = os.path.join(self.root, name, ".zarray")
            if name.isdigit() and os.path.exists(p):
                meta = json.load(open(p))
                if meta.get("compressor") is not None or meta.get("filters"):
                    raise RuntimeError("shim vc reads uncompressed zarr chunks only")
                self._levels[int(name)] = meta
        if not self._levels:
            raise RuntimeError(f"no zarr levels under {self.root}")
        dts = {np.dtype(m["dtype"]) for m in self._levels.values()}
        if len(dts) != 1:
            raise RuntimeError("mixed dtypes across levels")
        self.dtype = dts.pop()

    @classmethod
    def open(cls, path):
        return cls(path)

    @classmethod
    def open_url(cls, url):
        raise RuntimeError("shim vc: remote volumes are not supported; pass a local directory")

    def has_scale_level(self, level):
        return int(level) in self._levels

    def present_scale_levels(self):
        return sorted(self._levels)

    def shape_at(self, level):
        return tuple(int(v) for v in self._levels[int(level)]["shape"])

    def chunk_shape(self, level):
        return tuple(int(v) for v in self._levels[int(level)]["chunks"])

    def read_zyx(self, starts, shape, level=0, missing_policy="all_fill"):
        meta = self._levels[int(level)]
        vshape = np.array(meta["shape"], np.int64)
        ch = np.array(meta["chunks"], np.int64)
        sep = meta.get("dimension_separator", ".")
        fill = meta.get("fill_value") or 0
        start = np.array(starts, np.int64)
        size = np.array(shape, np.int64)
        out = np.full(tuple(size), fill, dtype=self.dtype)
        lo = np.maximum(start, 0)
        hi = np.minimum(start + size, vshape)
        if np.any(hi <= lo):
            return out
        c0, c1 = lo // ch, (hi - 1) // ch
        for cz in range(c0[0], c1[0] + 1):
            for cy in range(c0[1], c1[1] + 1):
                for cx in range(c0[2], c1[2] + 1):
                    key = (int(level), int(cz), int(cy), int(cx))
                    fn = os.path.join(self.root, str(level), sep.join(str(k) for k in key[1:]))
                    if not os.path.exists(fn):
                        MISSING.add(key)
                        continue
                    READ.add(key)
                    arr = np.memmap(fn, dtype=self.dtype, mode="r", shape=tuple(ch))
                    base = np.array(key[1:]) * ch
                    a, b = np.maximum(lo, base), np.minimum(hi, base + ch)
                    out[tuple(slice(a[i] - start[i], b[i] - start[i]) for i in range(3))] = arr[
                        tuple(slice(a[i] - base[i], b[i] - base[i]) for i in range(3))]
        return out
