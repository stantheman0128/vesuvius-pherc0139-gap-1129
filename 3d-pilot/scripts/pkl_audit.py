"""Audit a torch data.pkl without executing it.

1. Lists every GLOBAL / STACK_GLOBAL the pickle references (pickletools, no execution).
2. Unpickles with a find_class that returns inert stubs only, so no foreign code runs, and maps
   each tensor storage key to the top-level dict path that holds it.
Usage: python -I pkl_audit.py <data.pkl> <index_json> <member_prefix> <out_json>
"""
import collections
import json
import pickle
import pickletools
import sys


def list_globals(raw):
    found = collections.Counter()
    strings = []
    for op, arg, _ in pickletools.genops(raw):
        if op.name in ("SHORT_BINUNICODE", "BINUNICODE", "UNICODE", "BINUNICODE8"):
            strings.append(arg)
        if op.name == "GLOBAL":
            found[arg.replace(" ", ".")] += 1
        elif op.name == "STACK_GLOBAL":
            found[f"{strings[-2]}.{strings[-1]}"] += 1
    return found


class Stub:
    def __init__(self, qual, args=()):
        self.qual, self.args = qual, args

    def __repr__(self):
        return f"<{self.qual}>"

    def __setstate__(self, st):
        self.state = st


class TensorStub:
    def __init__(self, storage, offset, size, stride, *rest):
        self.storage, self.offset, self.size, self.stride = storage, offset, tuple(size), tuple(stride)


class R(pickle.Unpickler):
    def find_class(self, module, name):
        q = f"{module}.{name}"
        if q == "collections.OrderedDict":
            return collections.OrderedDict
        if q in ("torch._utils._rebuild_tensor_v2", "torch._utils._rebuild_tensor"):
            return lambda *a: TensorStub(*a[:4])
        if q == "torch._utils._rebuild_parameter":
            return lambda data, *a: data
        def mk(*a, **k):
            return Stub(q, a)
        mk.__qualname__ = q
        if q.endswith("Storage") and module == "torch":
            return Stub(q)
        return mk

    def persistent_load(self, pid):
        # ('storage', storage_type, key, location, numel)
        return {"key": str(pid[2]), "dtype": repr(pid[1]), "numel": int(pid[4])}


def walk(obj, path, out, other):
    if isinstance(obj, TensorStub):
        out.append((path, obj.storage["key"], obj.storage["dtype"], obj.size))
    elif isinstance(obj, dict):
        for k, v in obj.items():
            walk(v, path + [str(k)], out, other)
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            walk(v, path + [str(i)], out, other)
    elif isinstance(obj, Stub):
        other.append((path, obj.qual))
        st = getattr(obj, "state", None)
        if st is not None:
            walk(st, path + ["<state>"], out, other)
        for i, v in enumerate(obj.args):
            walk(v, path + [f"<arg{i}>"], out, other)


def main():
    pkl, index_json, prefix, out_json = sys.argv[1:5]
    raw = open(pkl, "rb").read()
    g = list_globals(raw)
    print("globals referenced:")
    for k, v in sorted(g.items()):
        print(f"  {k}: {v}")
    obj = R(open(pkl, "rb")).load()
    print("top-level type:", type(obj).__name__, "keys:", list(obj.keys()) if isinstance(obj, dict) else None)
    tensors, other = [], []
    walk(obj, [], tensors, other)
    idx = json.load(open(index_json))
    size_of = {e["name"]: e["usize"] for e in idx["entries"]}
    by_top = collections.defaultdict(lambda: [0, set()])
    for path, key, dt, shape in tensors:
        top = path[0] if path else ""
        by_top[top][1].add(key)
    res = {}
    for top, (_, keys) in by_top.items():
        b = sum(size_of.get(f"{prefix}/data/{k}", 0) for k in keys)
        res[top] = {"n_storages": len(keys), "bytes": b}
        print(f"top key {top!r}: {len(keys)} storages, {b} bytes")
    print("non-tensor objects (first 30):")
    for p, q in other[:30]:
        print("  ", "/".join(p), q)
    cfg = obj.get("config") if isinstance(obj, dict) else None
    model_keys = sorted(by_top.get("model", [0, set()])[1], key=int)
    json.dump({"globals": dict(g), "by_top": res, "model_storage_keys": model_keys,
               "config": cfg if isinstance(cfg, dict) else repr(cfg),
               "non_tensor_objects": [["/".join(p), q] for p, q in other]},
              open(out_json, "wt"), indent=1, default=repr)


if __name__ == "__main__":
    main()
