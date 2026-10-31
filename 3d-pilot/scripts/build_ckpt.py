"""Rebuild a model-only checkpoint from the fetched members of copy_displacement_latest.pth.

The full file is 1,698,235,447 bytes; 1,131,873,288 of them are AdamW optimizer state that inference never reads,
so only data.pkl and the 310 storages referenced under 'model' were fetched (each CRC32-checked against the zip
central directory). This script unpickles data.pkl with an allowlist of exactly the globals it references
(collections.OrderedDict, torch.FloatStorage, torch._utils._rebuild_tensor_v2; pkl_audit.py lists them and all
three are on torch's weights_only allowlist), fills the model tensors from the local storage files, and saves
{'model', 'config', 'step', 'wandb_run_id'} with torch.save. Optimizer tensors become None.
Usage (WSL venv): python -I build_ckpt.py <members_dir>/ckpt_332500 <out.pth>
"""
import collections
import os
import pickle
import sys

import torch


def main():
    src, out = sys.argv[1], sys.argv[2]
    data_dir = os.path.join(src, "data")
    loaded, skipped = [0], [0]
    cache = {}  # one storage object per key, so tensors that share a storage keep sharing it

    class U(pickle.Unpickler):
        def find_class(self, module, name):
            q = f"{module}.{name}"
            if q == "collections.OrderedDict":
                return collections.OrderedDict
            if q == "torch.FloatStorage":
                return torch.float32
            if q == "torch._utils._rebuild_tensor_v2":
                def rebuild(storage, offset, size, stride, requires_grad, backward_hooks, *rest):
                    if storage is None:
                        return None
                    t = torch.empty(0, dtype=torch.float32)
                    t.set_(storage, offset, torch.Size(size), tuple(stride))
                    return t
                return rebuild
            raise pickle.UnpicklingError(f"global not allowed: {q}")

        def persistent_load(self, pid):
            kind, dtype, key, location, numel = pid
            if kind != "storage" or dtype is not torch.float32:
                raise pickle.UnpicklingError(f"unexpected persistent id {pid!r}")
            if key in cache:
                return cache[key]
            fn = os.path.join(data_dir, str(key))
            if not os.path.exists(fn):
                skipped[0] += 1
                return None
            raw = open(fn, "rb").read()
            if len(raw) != numel * 4:
                raise ValueError(f"storage {key}: {len(raw)} bytes != numel {numel} * 4")
            loaded[0] += 1
            cache[key] = torch.frombuffer(bytearray(raw), dtype=torch.float32).untyped_storage()
            return cache[key]

    with open(os.path.join(src, "data.pkl"), "rb") as f:
        obj = U(f).load()
    model = obj["model"]
    bad = [k for k, v in model.items() if v is None]
    if bad:
        raise SystemExit(f"{len(bad)} model tensors missing, e.g. {bad[:3]}")
    n_params = sum(v.numel() for v in model.values())
    print(f"storages loaded {loaded[0]}, skipped (optimizer) {skipped[0]}; model tensors {len(model)}, "
          f"elements {n_params}")
    torch.save({"model": model, "config": obj["config"], "step": obj.get("step"),
                "wandb_run_id": obj.get("wandb_run_id")}, out)
    print("saved", out, os.path.getsize(out))


if __name__ == "__main__":
    main()
