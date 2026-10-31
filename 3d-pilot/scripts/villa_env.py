"""Shared setup for running villa origin/main code in the WSL venv (imported first by the run_* scripts).

- Caps numeric threads at 4 (env vars before torch is imported, then torch.set_num_threads(4)).
- Puts the villa checkout named by VILLA_SRC (we used origin/main a23ec7c) first on sys.path, ahead of the venv's
  editable install of another villa checkout, and the vc / vc_object_store shims after it.
- Forces torch.load(weights_only=True) for every call, whatever the caller passes, and counts the calls.
"""
import os
import sys

for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "NUMBA_NUM_THREADS"):
    os.environ[_k] = "4"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WT_SRC = os.environ["VILLA_SRC"]  # vesuvius/src of a villa checkout at origin/main a23ec7c
sys.path.insert(0, os.path.join(HERE, "shim"))
sys.path.insert(0, WT_SRC)

import torch  # noqa: E402

torch.set_num_threads(4)
_orig_load = torch.load
LOAD_CALLS = []


def _load_weights_only(*a, **k):
    LOAD_CALLS.append({"requested_weights_only": k.get("weights_only", "default")})
    k["weights_only"] = True
    return _orig_load(*a, **k)


torch.load = _load_weights_only

import vesuvius  # noqa: E402

assert os.path.realpath(vesuvius.__file__).startswith(os.path.realpath(WT_SRC)), vesuvius.__file__


def peak_rss_mb():
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
