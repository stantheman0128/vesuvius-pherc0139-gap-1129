"""Run villa's sheet tracer (tracers/sheets/tracing/inference.py, origin/main a23ec7c) unmodified on CPU through
villa_env (4 threads, weights_only=True). Records wall time, peak RSS, and CT chunks the tracer asked for that are
not on disk (they read as 0, so a too-small box would otherwise pass silently).
Usage (WSL venv): python -I run_sheet.py <record_json> -- <inference.py args...>
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # python -I does not add the script dir
# The tracer writes LZW tifs, which tifffile can only do with imagecodecs; villa's models extra declares it but the
# venv we ran in lacks it. EXTRA_SITE, if set, points at a separate folder holding imagecodecs (we used 2026.10.10);
# nothing new is installed here. Must precede the tifffile import.
if os.environ.get("EXTRA_SITE"):
    sys.path.append(os.environ["EXTRA_SITE"])

import villa_env  # noqa: E402

import torch  # noqa: E402
from vesuvius.neural_tracing.tracers.shared.crop_sampling import volume as V  # noqa: E402
from vesuvius.neural_tracing.tracers.sheets.tracing import inference as T  # noqa: E402

MISSING, READ = set(), set()
_orig = V.ChunkedArray._load


def _load(self, key):
    arr = _orig(self, key)
    (READ if arr is not None else MISSING).add(tuple(int(k) for k in key))
    return arr


V.ChunkedArray._load = _load


def main():
    i = sys.argv.index("--")
    rec_path, argv = sys.argv[1], sys.argv[i + 1:]
    t0 = time.time()
    status, err = "ok", None
    try:
        T.main(argv)
    except BaseException as exc:
        status, err = "error", f"{type(exc).__name__}: {exc}"
        raise
    finally:
        rec = {"status": status, "error": err, "wall_s": round(time.time() - t0, 2),
               "peak_rss_mb": round(villa_env.peak_rss_mb(), 1), "torch": torch.__version__,
               "threads": torch.get_num_threads(), "torch_load_calls": villa_env.LOAD_CALLS,
               "ct_chunks_read": sorted(READ), "ct_chunks_missing": sorted(MISSING), "argv": argv}
        os.makedirs(os.path.dirname(os.path.abspath(rec_path)), exist_ok=True)
        json.dump(rec, open(rec_path, "wt"), indent=1)
        print(json.dumps({k: rec[k] for k in ("status", "error", "wall_s", "peak_rss_mb", "threads")}))
        print("chunks read", len(READ), "missing", sorted(MISSING))


if __name__ == "__main__":
    main()
