"""Run villa's infer_rowcol_triplet_wraps (origin/main a23ec7c) unmodified on CPU, through villa_env
(4 threads, weights_only=True, vc shim). Writes <out_dir>/run_record.json with wall time, peak RSS,
the CT chunks read and any chunks the tool asked for that were not downloaded.
Usage (WSL venv): python -I run_copy.py <out_dir> -- <tool args...>
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # python -I does not add the script dir
import villa_env  # noqa: E402

import torch  # noqa: E402
import vc  # noqa: E402  (the shim)
from vesuvius.neural_tracing.inference import infer_rowcol_triplet_wraps as m  # noqa: E402


def main():
    i = sys.argv.index("--")
    out_dir, argv = sys.argv[1], sys.argv[i + 1:]
    os.makedirs(out_dir, exist_ok=True)
    args = m.parse_args(argv)
    t0 = time.time()
    status, err, outputs = "ok", None, None
    try:
        outputs = m.run(args)
    except BaseException as exc:  # record, then re-raise
        status, err = "error", f"{type(exc).__name__}: {exc}"
        raise
    finally:
        rec = {
            "status": status, "error": err, "wall_s": round(time.time() - t0, 2),
            "peak_rss_mb": round(villa_env.peak_rss_mb(), 1),
            "torch": torch.__version__, "threads": torch.get_num_threads(),
            "omp_env": os.environ.get("OMP_NUM_THREADS"),
            "torch_load_calls": villa_env.LOAD_CALLS,
            "ct_chunks_read": sorted(vc.READ), "ct_chunks_missing": sorted(vc.MISSING),
            "argv": argv, "outputs": outputs,
        }
        json.dump(rec, open(os.path.join(out_dir, "run_record.json"), "wt"), indent=1, default=str)
        print(json.dumps({k: rec[k] for k in ("status", "error", "wall_s", "peak_rss_mb", "threads")}))
        print("chunks read", len(rec["ct_chunks_read"]), "missing", len(rec["ct_chunks_missing"]))


if __name__ == "__main__":
    main()
