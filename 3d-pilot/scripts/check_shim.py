"""The vc shim against villa's own ChunkedArray: the copy model's crop box read both ways must be identical.
Usage (WSL venv): python -I check_shim.py"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import villa_env  # noqa: E402
import numpy as np  # noqa: E402
import vc  # noqa: E402
from vesuvius.neural_tracing.tracers.shared.crop_sampling.volume import ChunkedArray  # noqa: E402

root = os.path.join(villa_env.ROOT, "downloads", "ct", "vol.zarr")
plan = json.load(open(os.path.join(villa_env.ROOT, "out", "plan_ct_l2.json")))
(mn, mx) = plan["boxes"][0]
size = [b - a for a, b in zip(mn, mx)]
a = vc.Volume.open(root).read_zyx(mn, size, level=2)
b = ChunkedArray(os.path.join(root, "2")).read(np.array(mn), np.array(size))
print("shape", a.shape, "equal", bool(np.array_equal(a, b)), "mean", float(a.mean()), "nonzero frac",
      float((a > 0).mean()))
