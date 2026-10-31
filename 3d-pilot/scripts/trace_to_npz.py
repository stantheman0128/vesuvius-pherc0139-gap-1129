"""Convert a sheet-tracer output folder (LZW tifs) to trace.npz (x, y, z, conf) for analysis on Windows Python,
which has no imagecodecs. Run in the WSL venv with EXTRA_SITE on the path. Usage: python -I trace_to_npz.py <dir>"""
import os
import sys

if os.environ.get("EXTRA_SITE"):
    sys.path.append(os.environ["EXTRA_SITE"])
import numpy as np  # noqa: E402
import tifffile  # noqa: E402

d = sys.argv[1]
arrs = {c: tifffile.imread(os.path.join(d, f"{c}.tif")) for c in "xyz"}
np.savez(os.path.join(d, "trace.npz"), conf=np.load(os.path.join(d, "confidence.npy")), **arrs)
print(d, {c: a.shape for c, a in arrs.items()})
