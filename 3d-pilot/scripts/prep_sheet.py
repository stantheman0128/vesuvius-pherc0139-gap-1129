"""Inputs for the sheet tracer smoke test at C240.

The checkpoint traces at its volume_level 0 in 9.6 um voxels. The 1.129 um scan's level 3 is 9.032 um, so the
tracer is given a local zarr whose level "0" holds level-3 chunks (downloads/ct/vol_level3_as_level0.zarr), and
all tracer coordinates are level-3 voxels (level-0 / 8).

Writes out/inputs/w035_C240_seed_l3/{x,y,z}.tif (the w035 band in level-3 voxels; the tracer reads only these
three tifs), out/inputs/sheet_plan.json (seed row/cols, level-3 chunks to fetch).
Usage: python -I prep_sheet.py [seed_row_offset=-30] [half_cols=48] [steps=8]
"""
import json
import math
import os
import site
import sys

sys.path.append(site.getusersitepackages())
import numpy as np  # noqa: E402
import tifffile  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from netlog import ROOT  # noqa: E402

off = int(sys.argv[1]) if len(sys.argv) > 1 else -30
half = int(sys.argv[2]) if len(sys.argv) > 2 else 48
steps = int(sys.argv[3]) if len(sys.argv) > 3 else 8
STEP_L3, PAD_L3 = 20.0, 64.0   # tracer step (checkpoint data.step) and a generous read pad (view radius + step + 4)

pj = json.load(open(os.path.join(ROOT, "out", "inputs", "w035_C240_patch.json")))
d = np.load(os.path.join(ROOT, "downloads", "mesh_bands", "w035_band.npz"))
g, row0 = d["g"], int(d["row0"])
srow = pj["center_row"] + off - row0          # index into the band
c0, c1 = pj["center_col"] - half, pj["center_col"] + half
line = g[srow, c0:c1]
behind = g[srow - 1, c0:c1]
assert np.all(line != -1) and np.all(behind != -1), "seed line or the line behind it has invalid vertices"
g3 = np.where(np.all(g != -1, -1, keepdims=True), g / 8.0, -1.0).astype(np.float32)
out = os.path.join(ROOT, "out", "inputs", "w035_C240_seed_l3")
os.makedirs(out, exist_ok=True)
for k, c in enumerate("xyz"):
    tifffile.imwrite(os.path.join(out, f"{c}.tif"), np.ascontiguousarray(g3[..., k]))

l3 = line / 8.0
adv = (line - behind).mean(0)
adv = adv / np.linalg.norm(adv)
end = l3 + adv[None] * STEP_L3 * steps
pts = np.concatenate([l3, end])
lo = pts.min(0) - PAD_L3
hi = pts.max(0) + PAD_L3
loc = json.load(open(os.path.join(ROOT, "out", "locate.json")))
u = np.array(loc["umb11_xyz"]) / 8.0
chunks = []
for z in range(int(lo[2] // 128), int(hi[2] // 128) + 1):
    for y in range(int(lo[1] // 128), int(hi[1] // 128) + 1):
        for x in range(int(lo[0] // 128), int(hi[0] // 128) + 1):
            chunks.append([z, y, x])
seed_ang = [math.degrees(math.atan2(p[1] - u[1], p[0] - u[0])) % 360 for p in l3[[0, -1]]]
plan = {"level": 3, "seed_band_row": int(srow), "seed_tif_line": f"row:{srow}", "seed_range": f"{c0}:{c1}",
        "seed_row_global": int(srow + row0), "advance_dir_xyz": adv.tolist(), "steps": steps,
        "seed_len_l3": float(np.linalg.norm(np.diff(l3, axis=0), axis=1).sum()),
        "seed_angles_deg": seed_ang, "read_lo_xyz_l3": lo.tolist(), "read_hi_xyz_l3": hi.tolist(),
        "chunks": chunks}
json.dump(plan, open(os.path.join(ROOT, "out", "inputs", "sheet_plan.json"), "wt"), indent=1)
print(f"seed: band row {srow} (global {srow + row0}), cols {c0}:{c1}, length {plan['seed_len_l3']:.1f} L3 vox, "
      f"angles {seed_ang[0]:.1f}..{seed_ang[1]:.1f} deg, advance dir xyz {np.round(adv, 3)}")
print(f"level-3 chunks: {len(chunks)} ({len(chunks) * 2097152} bytes), lo {np.round(lo)} hi {np.round(hi)}")
