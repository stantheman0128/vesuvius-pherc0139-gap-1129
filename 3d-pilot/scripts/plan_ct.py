"""Which level-2 CT chunks will infer_rowcol_triplet_wraps read for a tifxyz input? Uses the tool's own input
loader and cover-bbox generator (same arguments as its defaults), so the download matches what it reads.
Usage (WSL venv): python -I plan_ct.py <tifxyz_dir> <level> <out_json> [crop_z crop_y crop_x] [overlap]
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # python -I does not add the script dir
import villa_env  # noqa: E402,F401  (paths, threads)
import numpy as np
from vesuvius.neural_tracing.inference import infer_rowcol_triplet_wraps as m


def main():
    tif, level, out_json = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    crop = tuple(int(v) for v in sys.argv[4:7]) if len(sys.argv) >= 7 else (128, 384, 384)
    overlap = float(sys.argv[7]) if len(sys.argv) >= 8 else 0.5
    surface, grid, valid = m._load_input_grid(tif, retarget_factor=float(2 ** level))
    world, uv = m._extract_surface_points_for_iteration(grid, valid, float(2 ** level), False)
    recs = m._generate_cover_bboxes_from_points(world, "patch", crop, overlap=overlap, prune_bboxes=True,
                                                band_workers=1)
    chunks = set()
    boxes = []
    for r in recs:
        mn, _ = m._bbox_to_min_corner_and_bounds_array(tuple(r["bbox"]))
        mx = mn + np.array(crop)
        boxes.append([mn.tolist(), mx.tolist()])
        lo = np.maximum(mn, 0) // 128
        hi = (mx - 1) // 128
        for z in range(lo[0], hi[0] + 1):
            for y in range(lo[1], hi[1] + 1):
                for x in range(lo[2], hi[2] + 1):
                    chunks.add((int(z), int(y), int(x)))
    print(f"points {len(world)}, crops {len(recs)}, chunks {len(chunks)} ({len(chunks) * 2097152} bytes)")
    for b in boxes:
        print("  crop zyx min", b[0], "max", b[1])
    json.dump({"level": level, "crop": crop, "overlap": overlap, "boxes": boxes, "chunks": sorted(chunks),
               "world_min_zyx": world.min(0).tolist(), "world_max_zyx": world.max(0).tolist()},
              open(out_json, "wt"), indent=1)


if __name__ == "__main__":
    main()
