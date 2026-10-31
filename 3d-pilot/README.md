# 3D pilot: can a sheet tracer follow the layers into the w035 gap?

This folder follows up on [ScrollPrize/villa#2167](https://github.com/ScrollPrize/villa/issues/2167). The main README of this repository reads the w034/w036 gap in the 1.129 µm scan one radial cut at a time. It cannot say which winding each band is. The thread proposed following the layers into the gap in 3D, from where they are still separate. This folder tests whether villa's sheet tracer can do that on this scan, before using it inside the gap.

Our answer is that it probably cannot, at the pitch the gap has. Where the published meshes are wider apart the tracer follows the CT well. Where two sheets come closer than about one sheet's width, the CT shows a single ridge, even at levels 2, 1 and 0 (2x, 4x and 8x finer than the level the tracer reads; level 0 is the scan's full 1.129 µm). A trace that starts on such a ridge has to choose a sheet when the two part, and the CT gives it nothing to choose by. In the gap the pitch is about 63 µm and several places already show w034 and w036 in one dense band. So we did not run the tracer at any gap place, and we did not run the 12 controls that were planned for that run.

No gap place was opened, downloaded or traced. Every CT chunk used here has all 8 corners at least 30° outside the gap sector (327° to 85°).

Frame: heights z are full-resolution voxels of volume `20260102150214` (2.399 µm), angles are degrees about the umbilicus at the same z, as in the main README. Tracer coordinates are level-3 voxels of the 1.129 µm scan `20260413113053`.

## What was run

The tool is villa's sheet tracer, `vesuvius/src/vesuvius/neural_tracing/tracers/sheets/tracing/inference.py` at villa origin/main `a23ec7c` (2026-10-30), unmodified. Checkpoint: Hugging Face `scrollprize/sheet_tracer_7pt_gap60`, commit `443497f`, file `last.pt`, sha256 `b1d133bcb8bf3595fb9f3676baca82ad64acae3010a55a65b1dda25b2cc853ec`.

Settings were the same at every place:

- Level 3 of the 1.129 µm zarr, 9.03 µm voxels. The checkpoint was trained at 9.6 µm.
- A 7-point front seeded on a line of w035 grid vertices, 8 steps, threshold 0.3. villa's own note says this checkpoint is best near 0.3. At its default of 0.5 the C240 trace stopped after 2 steps.
- The neighbour channel is left empty, which is what `inference.py` does when no neighbour is passed. So the tracer sees only the CT and the seed. The published w034 and w036 are not inputs, and comparing the trace with them is not circular.
- CPU, 4 threads. About 10 s per run.
- `torch.load` is forced to `weights_only=True`. The checkpoint's pickle uses only globals on torch's allowlist.

To check that the trace reads the CT, at each place the tracer ran three times: on the real CT, and on the CT shifted by one level-3 voxel outward and inward along the axis nearest the local w035 normal. A tracer that reads the CT should move with the shift. The measure is the mean displacement of the shifted trace from the real one along the outward normal, over steps 5 to 8, divided by the shift's component along the normal. 1.0 means it follows the CT fully.

The places:

- C240 (z 17,995, 240°), the control of the main README, used to pick the settings. It is a tuning place, not a test.
- Four calibration places, picked by a rule written into `pick_places.py` before any of their CT was read. Two are the narrowest eligible places (by the 10th percentile of the spacing from w035 to its nearer neighbour) and two are typical for their height. All four are at least 30° outside the gap sector, 49.9° at the closest, and their CT holds no vertex of the 12 planned controls.

| Place | z | Trace angles | Kind | Spacing to the nearer neighbour, p10 (µm) |
| --- | --- | --- | --- | ---: |
| cal1 | 17,000 | 248.0° to 260.5° | narrow | 45.6 |
| cal2 | 14,856 | 239.1° to 250.6° | narrow | 52.8 |
| cal3 | 12,400 | 173.1° to 180.8° | typical | 70.5 |
| cal4 | 19,500 | 154.4° to 161.8° | typical | 80.1 |

cal1 sits just below C240 and shares 8 CT chunks with it, so it is not independent of the tuning place.

cal2 failed the check, so we read its CT along the w035 normal at every trace point (`cal2_look*.py`). Then, with a rule fixed and hashed before any finer voxel was read (`FINER-RULE.md`), we read levels 1 and 2 at the seed (2.258 µm and 4.516 µm voxels, 14 chunks). Later, with an addendum fixed and hashed before any level-0 voxel was read (`LEVEL0-ADDENDUM.md`), we read level 0 there as well (1.129 µm voxels, 48 chunks).

villa's copy model (`infer_rowcol_triplet_wraps.py`, checkpoint `scrollprize/copy_displacement_latest`) was also run at C240 with the same kind of CT check. It is not part of the pilot, for the reason given under Results.

## Results

### C240

On the real CT the trace reached all 8 steps and kept 301 of 333 points (90%). Every kept point was nearest w035, at a median 10.4 µm (p90 31.6 µm).

With the CT shifted 26 µm outward along the normal, the trace moved over to the shifted sheet within about 5 steps: its mean displacement over steps 5 to 8 was 23.9, 25.5, 26.4 and 28.2 µm. With the same shift inward it traced 2 steps and stopped. With a constant CT it took no step. So the tracer does read the CT.

### The copy model reads the CT only weakly

At C240, with a constant CT, the copy model still put 605 of 625 front vertices nearest w034 and 620 of 625 back vertices nearest w036. Under a 9 µm CT shift its outputs followed by about 0.7 µm, or 8%. Its outputs come mostly from the input mesh and a learned spacing. Inside the gap they would restate w035's position, so we left it out.

### Four calibration places

All 12 runs reached 8 steps.

| Place | Kind | Follow ratio, steps 5 to 8, outward / inward shift | Real CT: kept points to w035, median / p90 (µm) | Real CT: control points nearest w035 |
| --- | --- | --- | --- | --- |
| cal1 | narrow | 0.60 / 0.73 | 12.8 / 32.2 | 53 of 54 |
| cal2 | narrow | 0.34 / -0.37 | 18.5 / 60.3 | 41 of 52 |
| cal3 | typical | 1.00 / 1.03 | 14.8 / 33.8 | 55 of 56 |
| cal4 | typical | 1.03 / 1.06 | 17.0 / 49.1 | 53 of 56 |

The shift was 9.03 µm. Its component along the normal was 8.5 to 8.9 µm at cal1 to cal3 and 6.1 µm at cal4, where the nearest axis is 45.9° from the normal.

At both typical places the trace follows the CT fully. At cal1 it follows most of the shift. At cal2 it follows a third of the outward shift, and the inward shift moves it outward. On the real CT, 61 of cal2's 282 kept points are nearest w034.

### cal2: the trace ends on the ridge that carries w034

The drift happens in trace columns 11 to 21. There, at steps 6 to 8, the w034 mesh sits on the trace's ridge in 10 of 11 columns and w035 in 1. In the other columns it is the other way round: w035 in 19 of 26.

The trace stays on one bright ridge from the seed to step 8. At the seed, in columns 11 to 17, w034 is only 27 to 43 µm from w035, and between -130 and +60 µm along the normal the CT shows that one ridge and no other. The ridge is 56 µm wide at half height. After step 4 the w035 mesh leaves that ridge.

Two readings fit:

- **A.** The seed sits where w034's and w035's sheets touch, and the tracer kept w034's sheet when they parted. This needs one error: a sheet choice at a contact the CT cannot resolve.
- **B.** The seed ridge is w035's sheet alone, and the w034 and w035 meshes both moved off it the same way. This needs two mesh errors.

A needs fewer assumptions. Following ridges in the CT cannot prove it, because at the seed the CT shows one ridge.

Across cal2, the CT separates w034 from w035 where the meshes are at least about 60 µm apart, and mostly does not below that:

| w034 to w035 mesh spacing | separate : merged (lines where both meshes are on a sheet) |
| --- | --- |
| below 60 µm | 5 : 9 |
| 60 to 80 µm | 26 : 4 |
| 80 to 100 µm | 14 : 0 |
| above 100 µm | 20 : 0 |

From the CT alone, adjacent ridges at cal2 are a median 94 µm apart (p10 49.6 µm), and a ridge is a median 64 µm wide at half height.

### Finer levels at the cal2 seed

The rule (`FINER-RULE.md`, sha256 `c7d471962f91e3f6e36e42a18ba1bf3efe3760f874ba55c13f58a23b0b2feea8`) classes each of the 7 seed columns as two sheets (reading A), one sheet (reading B) or other, and lets level 1 decide.

| Level | Voxel (µm) | Two sheets | One sheet | Other | Call |
| --- | ---: | ---: | ---: | ---: | --- |
| 0 (read later) | 1.129 | 0 | 4 | 3 | undecided |
| 1 (decides) | 2.258 | 0 | 4 | 3 | undecided |
| 2 | 4.516 | 0 | 5 | 2 | one sheet |
| 3 | 9.03 | 0 | 5 | 2 | one sheet |

No column shows two ridges at any level. The 3 "other" columns at levels 0 and 1 have their single ridge on w034's side of the midpoint, not two ridges.

The finer levels add no detail along the normal here. The seed ridge's half-height width has a median of 59 µm at level 0, 60 µm at level 1, 59 µm at level 2 and 63 µm at level 3. Its flanks rise from 10% to 90% over a median of 31, 30.5, 28 and 28 µm. Two sheets pressed together with no dip between them would look like this at all four levels, so finding no second ridge leans toward reading B only weakly.

Level 1 is the 2x2x2 mean of level 0 (largest difference 0.5 grey levels), so level 0 is the finest data this scan has. On the normal through each seed column, between the points where w034 and w035 cross it, the deepest dip at level 0 is 0.45 grey levels, under a tenth of the rule's cut, while the crest stands 106 to 136 grey levels above the profile's 10th percentile. Fine texture at level 0 loses its correlation within 2 to 3 voxels, so the scan does show detail of a few µm, and a gap of more than a few µm between two sheets would likely have shown as a dip. That last point is an inference, not a measurement. `LEVEL0-LOOK.md` has the numbers per column.

## What this means for the gap

- In the gap the pitch is about 63 µm, close to the single-ridge width above. The main README found w034 and w036 in one dense band at E1 and E5, and 31 µm apart in one band at E6. Those are the places where a 3D trace would have to choose a sheet, and the CT does not show which one is right.
- The tracer passed the CT check at typical pitch and at one of the two narrow places. It failed at the other, which is closer to the gap's situation. A gap trace could look clean and still be on the wrong winding.
- The planned controls next to the gap edges have the same problem. At every low-edge control, w034 runs within 21 µm of w035 over 23% to 56% of the planned region, so even outside the gap a nearest-mesh test there would be set by mesh placement.

So a trace into the gap, on this scan with this tool, would not settle the winding order. We stopped before the gap run.

## Limits

- Four calibration places and one tuning place. cal1 is not independent of C240.
- One shift size (9 µm) and one shift axis per place. Shifts between 9 and 27 µm were not tried.
- Level 0 (1.129 µm) was read at the cal2 seed only.
- Possible half-voxel offset between CT levels and the level-0 meshes, not corrected. Level-1 crests sit 0 to 8 µm inward of level-3 crests (median 6), and level-0 crests agree with level-1 crests to within 1 µm.
- The w034 and w036 mesh rows were fetched every 3rd row and interpolated in between.
- Both villa checkpoints list PHerc0139 in their training data, so agreement with the published meshes can be partly recall.
- The ridge followers in `cal2_look_ridge.py` and `cal2_look_split.py` are a heuristic, one level-3 voxel per step. They can hop where ridges touch.
- The tracer ran at 9.03 µm voxels against its 9.6 µm training size.
- We did not check whether w034, w035 and w036 were fitted independently. If they share a fitting step, reading B becomes less unlikely.

## How to rerun

Everything ran in a WSL Ubuntu venv with Python 3.14.6 and torch 2.12.1+cpu. The sheet tracer also needs `imagecodecs` to write its LZW tifs. Run Python with `-I`.

1. Check out villa at `a23ec7c` and set `VILLA_SRC` to its `vesuvius/src`. The `scripts/run_*.sh` files run from this folder and take the Python from `PY` (default `python3`); if `imagecodecs` sits in a separate folder, set `EXTRA_SITE` to it. `villa_env.py` puts that checkout first on `sys.path` and forces `weights_only=True`.
2. Fetch the checkpoint with `scripts/fetch_file.py` (the whole `last.pt`, checked against its sha256). All fetches go through `scripts/netlog.py`, which logs every byte and stops at a cap.
3. C240: `fetch_meta.py`, `locate.py`, `make_inputs.py`, `plan_ct.py`, `fetch_ct.py`, `prep_sheet.py`, then `run_sheet7.sh` (which calls `run_sheet.py` and `trace_to_npz.py`) and `analyze_sheet.py`. The CT shift runs: `ablate_build.py`, `ablate_build_small.py`, `run_ablate_sheet.sh`, `ablate_compare.py`.
4. Calibration places: `pick_places.py` (mesh rows only, no CT), `edge_scan.py`, `calib_prep.py` (CT fetch, shifted views, seeds), `run_calib.sh` (which calls `run_sheet_calib.py`), then `calib_compare.py`.
5. cal2: `cal2_look.py`, `cal2_look_png.py`, `cal2_look_ridge.py`, `cal2_look_split.py`. Finer levels: `cal2_look_finer_estimate.py` (HTTP HEAD only), `cal2_look_finer_fetch.py`, `cal2_look_finer.py`, `cal2_look_finer_width.py`. Level 0: `cal2_look_level0_estimate.py` (HTTP HEAD only), `cal2_look_level0_fetch.py`, `cal2_look_level0.py`, `cal2_look_level0_checks.py`.
6. Copy model (optional): `hf_zip_index.py`, `pkl_audit.py`, `hf_zip_member.py` and `build_ckpt.py` fetch only the model weights, without the optimizer state. `run_copy.py`, `analyze_copy.py`, `run_ablate.sh`, `run_ablate2.sh`, `ablate_fig.py`. The copy model imports volume-cartographer's compiled `vc` bindings; `scripts/shim/vc` replaces the calls it makes with numpy over a local zarr. It peaks at 16.8 GiB of RAM at batch 1.

The CT, checkpoints and mesh rows are not in this folder. A full rerun downloads about 1.57 GB, nearly all of it the two checkpoints and the CT chunks, from the public Vesuvius Challenge bucket and Hugging Face.

## Files

- `scripts/`: every script named above, plus `meshgeom.py`, `tif_rows.py`, `check_shim.py`, `sheet_cfg.py`, `rules_numbers.py`, `ray_check.py` and the `shim/` package.
- `out/calib/compare.json`: all calibration results (steps, per-step control points by nearest mesh, distances, shift displacements and ratios). `plans.json` and `views.json` hold the CT plans and the read-back checks of the shifted views. Per run: `trace_*_record.json` and `trace_*_reads.json`.
- `out/places/`: `places.json` (controls, candidates, picks and the rule), `edge_scan.json`, `control_share_below_21.json`.
- `out/cal2_look/`: `CAL2-LOOK.md`, `FINER-RULE.md`, `FINER-LOOK.md`, `LEVEL0-ADDENDUM.md`, `LEVEL0-LOOK.md` and their json files. `png/` holds ribbons along the trace front at steps 2, 5 and 8 and the column 15 sections at levels 3, 1 and 0.
- `out/ablate/`: the C240 shift checks (`compare.json`, `build*.json`, `trace_*_record.json`) and `out/abl_*_metrics.json`.
- `out/sheet7_s8_t03_*`, `out/sheet7_s8_t05_*`, `out/rules_numbers.json`, `out/locate.json`, `out/ray_check.json`: C240 runs and checks.
- `out/png/sheet7_s8_t03_slice.png` (the C240 trace on the CT) and `out/png/abl_copy_panels.png` (the CT each copy-model run saw, with its outputs).
- `logs/finer_rule_sha256.json` and `logs/level0_addendum_sha256.json`: the save time and hash of `FINER-RULE.md` and of `LEVEL0-ADDENDUM.md`. `logs/level0_fetch_record.json` records the level-0 fetch.
