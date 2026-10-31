# PHerc0139: the w034/w036 gap in the 1.129 µm scan

Code and results for one question from [ScrollPrize/villa#2167](https://github.com/ScrollPrize/villa/issues/2167).

## The question

In the published PHerc0139 segments, w035 covers only part of its turn. Between z 9,000 and 27,800 there is a sector, a median 98° wide, where w034 and w036 are traced and w035 is not. In that gap w036 lies close to w034's radius instead of two windings out, and on some heights it lies inside w034. So what is the order of w034 and w036 there, and is there a sheet between them where w035 should be?

The note in #2167 looked at the 2.399 µm CT. Nobody in the thread had used the 1.129 µm scan. This repository reads that scan at the places the note names that fall inside it.

Frame: positions are full-resolution voxels (x, y, z) of volume `20260102150214` (2.399 µm). Angles are degrees about the umbilicus `20260102150214-umbilicus-20260826112529.json` at the same z, atan2(y − y_u, x − x_u), 0 to 360. Radii r are µm from the umbilicus along the ray.

## The scan and how much of the gap it covers

The volume is `20260413113053-1.129um-0.2m-59keV-masked.zarr`, shape (z, y, x) 19,297 × 35,971 × 32,318, uint8, 128³ chunks stored uncompressed (2,097,152 bytes each). The bucket catalog maps it to `20260102150214` with the same matrix as the zarr's `transform.json`.

Through that matrix the scan covers z 11,394 to 20,498 of the #2167 frame. That is 9,104 of the gap's 18,800 z voxels (48%), or 46 of its 95 planes. Inside that z band it covers every angle of the sector: all 13,809 w034 vertices and all 11,406 w036 vertices in the angle envelope (327° to 85°) with z 11,420 to 20,460 fall inside the scan and inside non-empty chunks. Over the whole sector (z 9,000 to 27,800) that is 42.9% of w034's vertices and 42.1% of w036's.

Not covered: z 9,000 to 11,394 and z 20,498 to 27,800. Those heights include the bands where w036 sits furthest inside w034 (median −96 µm at z 23,000 to 24,800).

Places inside the scan: E1, E2, E5, E6, G2, the note's gap examples w034 434,248 (here X434) and w036 542,1016 (X542), and both polar unwraps (z 14,856 and 19,500). Outside: E4, G1, G3, G4, C2 and the #2167 body example at z 24,543.

A second registration of this scan, the one in [villa#1381](https://github.com/ScrollPrize/villa/issues/1381), agrees with the catalog chain to a median of 20 µm and at most 42 µm over the whole volume.

Script: [`coverage/cov2167.py`](coverage/cov2167.py) with its output in [`results/coverage/cov2167.txt`](results/coverage/cov2167.txt). `cov2167b.py` checks how far the gap lies from the edge of the scan, `cov2167c.py` sizes a full download, and `cov2167d.py` compares the two registrations.

## Method

At each place:

1. The ray starts at the umbilicus at the place's z and runs at the place's angle. Mesh crossings come from the published 2.399 µm tifxyz copies (48 µm grid) by ray-triangle intersection, every 10 µm in dz.
2. Radial cut: the 1.129 µm CT is sampled along the ray at 1 µm steps in r and dz ±200 µm, trilinear, each sample the mean over ±3 µm across the line.
3. Z slices at dz −150, 0 and +150 µm, each the mean of 5 planes (±2.3 µm).
4. The cut is smoothed (Gaussian, sigma 2 µm) and thresholded with Otsu over the whole cut. An opening follows, and bodies under 200 px are removed.
5. Along each of the 41 rows dz −200 to +200 µm, material runs are taken (runs under 6 µm dropped, air gaps under 3 µm closed). Each mesh crossing is assigned to the run it is on or the nearest one. The count is the number of runs strictly between the w034 run and the w036 run.
6. The count is repeated with every mesh moved −30 µm and +30 µm along r.
7. As a second reading on the ray, peaks are taken on the profile: the mean of dz −5 to +5 µm, sigma 2 µm, peaks with prominence at least 15% of the range and at least 10 µm apart.

The polar unwraps use the same thresholding and run counting on one z slice, per 1° angle row from 55° to 95° (41 rows) instead of per dz row.

In the `*_layers.json` files, `d46_0` is the run index of w036 minus that of w034 at one height: 0 means the same run, 1 adjacent runs (no sheet between), 2 one sheet between. `d46_hist_-30` and `d46_hist_30` are the same with the shift. A mesh's "end" is the outer of its two crossings on the ray where it crosses twice (a doubly traced end), and `enddiff_<mesh>` is the run index of the end minus that of the inner piece.

## Results

At E1, E5, E6, X434, X542 and G2 the count finds no sheet between w034 and w036 at any of the 41 heights. E2 is the one place with a sheet between them. The ±30 µm shift changes the count at only 2 to 5 heights, all at X542. At E1, E5 and E6 that zero means the bands are pressed together, not that a sheet is absent. An independent re-count of X542 and E2 (see "Independent check") confirms E2 and qualifies X542.

| Place | z | Angle | Sheets between w034 and w036 | With ±30 µm shift | Can the bands be separated? |
| --- | ---: | ---: | --- | --- | --- |
| E2 (pilot) | 17,995 | 359.805° | 1 at 36 of 41 heights, 0 at 5 | count unchanged | yes |
| E1 | 14,856 | 35.16° | 0 at 41 | 0 at 41 | no |
| E5 | 13,741 | 4.16° | 0 at 41 | 0 at 41 | no |
| E6 | 14,109 | 18.79° | 0 at 41 | 0 at 41 | one band, air on both sides |
| X434 | 12,837 | 12.66° | 0 at 41 | 0 at 41 | meshes 5 µm apart and crossing |
| X542 | 18,432 | 346.17° | 0 at 41 | −30: 1 at 2; +30: 1 at 5 | yes |
| G2 | 19,500 | 15.0° | 0 at 41 | 0 at 41 | yes, ray oblique to the sheets |

### E2, z 17,995 (pilot)

Five material bands, with air gaps 30 to 170 µm wide between them, cross the ray between r 9,380 and 10,140. Crossings at dz 0: w033 9519.1, w034 9629.3, w033 end 9784.3, w034 end 9971.0, w036 9997.8, w037 10092.7.

- w034's inner piece to w036: exactly one sheet between at 36 of 41 heights. That sheet (9768 to 9796 at dz 0) carries only w033's doubly traced end. At dz +100, +110, +130, +140 and +150 the count is 0, where w034's band touches the band carrying w033's end.
- w034's doubly traced end is in air 8 µm inside w036's band (9979 to 10021) at dz 0, and on the same band as w036 at all 41 heights.
- Order along the ray at every height: w033 < w034 < w033 end < w034 end < w036 < w037.
- Otsu threshold 68. Between w034 and w036 the profile drops to 12 to 30 on both sides of the middle band, which rises to 103.

Figures: [radial cut](png/pilot/E2_radialcut.png), [profile](png/pilot/E2_profile.png), [slice](png/pilot/E2_slice.png), [slices at three heights](png/pilot/E2_slices_3heights.png) and [zoomed](png/pilot/E2_slices_3heights_zoom.png), [single-line cut](png/pilot/look_E2_cut1.png). Table: [`results/pilot/E2_layers.json`](results/pilot/E2_layers.json).

### E1, z 14,856, 35.16°

Crossings at dz 0: w033 6219, w034 6306, w036 6327, w037 6402. On the ray w034, w036 and w037 lie in one dense run, 6301 to 6532 at dz 0. Between w034 and w036 the profile stays at 113 to 130 against an Otsu threshold of 89, and no air gap appears at any height. The 1.129 µm scan does not resolve this band better than the 2.399 µm scan did, so the 0 here means touching and unresolved.

Figures: [radial cut](png/places/E1_radialcut.png), [zoom](png/places/E1_radialcut_zoom.png), [profile](png/places/E1_profile.png), [slice](png/places/E1_slice.png), [three heights](png/places/E1_slices_3heights.png), [slices zoom](png/places/E1_slices_zoom.png). Table: [`results/places/E1_layers.json`](results/places/E1_layers.json).

### E5, z 13,741, 4.16°

Crossings at dz 0: w033 6656, w034 6814, w036 6844, w037 6941. w034 sits on the inner face of a dense run (6812 to 6959 at dz 0) just outside a 40 µm air gap, and w036 lies 30 µm further out in the same run. With the −30 µm shift the two are on adjacent runs at 26 heights and on the same run at 15. Nothing separates w034 from w036 at any height.

Figures: [radial cut](png/places/E5_radialcut.png), [zoom](png/places/E5_radialcut_zoom.png), [profile](png/places/E5_profile.png), [slice](png/places/E5_slice.png), [three heights](png/places/E5_slices_3heights.png), [slices zoom](png/places/E5_slices_zoom.png). Table: [`results/places/E5_layers.json`](results/places/E5_layers.json).

### E6, z 14,109, 18.79°

Crossings at dz 0: w033 6126, w034 6233, w036 6264, w037 6378. w034 and w036 lie 31 µm apart in one band 60 to 80 µm wide (6220 to 6294 at dz 0), bounded by air on both sides at all 41 heights. No seam is visible between the two meshes. Whether that band is one sheet or two pressed plies cannot be read from the cut.

Figures: [radial cut](png/places/E6_radialcut.png), [zoom](png/places/E6_radialcut_zoom.png), [profile](png/places/E6_profile.png), [slice](png/places/E6_slice.png), [three heights](png/places/E6_slices_3heights.png), [slices zoom](png/places/E6_slices_zoom.png). Table: [`results/places/E6_layers.json`](results/places/E6_layers.json).

### X434, the note's w034 example (row 434, column 248), z 12,837, 12.66°

Crossings at dz 0: w033 6183, w034 6281, w036 6287, w037 6425. w034 and w036 are in the same band at every height and 5 µm apart at dz 0. They cross each other at about dz +25 µm: w036 is outside w034 at 24 heights and inside at 17. The sheets are strongly inclined in the cut, so the count says little here.

Figures: [radial cut](png/places/X434_radialcut.png), [zoom](png/places/X434_radialcut_zoom.png), [profile](png/places/X434_profile.png), [slice](png/places/X434_slice.png), [three heights](png/places/X434_slices_3heights.png), [slices zoom](png/places/X434_slices_zoom.png). Table: [`results/places/X434_layers.json`](results/places/X434_layers.json).

### X542, the note's w036 example (row 542, column 1016), z 18,432, 346.17°

Crossings at dz 0: w033 9711, w034 9878, w033 end 9896, w034 end 9997, w036 10004, w037 10111.

- An air gap about 100 µm wide lies between w034's band (9811 to 9890 at dz 0) and w036's band (9993 to 10035). The profile drops to 22 there, against an Otsu threshold of 85.
- w034 and w036 are on adjacent runs at 37 heights. At the other 4 (dz −200 to −170) w036 runs 76 to 84 µm into air below the lower end of its band and the nearest-run rule puts it on w034's band.
- So on this ray w036 sits on the next band out from w034, with no sheet between. w036 is outside w034 at 41 of 41 heights.
- The independent check below found two things this count misses. A sheet enters the gap from just above the window: at dz +200 it is starting to cross between the two meshes, and half a degree lower in angle it lies between them from dz +110 up. Lower down, at dz −160 to −60, w034's band splits into two plies. So X542 does not show that a winding is missing here. The missing sheet could be the one entering the gap or one of the two plies, and only a 3D trace can tell.
- w034's band is 80 to 110 µm wide and fibrous in the slices, so it could be two pressed plies.
- w034's doubly traced end is on w036's band, one band out from w034's inner piece at 36 of 41 heights (same band at 5). w033's end runs along the outer face of w034's band; it is one band out from w033's inner piece at 30 of 41, two at 10, three at 1.

Figures: [radial cut](png/places/X542_radialcut.png), [zoom](png/places/X542_radialcut_zoom.png), [profile](png/places/X542_profile.png), [slice](png/places/X542_slice.png), [three heights](png/places/X542_slices_3heights.png), [slices zoom](png/places/X542_slices_zoom.png), [slices zoom, r 9,640 to 10,060](png/places/X542_slices_zoom_mid.png). Table: [`results/places/X542_layers.json`](results/places/X542_layers.json).

### G2, z 19,500, 15.0°

Crossings at dz 0: w033 9684, w034 9880, w036 9922, w037 10079. w034 and w036 share one band (9816 to 9969 at dz 0) above dz 0. From dz −10 down, a thin crack splits that band; w034 lies in the crack, 7 to 11 µm from material, and w036 is on the outer ply. The count is 0 either way (same run at 22 heights, adjacent at 19). Below dz 0 two plies between w033 and w036 carry no trace, which matches the note's 2.399 µm reading. The sheets run obliquely to the ray and branch near it.

Figures: [radial cut](png/places/G2_radialcut.png), [zoom](png/places/G2_radialcut_zoom.png), [profile](png/places/G2_profile.png), [slice](png/places/G2_slice.png), [three heights](png/places/G2_slices_3heights.png), [slices zoom](png/places/G2_slices_zoom.png). Table: [`results/places/G2_layers.json`](results/places/G2_layers.json).

### Polar unwraps at z 14,856 and 19,500, 55° to 95°

- z 14,856, gap part 55° to 75° (21 rows): w034 and w036 share a band at 17 rows, and w036 lies one band inside w034 at 4 rows (58°, 60°, 61° and 71°). There is never a sheet between them. From 81° w035 is on a band between them; by eye that band runs back to at least 76°, carrying w037's piece that ends at 78°.
- z 19,500, gap part 55° to 78° (24 rows): no sheet between w034 and w036 at 24 of 24 rows (same band at 18, adjacent at 6). Past 79° the bands are pressed even where w035 is present (84° to 95°, where w035 shares a run with w034 or w036 at 6 of 12 rows), so counts there are weak.
- Where w036 or w037 crosses a ray twice near the gap end, the two crossings are mostly two bands apart (z 14,856: w036 at 4 of 4 rows, w037 at 8 of 9; z 19,500: w036 at 2 of 3, w037 at 7 of 10), with the other mesh's piece on the band between.
- At z 19,500, 55° to 64°, w033 runs through an air pocket, up to about 550 µm from material.

Figures: [unwrap z 14,856](png/places/U14856_unwrap.png), [unwrap z 19,500](png/places/U19500_unwrap.png). Tables and logs: [`results/places/U14856_layers.json`](results/places/U14856_layers.json), [`U19500_layers.json`](results/places/U19500_layers.json), [`U14856_log.txt`](results/places/U14856_log.txt), [`U19500_log.txt`](results/places/U19500_log.txt).

## Control at 240°, z 17,995

The control lies outside the w035 gap at the height of E2, with all five meshes in the published order: w033 8130.8, w034 8257.8, w035 8382.6, w036 8550.9, w037 8679.3. Among the candidate angles, 240° had the most even spacing (127 to 168 µm).

- Runs at dz 0: 8125 to 8158, 8237 to 8292, 8319 to 8360, 8530 to 8688, 8749 to 8830. Otsu threshold 86.
- w034 to w036: one sheet between at 36 of 41 heights. That sheet is the one w035 runs beside; w035 lies 7 to 25 µm outside it, in air, from dz −80 to +80.
- w036 and w037 share one run at dz −200 to +100, with two peaks 78 µm apart, which looks like two sheets pressed together.

So where w035 is present, the 1.129 µm scan shows one sheet between w034 and w036, with each mesh on its sheet or up to about 25 µm beside it. Some sheets are pressed into one body even there. E2 has the same count. The other gap places do not.

Figures: [radial cut](png/pilot/C240_radialcut.png), [profile](png/pilot/C240_profile.png), [slice](png/pilot/C240_slice.png), [single-line cut](png/pilot/look_C240_cut1.png). Table: [`results/pilot/C240_layers.json`](results/pilot/C240_layers.json).

## Independent check

A separate implementation re-counted X542 and E2. It shares no code with this repository: its own ray and triangle crossings, its own byte-range CT reader, the zarr's `transform.json` inverted (landmark residuals 0.3 to 1.95 voxels), and a different threshold rule, T = a + f(p − a) with a and p the 10th and 90th percentiles of the smoothed cut, at f = 0.5, 0.35, 0.2 and 0.1. It also starts each height's ray at the umbilicus at that height, so its crossings differ from the ones above by up to 11 µm at X542 and 47 µm at E2 away from dz 0. It downloaded 45,137,206 bytes. Its code is not included here.

- E2: one sheet between w034 and w036 at 41 of 41 heights at f = 0.5, 36 of 41 at f = 0.35 (zeros at dz +100 to +140, where w034's band touches the middle band) and 30 of 41 at f = 0.2. Moving the meshes 42 µm out or the ray 0.5° raises the count at some heights. It never drops to zero everywhere, so the method can see the sheet that is there.
- X542 on the nominal ray: no sheet between at 40 of 41 heights at f = 0.5 and at 41 of 41 at f = 0.35 and 0.2. Inside the gap the largest peak is 4.8 grey levels at f = 0.1, against air at about 28 and papyrus at about 143. Around dz 0 the gap is real air 85 to 93 µm wide.
- X542 at the top: the gap narrows to 24 to 40 µm at dz +160 to +190 and 11 µm at dz +200. A small peak at r 9944 to 9956 rises with height (86 at dz +200, 92 at +204, 99 at +208), so a sheet is entering from just above the window. At −0.5° it lies between w034 and w036 at 9 heights, dz +110 to +200. At +0.5° the count stays 0 at 41.
- X542 lower down: w034's band splits into two plies at dz −160 to −60, and w034 sits on the outer one. Moving the meshes 30 or 42 µm inward puts w034 on the inner ply, and the outer ply then counts as one sheet between, at 2 and 5 heights.
- Its angle for X542 is 346.1619° against 346.1683° here, about 1 µm sideways at this radius. Counts agree with `results/places/X542_layers.json` except at three heights: dz −180 and −170 (same band here, adjacent there) and dz +200 (0 here, 1 there).

Figures: [X542](png/verify/X542_verify.png), [X542 at ±0.5°](png/verify/X542_verify_angles.png), [X542 gap zoom](png/verify/X542_verify_slices_zoom.png), [E2](png/verify/E2_verify.png), [E2 at ±0.5°](png/verify/E2_verify_angles.png).

## What this cannot settle

- A band can be two sheets pressed together, and the method cannot tell. That applies at E1, E5 and E6, and to w034's band at X542. Where it applies, "0 sheets between" could be "1 sheet, unresolved".
- The count does not say which winding a band is. That would need the bands followed in 3D over a wider area than a radial cut and a few slices per place.
- The catalog registration and the #1381 registration differ by 20 to 42 µm, against a layer pitch of about 63 µm in the gap. The ±30 µm check covers only a common shift along r, not a local warp or a mosaic tile seam.
- At G2 and X434 the ray crosses the sheets at a slant, so counts along it mix in sheets that bend across it.
- In the unwraps, which crossing counts as the "end" follows the outer-crossing convention, which does not always match the piece that stops.

## 3D follow-up

Wittmoon asked in #2167 whether the entering sheet and the two plies at X542 can be followed in 3D ([comment](https://github.com/ScrollPrize/villa/issues/2167#issuecomment-6537933181)). [`3d-pilot/`](3d-pilot/README.md) is a pilot with villa's sheet tracer, run at the 240° control and at four places outside the gap to test whether its trace follows the CT. It stopped before the gap. At one calibration place the trace drifted from w035 onto the ridge that carries w034; the two meshes are 27 to 43 µm apart at its seed, and the CT there shows a single ridge at levels 3, 2 and 1. In the gap the pitch is about 63 µm and the bands at E1, E5 and E6 are already pressed together, so a trace there meets the same question with nothing in the CT to settle it.

## Data sources

All CT and mesh files are public, from the bucket `s3://vesuvius-challenge-open-data` (HTTPS: `https://vesuvius-challenge-open-data.s3.amazonaws.com/<key>`). No credentials are needed.

| What | Key |
| --- | --- |
| Bucket catalog (stored gzip-encoded) | `metadata.json` |
| Umbilicus | `PHerc0139/representations/umbilicus/20260102150214-umbilicus-20260826112529.json` |
| 1.129 µm CT, level 0 chunks | `PHerc0139/volumes/20260413113053-1.129um-0.2m-59keV-masked.zarr/0/<z>/<y>/<x>` |
| Its metadata | same zarr: `0/.zarray`, `metadata.json`, `transform.json` |
| Level 5 chunks (coverage sanity check) | same zarr: `5/1/5/5`, `5/3/4/6` |
| Level 3 key listing (coverage) | ListObjectsV2 of the same zarr's `3/` prefix |
| Meshes, 9.362 µm copies | `PHerc0139/segments/<id>/mesh/<ts>-on-20250728140407-9.362um.tifxyz/{x,y,z}.tif` and `meta.json` |
| Meshes, 2.399 µm copies | `PHerc0139/segments/<id>/mesh/<ts>-on-20260102150214-2.399um.tifxyz/{x,y,z}.tif` and `meta.json` |

Segment ids (`<ts>` is the part before the first `-`):

| Mesh | `<id>` |
| --- | --- |
| w033 | `20260219000000-w033_2026021913` |
| w034 | `20260303000000-w034_2026030317` |
| w035 | `20260317000000-w035_2026031718` |
| w036 | `20260316000000-w036_2026031607` |
| w037 | `20260310000000-w037_2026031015` |

The #1381 transform used by `coverage/cov2167d.py` is `20260413113053-to-20250728140407.json` at the root of [7jycwjmbfn-eng/pherc0139-physical-audit](https://github.com/7jycwjmbfn-eng/pherc0139-physical-audit) at commit `3f2bc75e63fbdabb3926e4d550b3a5851a231dfb`.

## Reproducing

Tested with Python 3.13.3, numpy 2.3.1, scipy 1.16.2, matplotlib 3.10.6 and tifffile 2026.4.11 (`requirements.txt`), on CPU. The scripts were run with `python -I`; they put the repository root and the user site directory on the path themselves.

Two environment variables set the folders:

- `PHERC0139_DATA`: download cache. Default `data/` at the repository root.
- `PHERC0139_OUT`: outputs (`.npz` caches, PNGs, `*_layers.json`). Default `out/` at the repository root.

The `.npz` files (134,953,336 bytes for all places) are caches of the resampled CT and the mesh crossings. They are not in git. Running the sampling step for a place (`p3_sample.py`, `f2_sample.py`, `f5_unwrap.py`) writes them again.

One place, X542:

```bash
python -I fetch_inputs.py              # catalog, umbilicus, 9.362 µm meshes
python -I places/f1_geom.py            # crossings at all six places -> out/places/geo.json
python -I places/f2_sample.py X542     # CT planes and mesh rows for X542 -> out/places/X542.npz
python -I places/f3_analyze.py X542    # -> out/places/X542_layers.json and the figures
python -I places/f4_zoom.py X542
python -I places/f4_zoom.py X542 9640 10060 _mid
```

The other places: `f2_sample.py` and `f4_zoom.py` with `E1`, `E5`, `E6`, `X434` or `G2`, then `f3_analyze.py E1 E5 E6 X434 X542 G2`. The unwraps: `f5_unwrap.py U14856 14856` and `f5_unwrap.py U19500 19500`.

The pilot (E2 and the 240° control):

```bash
python -I pilot/p1_geom.py
python -I pilot/p2_fetch.py            # -> out/pilot/plan.json and the CT chunks
python -I pilot/p3_sample.py
python -I pilot/p4_look.py E2 C240
python -I pilot/p5_analyze.py E2 C240
python -I pilot/p6_levels.py           # two more z slices at E2
python -I pilot/p7_levels_fig.py
python -I pilot/p7b_zoom.py
```

The coverage check: `python -I fetch_inputs.py --coverage`, then `python -I coverage/cov2167.py` (and `cov2167b.py`, `cov2167c.py`, `cov2167d.py`).

A run of X542 alone into an empty data folder gave a `X542_layers.json` byte-identical to the committed one and the same PNGs.

`p4_look.py` also writes `look_*_cut.png` and `look_*_slice.png`. They show the same as the top panel of `*_radialcut.png` and as `*_slice.png`, so they are left out of `png/`.

### Bytes downloaded

The CT is fetched as byte ranges of the uncompressed chunks, only the z planes each place needs. The 2.399 µm mesh tifs are uncompressed too, so only a band of grid rows around each place is fetched. The band is picked from the 9.362 µm copy and checked against the full-tif crop.

| Place | CT |
| --- | ---: |
| E1 | 47.7 MB |
| E5 | 55.1 MB |
| E6 | 48.2 MB |
| X434 | 47.6 MB |
| X542 | 66.7 MB |
| G2 | 55.7 MB |
| U14856 | 31.9 MB |
| U19500 | 47.9 MB |
| E2 and the 240° control (pilot) | 233.1 MB |

On top of the CT, each place fetches 4.7 to 6.0 MB of mesh rows, and `fetch_inputs.py` fetches 8,205,629 bytes once. A run of X542 alone into an empty `PHERC0139_DATA` downloaded 80,907,163 bytes: 8,205,629 from `fetch_inputs.py`, 66,732,032 of CT planes, 5,967,552 of mesh rows and 1,950 of tif headers.

## Layout

```
coverage/   coverage of the gap sector by the scan (cov2167*.py)
3d-pilot/   the 3D pilot with villa's sheet tracer (its own README)
pilot/      E2 and the 240° control (pcommon.py, p1 to p7b)
places/     the other places and the unwraps (fcommon.py, f1 to f5)
results/    per-height tables (*_layers.json), logs, geo.json, plan.json, coverage output
png/        figures
fetch_inputs.py, repo_paths.py
```

`results/places/analyze_log.txt` is the printed output of `f3_analyze.py` for E5, E6, X434, X542 and G2.

## Credit

The places and the question come from ChaseHendrick's note in #2167 ([the note at commit 1243518](https://github.com/ChaseHendrick/Scrolls/blob/1243518dae102dd334ff1d7a63995f323c1095e4/docs/logs/2026-10-14-pherc0139-segment-defects.md), [his comment](https://github.com/ScrollPrize/villa/issues/2167#issuecomment-6326350177)) and Wittmoon's audit in the same thread ([comment](https://github.com/ScrollPrize/villa/issues/2167#issuecomment-6313289842)).

## License

MIT, see [LICENSE](LICENSE).
