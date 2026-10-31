# LEVEL0-LOOK: the cal2 seed at level 0 (1.129 µm)

Rule: out/cal2_look/FINER-RULE.md, unchanged, sha256 c7d471962f91e3f6e36e42a18ba1bf3efe3760f874ba55c13f58a23b0b2feea8.
Addendum: out/cal2_look/LEVEL0-ADDENDUM.md, saved 2026-11-01 04:37:02 +0800, sha256 64e021755e8815ba489445a84702a8287a860c12ef7d7e0cb43bf4c1757086bd (logs/level0_addendum_sha256.json). The first level-0 request started at 04:39:08 +0800, and its ledger line is stamped 04:39:12. The analysis script checks both hashes before it reads any data.

## Result: one ridge at level 0, as at levels 1 to 3

Level 0 gives 0 A-like columns, 4 B-like and 3 other. Under the rule's call thresholds that is undecided, the same line as level 1. Column by column, the level-0 classes are the same as level 1's.

| level | voxel | cut (prominence) | A-like | B-like | other | call |
|---|---|---|---|---|---|---|
| **0 (this run, its own line)** | 1.129 µm | 5.58 | 0 | 4 | 3 | **undecided** |
| 0 at twice the cut | | 11.15 | 0 | 4 | 3 | undecided |
| 1 (decides under FINER-RULE.md) | 2.258 µm | 5.32 | 0 | 4 | 3 | undecided |
| 1 at twice the cut | | 10.64 | 0 | 4 | 3 | undecided |
| 2 (reported) | 4.516 µm | 6.67 | 0 | 5 | 2 | B |
| 3 (premise check) | 9.032 µm | 17.79 | 0 | 5 | 2 | B |

- **No column shows two ridges at level 0.** Every column has exactly one ridge inside W0, with a prominence of 65 to 120 against a cut of 5.6.
- **The 3 "other" columns have one ridge each, with the crest on w034's half of the window:**
  - Column 11: crest at -30 µm, mid -12.1.
  - Column 16: crest at -16 µm, mid -11.9. Level 2 put this crest at -11 and called the column B-like. Level 0 agrees with level 1 (-15).
  - Column 17: crest at -23 µm, mid -18.2.
- **Rows 1 and 2:** the followers run only for A-like columns. Level 0 has none, so they did not run. All 21 level-0 profiles, and the column paths for all seven columns, lie inside the fetched chunks. No profile is "not checked".

## Width at half height and edge rise

These are descriptive and outside the rule. The method is cal2_look_finer_width.py's.

| level | median half-height width | median 10% to 90% edge | flanks measured |
|---|---|---|---|
| 0 | 59 µm | 31 µm | 10 |
| 1 | 60 µm | 30.5 µm | 10 |
| 2 | 59 µm | 28 µm | 9 |
| 3 | 63 µm | 28 µm | 9 |

Per column at row 0. A dash means the outer flank does not fall to 10% before +60 µm.

| col | L0 crest | L0 width | L0 edge, inner / outer | width L1 / L2 / L3 | L0 level at w034 / w035 |
|---|---|---|---|---|---|
| 11 | -30 | 45 | 23 / - | 46 / 47 / 49 | 0.99 / 0.38 |
| 12 | -14 | 59 | 22 / 57 | 60 / 61 / 63 | 0.63 / 0.86 |
| 13 | 13 | 62 | 27 / 40 | 61 / 62 / 63 | 0.21 / 0.99 |
| 14 | 6 | 60 | 31 / 44 | 60 / 59 / 56 | 0.30 / 1.00 |
| 15 | 4 | 86 | 31 / - | 85 / 84 / 82 | 0.51 / 0.99 |
| 16 | -16 | 54 | 20 / - | 54 / 54 / 76 | 0.45 / 0.84 |
| 17 | -23 | 46 | 31 / - | 45 / 45 / 46 | 0.66 / 0.31 |

Level 0 and level 1 crests agree to within 1 µm in every column. The ridge is no narrower and its edges are no sharper at 1.129 µm than at 9 µm.

## Dips and second ridges

- **Second ridge:** none. No column has a second maximum that passes the cut anywhere in W0, at either cut.
- **Dip between the w034 and w035 positions** (level0_checks.json): the deepest dip strictly between the two crossings is 0.45 grey levels in column 15 (at -9 µm) and 0.37 in column 16 (at -10 µm), which is 0.08 and 0.07 of the cut. In the other five columns the profile has no dip at all between the crossings. For scale, the crest stands 106 to 136 grey levels above the profile's 10th percentile, and σ_prof is 1.31.
- **Two sub-cut bumps lie outside that interval.** Neither one is a dip between w034 and w035:
  - Column 13 has a flat ridge top from about -4 to +14 µm, with a second maximum at -3 µm (prominence 3.2, 0.58 of the cut). The sag that separates it from the crest at +13 sits at +8 µm. That is outward of the w035 crossing (-0.7), not between w034 (-29.1, level 0.21) and w035.
  - Column 17 has a 2.9 bump at +14 µm (0.52 of the cut) in the low ground outward of the ridge. The profile there is about 76 to 79 grey levels, against a crest of 171.
  - Both bumps are also present at level 1 (2.7 each) and smaller at level 2 (0.7 and 1.5). Level 3 does not show them.

## What level 0 adds

This section is descriptive and outside the rule.

- **Level 0 is the finest data there is for this scan.** Level 1 equals the 2x2x2 block mean of level 0. On level-1 chunk (27, 49, 59), the mean difference is 0.25 grey levels and the largest is 0.5, which is rounding. Even-voxel decimation gives a mean difference of 0.85 and a largest of 20. This also confirms the level-0 keys, chunk layout and axis order used here.
- **The level-0 data do carry fine detail:**
  - In level-0 chunk (54, 98, 119), 41% of neighbour differences along x are 0 and 87% are within 1 grey level. Over all 48 chunks σ_vox is 1.05 along x and along z, so the voxel noise is about one grey level.
  - A high-pass residual decorrelates within 2 to 3 voxels along x. The lag correlations are 0.76 at 1.1 µm, 0.46 at 2.3 µm, 0.26 at 3.4 µm and 0.13 at 6.8 µm.
  - The inference, not a resolution measurement: the scan resolves structure of a few µm, so the 30 µm edges and the 59 µm ridge are the density profile at this spot, not voxel blur.
- **What follows for A against B (also an inference):**
  - The meshes put w034 and w035 27 to 43 µm apart at the seed. A lower-density gap of more than a few µm between two sheets there should have shown as a dip at level 0, and none did.
  - So the level-0 evidence leaves two readings. One is a single sheet (B). The other is two sheets pressed together with no gap this CT can see.
  - The rule cannot tell these apart, and its count stays undecided. The reason is that 3 of 7 crests sit on w034's half of the window, not that any column shows two ridges.

## Fetch

- **Estimate** (out/cal2_look/level0_ct_estimate.json, 04:37:39 to 04:38:21): HTTP HEAD only, with the method of cal2_look_finer_estimate.py, whose functions were imported. The estimate found 48 chunks, all HTTP 200 with Content-Length 2,097,152, for 100,663,296 bytes. These are the same 48 keys as the level-0 entry of finer_ct_estimate.json. Every chunk's level-3 parent is in cal2's box, and every chunk passes fetch_ct.outside_gap. The estimate was within the 160,000,000-byte budget.
- **Fetch** (scripts/cal2_look_level0_fetch.py, Windows Python with `py -I`, below-normal priority, 04:39:08 to 04:41:44):
  - All 48 chunks were fetched, none were absent, and the total was 100,663,296 bytes. They went into the new directory downloads/ct_level0/0/z/y/x.
  - The script recomputed the ledger from fetch_log.jsonl before every chunk and checked this run's own tagged bytes (tag `ct-L0`) against its in-process count after every chunk. There was no mismatch, so the WSL stale read did not occur.
- **Ledger:** 1,473,529,433 bytes before and 1,574,192,729 after, with 48 new lines. That leaves 59,336,704 bytes under the 1,633,529,433 limit. An independent re-read of the file gives the same total.
- **The ledger is now 74,192,729 bytes over netlog.CAP_BYTES (1,500,000,000).** This step's budget allowed that. netlog.py was not changed.

## Checks on the code

- The level-0 script imports the rule's functions from cal2_look_finer.py. It recomputes levels 1, 2 and 3 in the same run and gets finer_look.json (noise, cuts, classes, calls, every ridge and prominence) and finer_width.json back exactly. The level-0 numbers therefore come from code that reproduces the earlier results.
- Analysis ran in the WSL venv with `python -I`, `nice -n 10` and 4 threads. It took 12.5 s and peaked at 3.6 GiB resident memory.

## Not checked

- Whether levels 2 and 3 are block means of the level below. Only level 1 against level 0 was checked. If they are, levels 1 and 3 are read 0.5 and 3.5 level-0 voxels off per axis by the rule's level-0 / 2^L sampling, while level 0 has no offset. In six columns the level-3 crest sits 3 to 9 µm outward of the level-0 crest. Column 13 is the exception, at 1 µm inward. The median over all seven is 7 µm outward.
- One sheet's two fibre layers against two touching sheets. This did not come up, because no column shows two ridges.
- The followers for rows 1 and 2. They did not run, because no column is A-like.
- Any place other than cal2.

## Files

| file | contents |
|---|---|
| out/cal2_look/LEVEL0-ADDENDUM.md | The addendum, written and hashed before any level-0 byte was fetched. |
| logs/level0_addendum_sha256.json | The addendum's save time and sha256. |
| out/cal2_look/level0_ct_estimate.json | The HEAD estimate: keys, sizes and the budget check. |
| logs/level0_fetch_record.json | The fetch: start, end, bytes and ledger before and after. |
| out/cal2_look/level0_look.json | Rule numbers for levels 0 to 3, the level-0 followers block, the descriptive width and dip numbers, the reproduction flags, the timeline, the fetch and the checks. |
| out/cal2_look/level0_checks.json | The pyramid check, level-0 texture, the sub-cut bumps and the deepest dip between the crossings. |
| out/cal2_look/png/cal2_col15_L0.png | Column 15, level 0, steps 0 to 2, s -130 to +60 µm, drawn as cal2_col15_L1.png. Cyan is w034, yellow w035 and red + the trace. |
| out/cal2_look/level0_estimate_console.txt, level0_fetch_console.txt, level0_console.txt, level0_checks_console.txt, level0_stderr.txt | Console output. |
| scripts/cal2_look_level0_estimate.py, cal2_look_level0_fetch.py, cal2_look_level0.py, cal2_look_level0_checks.py | The scripts. The fetch script runs under Windows `py -I` and the others in the WSL venv with `python -I`. |
