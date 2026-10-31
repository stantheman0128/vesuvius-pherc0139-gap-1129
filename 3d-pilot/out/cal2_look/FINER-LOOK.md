# FINER-LOOK: the cal2 seed at level 1 (2.258 µm) and level 2 (4.516 µm)

Rule: out/cal2_look/FINER-RULE.md, saved 2026-11-01 02:00:57 +0800, sha256 c7d471962f91e3f6e36e42a18ba1bf3efe3760f874ba55c13f58a23b0b2feea8 (logs/finer_rule_sha256.json). It was written before any level-1 or level-2 voxel was fetched, and it was not changed afterwards.

## Result: undecided

Level 1 decides. It gives 0 A-like columns, 4 B-like and 3 other, and a B call needs 5 B-like.

| level | cut (prominence) | A-like | B-like | other | call |
|---|---|---|---|---|---|
| 1 (decides) | 5.32 | 0 | 4 | 3 | undecided |
| 1 at twice the cut | 10.64 | 0 | 4 | 3 | undecided |
| 2 (reported) | 6.67 | 0 | 5 | 2 | B |
| 3 (premise check) | 17.79 | 0 | 5 | 2 | B |

- **No column shows two sheets at any level.** In every column at the seed there is exactly one ridge inside the window W0 = [s34 - 15, s35 + 15] µm. Prominences run from 62 to 120, against a level-1 cut of 5.3.
- **The 3 "other" columns are not two-sheet columns.** Their single ridge sits on w034's half of the window:
  - Column 11: the crest is at -29 µm, with w034 at -26.8 and w035 at +2.6. The CT level is 0.99 at w034's crossing and 0.41 at w035's.
  - Column 17: the crest is at -24, with w034 at -42.7 and w035 at +6.3. The levels are 0.64 and 0.32.
  - Column 16: the crest is at -15, against a midpoint of -11.9. This is borderline, and level 2 puts the crest at -11 and calls the column B-like.
- **Rows 1 and 2:** the follower part runs only for A-like columns. There are none, so it did not run, and the rule says nothing about which sheet the trace is on at rows 1 and 2.

### Why the finer levels could not settle it

This is descriptive, outside the rule; the numbers are in finer_width.json.

Along the normal at this place, levels 1 and 2 carry no more detail than level 3:
- The single seed ridge's half-height width has a median of 60 µm at level 1, 59 µm at level 2 and 63 µm at level 3.
- Its flanks rise from 10% to 90% over 20 to 31 µm at every level. The medians are 30.5, 28 and 28 µm. Column 14 gives 30, 30 and 30 µm on the inner flank.
- The voxel noise measured at level 1 is lower than at level 3: σ_vox is 1.05 (x) and 2.1 (z) grey levels, against 4.19 at level 3. The data are smooth at fine scales, so the finer voxels add samples but no detail.

Two sheets pressed together with no density dip between them would look like one thick sheet at all three levels. A dip between two sheets 27 to 43 µm apart was the signal that could have shown A. On this evidence it would be no easier to see at level 1 than at level 3. Finding zero two-sheet columns therefore leans toward B only weakly. Under the rule the call is undecided.

## Per column (level 1 and level 2)

s is in µm along the local w035 normal, + outward, 0 at the trace point. The table lists ridges between -60 and +45 µm with their prominence. Ridges over the full -130 to +60 µm range are in finer_look.json.

| col | row | s34 | s35 | level-1 ridges (prominence) | level-2 ridges (prominence) | row-0 class, L1 / L2 / L3 |
|---|---|---|---|---|---|---|
| 11 | 0 | -26.8 | 2.6 | -29 (88.8) | -27 (86.9) | other / other / other |
| 11 | 1 | -33.0 | 11.1 | -15 (93.6), 41 (44.9) | -15 (89), 41 (42.1) | |
| 11 | 2 | -45.6 | 21.5 | 21 (21.6) | 23 (18.4) | |
| 12 | 0 | -29.1 | 0.0 | -14 (97.6) | -7 (96.2) | one sheet / one sheet / one sheet |
| 12 | 1 | -35.5 | 6.9 | 7 (85.9) | 7 (83.6) | |
| 12 | 2 | -42.7 | 16.2 | -20 (38.9), 6 (9.2), 36 (17.2) | -18 (38.7) | |
| 13 | 0 | -29.1 | -0.7 | 13 (114.8) | 13 (113.7) | one sheet / one sheet / one sheet |
| 13 | 1 | -34.0 | 1.9 | -24 (9.9), 9 (75.7) | 5 (73.1) | |
| 13 | 2 | -35.6 | 12.5 | -36 (57.5), 1 (7.6) | -34 (57.9), 0 (7) | |
| 14 | 0 | -28.2 | 4.0 | 5 (120.3) | 6 (119.2) | one sheet / one sheet / one sheet |
| 14 | 1 | -33.4 | 0.8 | 0 (116.2) | 1 (115.5) | |
| 14 | 2 | -33.9 | 12.4 | -8 (92.2) | -7 (92.3) | |
| 15 | 0 | -27.4 | 7.4 | 5 (65), 31 (24.6) | 6 (64), 31 (18.7) | one sheet / one sheet / one sheet |
| 15 | 1 | -32.0 | 2.1 | -26 (7.5), 5 (100) | 6 (97.6) | |
| 15 | 2 | -31.5 | 13.2 | -21 (98), 17 (14.9) | -20 (98.3), 17 (12.4) | |
| 16 | 0 | -31.4 | 7.6 | -15 (111.4) | -11 (108.8) | other / one sheet / one sheet |
| 16 | 1 | -37.2 | 4.9 | -24 (116), 1 (7.1) | -24 (118.3) | |
| 16 | 2 | -33.5 | 15.4 | -41 (111.5), 15 (24.1) | -40 (110), 15 (20.4) | |
| 17 | 0 | -42.7 | 6.3 | -24 (113.6), 33 (20.7) | -22 (112.4), 33 (17.2) | other / other / other |
| 17 | 1 | -46.4 | 8.2 | -36 (107.2), 23 (34.8) | -35 (105.6), 24 (30.8) | |
| 17 | 2 | -38.0 | 18.4 | -45 (93.7), -11 (15.4) | -44 (90.9), -10 (11.6) | |

Noise and cuts:
- Level 1: σ_half 1.25, σ_white 0.52, so σ_prof is 1.25 and the cut is 5.32.
- Level 2: σ_half 1.13, σ_white 1.57, so σ_prof is 1.57 and the cut is 6.67.
- Level 3: σ_prof is cal2_look.py's 4.19, giving a cut of 17.79.

## Fetch

- Exactly the 14 chunks of finer_ct_estimate.json were fetched: 10 at level 1 and 4 at level 2, 29,360,128 bytes. They went through netlog.fetch with the fetch_ct.py URL and tags.
- The ledger moved from 1,444,169,305 to 1,473,529,433 bytes, which leaves 26,470,567 under the 1.5 GB cap.
- The per-chunk cap check read a stale total during the run. Inside the WSL process, netlog.total_fetched() kept returning the pre-run total while lines were being appended to logs/fetch_log.jsonl on /mnt/c. The file itself is correct (14 new lines). The fetch could not pass the cap in any case, because 29.4 MB was below the 55.8 MB that was left. A later fetch run from WSL should not rely on that in-process check.

## Not checked

- **Level 0** (not fetched, because of the cap). Level 1 is no sharper than level 3 here, so level 0 may not be sharper either. That is an inference, not a measurement.
- **Half-voxel offset between levels**, not corrected, as at level 3. Level-1 crests sit 0 to 8 µm inward of level-3 crests (median 6), with level 2 in between.
  - That order and rough size fit sampling at level-0 coordinates divided by 2^L without the centre correction. For a seed normal near (0.22, -0.92, -0.33), the correction predicts level 3 reading about 3.5 to 4 µm outward of level 1. This was not computed per column.
  - If the meshes are level-0 exact, the level-3 sheet positions in CAL2-LOOK.md may sit about 4 µm too far outward. That is small next to the 15 µm windows used there.
- **One sheet's two fibre layers versus two touching sheets.** The rule cannot tell these apart (stated in advance). It did not come up, because no column showed two ridges.
- **The followers for rows 1 and 2.** They did not run, because there was no A-like column.

## Files

| file | contents |
|---|---|
| out/cal2_look/FINER-RULE.md | The rule, fixed before the data were read. |
| logs/finer_rule_sha256.json | The rule's save time and sha256. |
| out/cal2_look/finer_look.json | The rule's numbers: noise, cuts, per level, column and row s34, s35, ridges and prominences, classes, calls. |
| out/cal2_look/finer_width.json | Descriptive, outside the rule: ridge half-height width, 10% to 90% edge, CT level at the w034 and w035 crossings. |
| out/cal2_look/png/cal2_col15_L1.png | Column 15, level 1, steps 0 to 2, s -130 to +60 µm. Cyan is w034, yellow w035, red + the trace. No follower paths are drawn, because the column is not A-like. |
| out/cal2_look/finer_console.txt, finer_width_console.txt, finer_fetch_console.txt | Console output. |
| scripts/cal2_look_finer.py, cal2_look_finer_width.py, cal2_look_finer_fetch.py | The scripts. Run them in the WSL venv with `python -I`. |
