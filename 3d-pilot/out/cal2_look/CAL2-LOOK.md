# cal2 look: what happens where the trace drifts toward w034

Place cal2_z14856_a245, real-CT trace, steps 0 to 8. Level-3 CT (9.03 µm voxels), trilinear, profiles along the local w035 outward normal (+ is toward w036) over ±200 µm. CT read: cal2's own 30 chunks, already on disk. Nothing was fetched, and no control or gap-envelope CT was read. 2026-11-01.

## 1. Did the trace leave the sheet, or is the w035 mesh off it?

The CT cannot settle this. The trace does not leave papyrus. It ends on the sheet that the w034 mesh follows, and the CT cannot tell whether that sheet is w034's or w035's, because at the seed the two show as one ridge.

The numbers below come from columns 11 to 21, where the drift happens, compared with the other 26 columns.

- The trace stays on one bright ridge from the seed to step 8. All 37 of its on-sheet points at steps 5 to 8 lead back to the seed's ridge. It runs at gap level for 1.9% of its path (median per column), against 13.8% in the other columns.
- That ridge also carries the w034 mesh. At steps 6 to 8 the w034 mesh sits on the trace's ridge in 10 of 11 drift columns and w035 in 1. In the other columns it is the other way round: w035 in 19 of 26, w034 in 7. The w034 mesh's on-sheet points at steps 1 to 4 lead back to the seed's ridge in 20 of 25 cases, against 20 of 54 elsewhere.
- At the seed, w034 is only 27 to 43 µm from w035 in columns 11 to 17. Between -130 and +60 µm the CT shows one ridge and no other. w034 sits on that ridge in columns 11 and 17. In columns 12 to 16 it sits on the ridge's lower flank, below half height (levels 0.13 to 0.45). The ridge is 56 µm wide at half height, the same as elsewhere (53.5). If the w034 mesh is within about 100 µm of its sheet, w034's sheet is inside the seed's ridge.
- After step 4, w035 leaves that ridge. Of its on-sheet points at steps 5 to 8, 8 of 21 lead back to the seed's ridge.
- The w035 mesh is not more off-sheet in the drift columns than elsewhere: it is at gap level for 21.2% of its path there and 20.5% in the other columns. Gap time alone does not mark w035 as the faulty surface.

Two readings fit these numbers:

- **A.** The seed sits where the w034 and w035 sheets touch, and the tracer kept to w034's sheet when they parted. That needs one error, a sheet choice at a contact that the level-3 CT cannot resolve.
- **B.** The seed ridge is w035's sheet alone, and the w034 and w035 meshes both moved off it in the same direction. That needs two mesh errors, and the w034 mesh would be more than 100 µm from its own sheet at the seed.

A needs fewer assumptions. It is not proven, because following ridges in the CT cannot tell A from B when the two sheets form one ridge at the seed.

The earlier follower verdict in ridge_follow.json (5 columns where only the trace's ridge reaches the seed's ridge, 0 where only w035's does) does not support B. The seed's ridge is also w034's.

Per step (all 37 columns; steps 5 and 6 keep 30 points):

| step | trace on a sheet | w035 on a sheet | trace on w034's ridge | w034/w035 separate : merged : not judged | median s of w034 / w036 (µm) |
|---|---|---|---|---|---|
| 0 | 25 | 25 | 0 | 5 : 2 : 30 | -66 / 99 |
| 1 | 32 | 31 | 0 | 6 : 4 : 27 | -69 / 102 |
| 2 | 26 | 26 | 1 | 7 : 4 : 26 | -64 / 111 |
| 3 | 21 | 22 | 0 | 13 : 3 : 21 | -62 / 123 |
| 4 | 18 | 20 | 3 | 10 : 0 : 27 | -58 / 132 |
| 5 | 19 of 30 | 14 of 30 | 6 | 4 : 0 : 26 | -42 / 147 |
| 6 | 21 of 30 | 16 of 30 | 6 | 5 : 2 : 23 | -36 / 160 |
| 7 | 30 | 20 | 10 | 7 : 0 : 30 | -44 / 130 |
| 8 | 30 | 23 | 11 | 13 : 0 : 24 | -48 / 108 |

At steps 3 and 4 the trace is off any sheet at about half the points (16 and 19). In columns 12 and 15 it crosses a dark lens there, where the seed's band splits into two strands (see png/cal2_col15.png).

## 2. Can the CT tell w034, w035 and w036 apart at about 63 µm?

From about 60 µm upward, yes, in most lines where both meshes sit on a sheet. Below 60 µm, mostly not.

- w034/w035 at steps 1 to 8, separate : merged by the meshes' own spacing:

  | spacing | separate : merged |
  |---|---|
  | below 60 µm | 5 : 9 |
  | 60 to 80 µm | 26 : 4 |
  | 80 to 100 µm | 14 : 0 |
  | above 100 µm | 20 : 0 |

- w035/w036 at 60 to 80 µm: 14 : 1.
- From the CT alone, adjacent ridges are p10 49.6, p25 67 and median 94 µm apart. 18.7% of these spacings are below 60 µm and 4.3% below 40 µm. Ridge width at half height has a median of 64 µm.
- Most lines could not be judged because one mesh crossing sits in a gap. For w034/w035 at 60 to 80 µm, 30 lines were judged and 71 were not.

The seed in the drift columns, where w034 and w035 are 27 to 43 µm apart, falls in the range where the CT shows one ridge.

Cut-offs: the main cut is ridge prominence ≥ 17.79 counts, which is 3·√2·σ. σ is the larger of the in-plane voxel noise values (MAD-based: x 3.15, z 4.19). noise2 (11.86) and the tangent-line null (p90 32.34, p95 48.02) are reported in profiles.json too. The null runs high because tangent lines cross real thin sheets, so it is only a conservative bound. "On a sheet" means at or above half of the nearest ridge's height over the line's p10.

## Not checked

- A possible half-voxel offset (up to about 4 µm) between level-3 CT and the level-0 meshes.
- The ±1 shift runs. Only the real-CT trace was examined.
- w034 and w036 mesh rows exist only every 3rd row, and are interpolated in between.
- Whether w034, w035 and w036 are independent fits. If they share a fitting step, reading B's two-mesh error becomes less unlikely.
- Finer CT (levels 0 to 2) at the seed of columns 11 to 17. That would separate A from B, and it is the one check worth doing next.
- The ridge followers are a heuristic: one level-3 voxel per step, slope at most 1. They can hop where ridges touch. The ridge_follow.json counts hold for slope 0.5 and for cut noise2. split_check.json was run at noise3 only, with slopes 1.0 and 0.5; the shared-ridge counts (10 vs 1, and 7 vs 19) are the same at both slopes.

## Files (all under the place folder)

| file | contents |
|---|---|
| out/cal2_look/profiles.json | Per step and point: ridges (position, prominence), w034/w035/w036 crossings, trace position, distances and levels, pair calls. Also per-step summaries for 4 cuts, the cuts and why, CT ridge spacing and width. |
| out/cal2_look/profiles_raw.npz | The raw profiles. |
| out/cal2_look/ridge_follow.json | Per column: the CT ridge followers F1 (from the seed), F2 (back from the trace) and F3 (back from w035), and how dark each surface's path is. |
| out/cal2_look/split_check.json | Per column and row: backward follows from the trace, w035 and w034, and which mesh shares the trace's ridge. |
| out/cal2_look/png/cal2_step2.png, cal2_step5.png, cal2_step8.png | Ribbons along the trace front. Each image column is the CT along that point's own normal, up is outward. Red + is the trace, cyan w034, yellow w035, magenta w036. 1.13 µm per pixel. |
| out/cal2_look/png/cal2_col15.png | Trace column 15 through steps 0 to 8, with F1 (green), F2 (pink) and F3 (orange). |
| out/cal2_look/png/figures.json | Figure metadata. |
| out/cal2_look/run_console.txt, ridge_console.txt, split_console.txt | Console output. run1_console.txt is from a superseded first version (±150 µm lines); ignore it. |
| scripts/cal2_look.py, cal2_look_png.py, cal2_look_ridge.py, cal2_look_split.py | The scripts. Run them in the WSL venv with `python -I`. |
