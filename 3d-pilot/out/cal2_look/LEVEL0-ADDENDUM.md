# LEVEL0-ADDENDUM: reading the cal2 seed at level 0 with FINER-RULE.md

Written before any level-0 byte was fetched or read. It stays fixed once its hash is saved to logs/level0_addendum_sha256.json.

Why: Stan asked whether a finer CT could change the result. Levels 3, 2 and 1 show one ridge about 60 µm wide at half height at all seven seed columns. Level 0 (1.129 µm voxels) is the finest level of the scan.

## The rule is unchanged

Level 0 is read with out/cal2_look/FINER-RULE.md as it stands (sha256 c7d471962f91e3f6e36e42a18ba1bf3efe3760f874ba55c13f58a23b0b2feea8). Nothing in it is edited, and every distance stays in µm, not voxels.

- **Places:** the same seven seed columns, 11 to 17, at rows 0, 1 and 2 of the real-CT trace, with the same w035 normals (normal_at) and the same tangent axes t1 and t2.
- **Profiles:** s from -130 to +60 µm in 1.0 µm steps. The in-plane square follows the rule's k = 2^(3-L), so k = 8 at level 0: 8 x 8 points 1.129 µm apart, one level-3 voxel (9.03 µm) wide, as at the other levels. Sampling is trilinear (map_coordinates, order 1) at level-0 xyz / 2^0.
- **Noise and cut:** sigma_half from the two halves of the square along t1 (i < 4 and i >= 4), sigma_white = sigma_vox / 8, sigma_prof = the larger of the two, cut = 3 x sqrt(2) x sigma_prof. sigma_vox is measured inside the fetched level-0 chunks as the rule says.
- **Classes:** the same window W0 = [s34 - 15, s35 + 15] µm and the same thresholds: |R_up - s35| <= 15 µm, R_lo <= R_up - 12 µm and R_lo <= mid for two sheets; one ridge with R > mid and |R - s34| > 15 µm for one sheet; everything else is other.
- **Call thresholds:** A needs at least 5 A-like and at most 1 B-like column, B needs at least 5 B-like and at most 1 A-like, anything else is undecided.
- **Also run at level 0, as the rule runs them at level 1:** the same classes at twice the cut, as a robustness check, and the rows 1 and 2 followers for any level-0 column that is A-like at row 0 (path step 4.516 µm, reach 4.516 + 2 µm, match 3 µm, unchanged).

## How level 0 is reported

- Level 0's verdict is its own line, next to level 1's. FINER-RULE.md's call (level 1 decides) is not changed or replaced by it.
- The level-0 line is reported whatever it says. If it is A, the earlier statement that levels 2 and 1 show one ridge stays true, but the claim that finer CT shows one ridge is overturned at level 0.

## Data

- **Chunks:** the level-0 keys of a new HTTP HEAD estimate made the same way as scripts/cal2_look_finer_estimate.py (the patch of trace rows 0 to 2 and columns 11 to 17, normals from -130 to +60 µm, the same margins). Nothing else is fetched.
- **Budget:** 160,000,000 bytes. The ledger (logs/fetch_log.jsonl) may reach 1,633,529,433 bytes and no more. If the estimate is over budget, nothing is fetched.
- **Where:** a new, empty directory, downloads/ct_level0/, as 0/z/y/x. A sample that falls in an unfetched chunk is NaN, and a profile that contains a NaN is "not checked", as in the rule.

## Descriptive, outside the rule

These do not enter the classes or the call.

- **Width and edge:** the half-height width of the single seed ridge and its 10% to 90% edge rise on each flank, measured as scripts/cal2_look_finer_width.py measures them, next to the level 1, 2 and 3 values.
- **Dip check:** for each column at row 0, every local maximum inside W0 with its prominence and no cut, and the largest prominence of any maximum other than the main ridge. That number is the depth of the deepest dip next to a second peak, and it is given in grey levels and as a fraction of the level-0 cut. A second ridge in the rule's sense is one that passes the cut.
- **Figure:** column 15 at level 0, rows 0 to 2, drawn the same way as png/cal2_col15_L1.png.

## Known limits, written down in advance

- At 1.129 µm the two fibre layers of a single sheet may be resolved. As FINER-RULE.md says, two ridges where the meshes put two sheets would mean "two ridges are resolved there", not more than that.
- Level-0 voxels are noisier. The 64-voxel square reduces the noise but does not remove it.
