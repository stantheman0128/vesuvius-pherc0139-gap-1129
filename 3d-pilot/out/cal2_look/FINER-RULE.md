# FINER-RULE: reading A or B at the cal2 seed, from level-1 CT

Written before any level-1 or level-2 voxel was fetched or read. The rule stays fixed once written. If it proves unworkable, the report says so and adds a second reading under a separate label.

Readings (from CAL2-LOOK.md):
- **A.** At the seed of columns 11 to 17, the single level-3 ridge is two sheets touching, w034's and w035's.
- **B.** That ridge is w035's sheet alone, and nothing sits at the w034 mesh's crossing.

## Data

- **Fetched chunks:** the 14 chunks listed in finer_ct_estimate.json, nothing else. Level 1 has 10 chunks (2.258 µm voxels) and level 2 has 4 (4.516 µm voxels). The level-1 set does not fill its 3x2x2 key box: the keys (28,48,59) and (28,48,60) are not fetched. Any sample that falls in an unfetched chunk is NaN, and a profile containing a NaN is "not checked".
- **Coordinates:** the same as cal2_look.CT. The level-L voxel index (z, y, x) comes from level-0 xyz divided by 2^L, and sampling is trilinear (map_coordinates, order 1). The possible centre/corner offset between levels is not corrected, the same as at level 3.
- **Level 3:** level-3 data (cal2's box, already on disk) goes through the same rule as a reference only.

## Profiles

- **Where:** at each trace point P[r, c] of the real-CT trace (cal2_look.load_trace), for rows r = 0, 1, 2 and columns c = 11 to 17.
- **Normal:** n is the local w035 outward normal, normal_at(P[r, c]) from cal2_look.load_meshes, the same normals as cal2_look.py.
- **Offsets:** s runs from -130 to +60 µm in 1.0 µm steps. 1.0 µm is cal2_look.py's DS. The brief said 2 µm, but cal2_look.py uses 1.0, and the brief asks for the same as cal2_look.py.
- **In-plane averaging:** each value is the mean CT over a small square in the sheet plane. The square has k x k points, with k = 2^(3-L): 4 at level 1, 2 at level 2, 1 at level 3. The points are one level-L voxel apart and centred on the line, so at every level the square is one level-3 voxel (9.03 µm) wide. This keeps the in-plane footprint of the level-3 data and adds resolution only along the normal.
  - Tangent axes: t1 is P[r, c+1] - P[r, c-1] made orthogonal to n and normalised. t2 = n x t1.
- **Mesh crossings:** s34 and s35 are the first crossings of w034 and w035 along the profile line (cal2_look Mesh.cross, nearest |t|).

## Noise and cut

- **Profile noise, sigma_prof (decides the cut), per level:** split each in-plane square into two halves along t1, columns i < k/2 and i >= k/2, and take d(s) = mean(half 1) - mean(half 2). Pool d over all samples of the 21 profiles (rows 0 to 2, columns 11 to 17) at that level.
  - sigma_half = 1.4826 x MAD(d) / 2, the noise of the full-square mean if the two halves are independent.
  - sigma_white = sigma_vox / k, the white-noise value for a mean of k^2 voxels.
  - sigma_prof = max(sigma_half, sigma_white). Correlated noise makes both estimates too low, so the larger one is used.
- **Voxel noise, sigma_vox (reported):** as in cal2_look.py, 1.4826 x MAD of neighbour differences / sqrt(2), along x and along z inside each fetched chunk, using only pairs where both voxels are nonzero. Take the larger of the x and z values.
- **Level 3:** k = 1, so there are no halves. sigma_prof is the cal2_look.py value, 4.19.
- **Ridge:** a peak of the profile (scipy find_peaks) whose prominence over the whole -130 to +60 µm profile is at least cut = 3 x sqrt(2) x sigma_prof, the same form as cal2_look's main cut noise3.

## Column class at the seed (row 0)

The window is W0 = [s34 - 15, s35 + 15] µm, and mid = (s34 + s35) / 2. Only ridges inside W0 count.

- **two sheets (A-like):** there are ridges R_up and R_lo with |R_up - s35| <= 15 µm, R_lo <= R_up - 12 µm and R_lo <= mid.
- **one sheet (B-like):** exactly one ridge in W0, with R > mid and |R - s34| > 15 µm.
- **other:** everything else. This covers no ridge in W0, one ridge on w034's half, or two ridges that do not meet the A-like conditions. It also covers a column with no w034 or w035 crossing, or a not-checked profile.

## Call (level 1 decides)

Level 1 has 7 columns (11 to 17).

- **A:** at least 5 columns are A-like and at most 1 is B-like.
- **B:** at least 5 columns are B-like and at most 1 is A-like.
- **undecided:** anything else.

Level 2 is classified the same way and reported beside level 1, but it does not decide. The same goes for level 3 (the premise check: it is expected to give mostly one sheet or other) and for level 1 at twice the cut (a robustness check).

## Rows 1 and 2: which sheet the trace is on

This part runs only for the level-1 columns that are A-like at row 0.

- **Path:** the column's own path, the polyline P[0..2, c], sampled every 4.516 µm. Each path point gets a profile as above, with normal_at at that point and t1 from the same row-interpolated neighbour columns.
- **Followers:** two followers start at row 0, one on R_up (w035's side) and one on R_lo (w034's side). At each step a follower moves to the nearest ridge within 4.516 + 2 µm of its last position. If there is no ridge in reach, it takes the brightest sample in reach and counts the step as lost.
- **Trace position at rows 1 and 2:** the trace's ridge is the ridge nearest s = 0. The trace is "on a sheet" if I(0) - p10 >= 0.5 x (I(ridge) - p10), with p10 taken over that profile.
- **Labels:**
  - "on w035's side" if only the up-follower is within 3 µm of the trace's ridge.
  - "on w034's side" if only the low-follower is.
  - "merged" if both followers are.
  - "neither" if neither is.
  - "off sheet" if the trace is not on a sheet.
- **Report:** the counts per row. No overall call is made from this part.

## Known limits, written down in advance

- A paper sheet is made of two crossed fibre layers. At fine voxels, one thick sheet can show two ridges. The rule cannot tell two touching sheets from one sheet's two layers when the two ridges sit where the w034 and w035 meshes put two sheets. In that case an A result means "two ridges are resolved where the meshes put two sheets", not more than that.
- Level 1 has the noise of finer voxels. The in-plane square averages 16 voxels, which reduces it but does not remove it.
- The followers at rows 1 and 2 can hop where ridges touch, as at level 3.
