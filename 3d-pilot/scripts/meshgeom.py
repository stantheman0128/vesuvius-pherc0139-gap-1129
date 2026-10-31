"""Grid-mesh helpers shared by the analysis scripts (same definitions as in analyze_copy.py)."""
import numpy as np


def densify(g, v, k=8):
    """Bilinear points inside every grid cell whose four corners are valid, plus the valid vertices."""
    q = v[:-1, :-1] & v[1:, :-1] & v[:-1, 1:] & v[1:, 1:]
    ii, jj = np.nonzero(q)
    t = (np.arange(k) + 0.5) / k
    a, b = np.meshgrid(t, t, indexing="ij")
    a, b = a.ravel()[None, :, None], b.ravel()[None, :, None]
    p00, p10 = g[ii, jj][:, None], g[ii + 1, jj][:, None]
    p01, p11 = g[ii, jj + 1][:, None], g[ii + 1, jj + 1][:, None]
    pts = (1 - a) * (1 - b) * p00 + a * (1 - b) * p10 + (1 - a) * b * p01 + a * b * p11
    return np.concatenate([pts.reshape(-1, 3), g[v]])


def cut(g, v, z):
    """Points where grid edges cross the plane z (level-0 voxels)."""
    pts = []
    for ax in (0, 1):
        a = g[:-1] if ax == 0 else g[:, :-1]
        b = g[1:] if ax == 0 else g[:, 1:]
        va = v[:-1] if ax == 0 else v[:, :-1]
        vb = v[1:] if ax == 0 else v[:, 1:]
        ok = va & vb & ((a[..., 2] - z) * (b[..., 2] - z) <= 0) & (a[..., 2] != b[..., 2])
        t = (z - a[..., 2][ok]) / (b[..., 2][ok] - a[..., 2][ok])
        pts.append(a[ok] + t[:, None] * (b[ok] - a[ok]))
    return np.concatenate(pts) if pts else np.zeros((0, 3))
