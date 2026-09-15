"""Convert an Ansys Maxwell field-to-grid export (.fld, vector J on a regular grid) to a compact .npz.

Usage: python yalaz/fld_to_npz.py DBS_export.fld r0_ansys.npz

The grid is rebuilt from the coordinate columns, so rows can be in any order and the header is not
trusted. Points outside the saline are exported as 'Nan' by Ansys and become zero current.
"""
import sys
import numpy as np
import pandas as pd


def convert(src, dst):
    d = pd.read_csv(src, sep=r'\s+', skiprows=2, header=None, na_values=['Nan', 'NaN', 'nan']).to_numpy(np.float64)
    xyz = np.rint(d[:, :3] * 1e3).astype(int)                     # integer mm
    lo = xyz.min(0); shape = xyz.max(0) - lo + 1
    J = np.zeros((*shape, 3), np.float32)
    idx = xyz - lo
    J[idx[:, 0], idx[:, 1], idx[:, 2]] = np.nan_to_num(d[:, 3:6])
    axes = [(lo[i] + np.arange(shape[i])) * 1e-3 for i in range(3)]  # metres
    np.savez_compressed(dst, xs=axes[0], ys=axes[1], zs=axes[2], Jx=J[..., 0], Jy=J[..., 1], Jz=J[..., 2])
    mag = np.linalg.norm(J, axis=-1); i = np.unravel_index(mag.argmax(), mag.shape)
    print(f'grid {tuple(shape)}, nonzero voxels {int((mag > 0).sum())}, max |J| {mag.max():.0f} A/m2 at '
          f'({axes[0][i[0]] * 1e3:.0f}, {axes[1][i[1]] * 1e3:.0f}, {axes[2][i[2]] * 1e3:.0f}) mm -> {dst}')


if __name__ == '__main__':
    convert(sys.argv[1], sys.argv[2])
