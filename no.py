# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 15:28:31 2026

@author: trini
"""


import pandas as pd
import numpy as np
import re

df = pd.read_csv('C:/Users/nebula/Desktop/r0_2a_3a.fld', sep=r'\s+', comment='#')


filepath = "C:/Users/nebula/Desktop/r0_2a_3a.fld"
with open(filepath, "r") as f:
    lines = f.readlines()
    print(lines[0:10])
    
    data_lines = lines[3:100]
    cleaned = "\n".join(
        re.sub(r"\bNan\b", "nan", line, flags=re.IGNORECASE) for line in data_lines
    )
    data = np.fromstring(cleaned, sep=" ")
    data = data.reshape(-1, 6)
    print(data[0:20])
    data[:, :3] = data[:, :3]*1000
    print(data[0:20])
    data[:, :3] = np.rint(data[:, :3])
    print(data[0:20])
    
    






import re
import numpy as np
 
MU0 = 4 * np.pi * 1e-7  # H/m
 
 
# ---------------------------------------------------------------------
# 1. Parse the .fld file: header (grid bounds/spacing) + data rows
# ---------------------------------------------------------------------
def parse_fld(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    print('total lines in file=',len(lines))
    print('last line=',lines[-1])
    header1 = lines[0]  # "Grid Output Min: [...] Max: [...]"
    header2 = lines[1]  # "Grid Size: [...] Unit: "mm""
 
    nums = lambda s: [float(x) for x in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)]
 
    min_max_vals = nums(header1)
    grid_min = np.array(min_max_vals[0:3])   # in mm
    grid_max = np.array(min_max_vals[3:6])   # in mm
 
    #size_vals = nums(header2)
    grid_size = np.array(min_max_vals[6:9])     # spacing in mm
 
    unit_match = re.search(r'Unit:\s*"(\w+)"', header2)
    unit = unit_match.group(1) if unit_match else "mm"
    unit_scale = {"mm": 1e-3, "m": 1.0, "cm": 1e-2}.get(unit, 1e-3)
 
    # data starts after the 3rd line (header1, header2, column-label line)
    data_lines = lines[2:]
    # replace Nan/NaN with 'nan' so numpy can parse it
    cleaned = "\n".join(
        re.sub(r"\bNan\b", "nan", line, flags=re.IGNORECASE) for line in data_lines
    )
    data = np.fromstring(cleaned, sep=" ")
    data = data.reshape(-1, 6)  # X, Y, Z, Jx, Jy, Jz
    print('length of data:',len(data))
    print('first data line = ',data[0])
    print('last data line = ',data[-1])
    
    data[:, :3] = data[:, :3]*1000
    print(data[0:20])
    data[:, :3] = np.rint(data[:, :3])
 
    nx = int(round((grid_max[0] - grid_min[0]) / grid_size[0])) + 1
    ny = int(round((grid_max[1] - grid_min[1]) / grid_size[1])) + 1
    nz = int(round((grid_max[2] - grid_min[2]) / grid_size[2])) + 1
    
    print(nx,ny,nz)
 
    # expected = nx * ny * nz
    # if data.shape[0] != expected:
    #     raise ValueError(
    #         f"Parsed {data.shape[0]} rows but expected {expected} "
    #         f"({nx} x {ny} x {nz}) -- check ordering/header parsing."
    #     )
 
    # reshape assuming Z fastest, then Y, then X (row-major / C order)
    J = data[:, 3:6].reshape(nx, ny, nz, 3)
    J = np.nan_to_num(J, nan=0.0)  # no current outside modeled region
 
    Jx = J[..., 0]
    Jy = J[..., 1]
    Jz = J[..., 2]
 
    xs = (grid_min[0] + np.arange(nx) * grid_size[0]) * unit_scale
    ys = (grid_min[1] + np.arange(ny) * grid_size[1]) * unit_scale
    zs = (grid_min[2] + np.arange(nz) * grid_size[2]) * unit_scale
    spacing = grid_size * unit_scale  # (dx, dy, dz) in meters
 
    return (xs, ys, zs), tuple(spacing), (Jx, Jy, Jz)
 
 
# ---------------------------------------------------------------------
# 2. Biot-Savart kernel K(r) = r / |r|^3 (double-sized for linear conv)
# ---------------------------------------------------------------------
def build_kernel(grid_shape, spacing):
    nx, ny, nz = grid_shape
    dx, dy, dz = spacing
 
    kx = (np.arange(2 * nx - 1) - (nx - 1)) * dx
    ky = (np.arange(2 * ny - 1) - (ny - 1)) * dy
    kz = (np.arange(2 * nz - 1) - (nz - 1)) * dz
 
    Kx, Ky, Kz = np.meshgrid(kx, ky, kz, indexing="ij")
    r2 = Kx**2 + Ky**2 + Kz**2
    r2[r2 == 0] = np.inf
    r3 = r2 * np.sqrt(r2)
 
    return Kx / r3, Ky / r3, Kz / r3
 
 
def fft_convolve_full(a, b):
    """Linear (non-circular) 3D convolution via zero-padded FFT."""
    out_shape = np.array(a.shape) + np.array(b.shape) - 1
    # pad each axis up to a fast FFT length (power of 2 here for simplicity)
    fshape = [int(2 ** np.ceil(np.log2(s))) for s in out_shape]
 
    A = np.fft.rfftn(a, fshape)
    B = np.fft.rfftn(b, fshape)
    conv = np.fft.irfftn(A * B, fshape)
 
    slices = tuple(slice(0, s) for s in out_shape)
    return conv[slices]
 
 
# ---------------------------------------------------------------------
# 3. Compute B from the 6 cross-product convolutions
# ---------------------------------------------------------------------
def compute_B(J_components, kernel, spacing, grid_shape):
    Jx, Jy, Jz = J_components
    Rx, Ry, Rz = kernel
    dx, dy, dz = spacing
    dV = dx * dy * dz
    nx, ny, nz = grid_shape
 
    JyKz = fft_convolve_full(Jy, Rz)
    JzKy = fft_convolve_full(Jz, Ry)
    JzKx = fft_convolve_full(Jz, Rx)
    JxKz = fft_convolve_full(Jx, Rz)
    JxKy = fft_convolve_full(Jx, Ry)
    JyKx = fft_convolve_full(Jy, Rx)
 
    Bx_full = JyKz - JzKy
    By_full = JzKx - JxKz
    Bz_full = JxKy - JyKx
 
    start = (nx - 1, ny - 1, nz - 1)
    sl = tuple(slice(start[i], start[i] + grid_shape[i]) for i in range(3))
 
    Bx = MU0 / (4 * np.pi) * Bx_full[sl] * dV
    By = MU0 / (4 * np.pi) * By_full[sl] * dV
    Bz = MU0 / (4 * np.pi) * Bz_full[sl] * dV
    return Bx, By, Bz
 
 
# ---------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------
if __name__ == "__main__":
    filepath = "C:/Users/nebula/Desktop/r0_2a_3a.fld"
 
    (xs, ys, zs), spacing, J_components = parse_fld(filepath)
    grid_shape = (len(xs), len(ys), len(zs))
    print("Grid shape:", grid_shape, "-- spacing (m):", spacing)

    kernel = build_kernel(grid_shape, spacing)
    Bx, By, Bz = compute_B(J_components, kernel, spacing, grid_shape)

    cx, cy, cz = grid_shape[0] // 2, grid_shape[1] // 2, grid_shape[2] // 2
    print("Sample B at grid center (T):", Bx[cx, cy, cz], By[cx, cy, cz], Bz[cx, cy, cz])

    np.savez("B_field_result_r0_2a_3a.npz", xs=xs, ys=ys, zs=zs, Bx=Bx, By=By, Bz=Bz)
    print("Saved B_field_result_r0_2a_3a.npz")


# ----------------------------------------------------------------------
# measuring at coordinates

import numpy as np

raw = """
-0.1066 0.0464 -0.0604
-0.1020 0.0631 -0.0256
-0.1085 0.0302 -0.0266
-0.1099 0.0131 -0.0627
-0.1074 0.0329 0.0080
-0.0989 0.0403 0.0413
-0.1011 0.0044 0.0408
-0.1083 -0.0011 0.0071
-0.0861 0.0988 0.0090
-0.0887 0.0757 0.0412
-0.0702 0.0758 0.0707
-0.1003 0.0659 0.0081
-0.0808 0.0413 0.0720
-0.0526 0.0406 0.0952
-0.0537 0.0059 0.0969
-0.0829 0.0062 0.0728
-0.0637 0.1254 0.0136
-0.0332 0.1397 0.0174
-0.0337 0.1274 0.0485
-0.0672 0.1089 0.0443
-0.0358 0.1048 0.0750
0.0001 0.0775 0.0967
-0.0184 0.0440 0.1063
-0.0368 0.0753 0.0922
-0.0185 0.0105 0.1096
0.0186 0.0105 0.1096
0.0186 -0.0233 0.1059
-0.0185 -0.0237 0.1058
0.0001 0.1445 0.0187
0.0001 0.1316 0.0500
0.0331 0.1397 0.0173
0.0638 0.1253 0.0135
0.0671 0.1088 0.0444
0.0338 0.1273 0.0486
0.0001 0.1093 0.0771
0.0358 0.1048 0.0750
0.0368 0.0752 0.0923
0.0184 0.0442 0.1062
0.0525 0.0406 0.0953
0.0809 0.0413 0.0721
0.0828 0.0061 0.0728
0.0535 0.0062 0.0970
0.0862 0.0986 0.0089
0.1003 0.0660 0.0082
0.0887 0.0757 0.0412
0.0699 0.0758 0.0709
0.0989 0.0404 0.0413
0.1074 0.0329 0.0081
0.1083 -0.0011 0.0068
0.1010 0.0044 0.0410
0.1020 0.0630 -0.0260
0.1065 0.0469 -0.0600
0.1098 0.0131 -0.0622
0.1083 0.0301 -0.0262
-0.1088 -0.0032 -0.0284
-0.1017 -0.0360 -0.0281
-0.0951 -0.0524 -0.0623
-0.1068 -0.0205 -0.0625
-0.1017 -0.0339 0.0056
-0.0952 -0.0308 0.0391
-0.0781 -0.0628 0.0394
-0.0866 -0.0640 0.0055
-0.0758 -0.0797 -0.0621
-0.0861 -0.0660 -0.0282
-0.0632 -0.0905 -0.0278
-0.0489 -0.0994 -0.0621
-0.0786 -0.0287 0.0696
-0.0518 -0.0277 0.0927
-0.0181 -0.0542 0.0923
-0.0552 -0.0627 0.0707
-0.0513 -0.0861 0.0397
-0.0335 -0.1033 0.0062
-0.0331 -0.1051 -0.0278
-0.0636 -0.0884 0.0060
-0.0186 -0.0801 0.0690
0.0186 -0.0802 0.0689
0.0169 -0.0972 0.0397
-0.0170 -0.0972 0.0397
0.0000 -0.1086 0.0070
0.0000 -0.1106 -0.0275
0.0170 -0.1098 -0.0618
-0.0171 -0.1098 -0.0619
0.0517 -0.0276 0.0928
0.0786 -0.0284 0.0698
0.0553 -0.0628 0.0706
0.0182 -0.0542 0.0922
0.0513 -0.0861 0.0397
0.0637 -0.0884 0.0062
0.0330 -0.1051 -0.0277
0.0333 -0.1034 0.0063
0.0952 -0.0306 0.0392
0.1017 -0.0338 0.0058
0.0866 -0.0639 0.0056
0.0781 -0.0629 0.0394
0.0630 -0.0906 -0.0276
0.0861 -0.0661 -0.0280
0.0757 -0.0798 -0.0620
0.0488 -0.0994 -0.0621
0.1087 -0.0033 -0.0280
0.1068 -0.0206 -0.0621
0.0951 -0.0524 -0.0620
0.1017 -0.0361 -0.0278

"""

sensors = np.fromstring(raw, sep=" ").reshape(-1, 3)

from scipy.interpolate import RegularGridInterpolator

interp_Bx = RegularGridInterpolator((xs, ys, zs), Bx)
interp_By = RegularGridInterpolator((xs, ys, zs), By)
interp_Bz = RegularGridInterpolator((xs, ys, zs), Bz)


Bx_s = interp_Bx(sensors)
By_s = interp_By(sensors)
Bz_s = interp_Bz(sensors)

B_sensors = np.column_stack([Bx_s, By_s, Bz_s])


print(xs.min(), xs.max())
print(ys.min(), ys.max())
print(zs.min(), zs.max())
print(spacing)
















 