import numpy as np, pydicom
from pathlib import Path

files = [pydicom.dcmread(p) for p in Path("e:/noah_noah/Bucket_phantom/CT LOCALIZER/43059617/89861227").iterdir() if p.is_file()]
files = [d for d in files if "PixelData" in d]

# sort along the slice axis, not by filename
normal = np.cross(*np.array(files[0].ImageOrientationPatient, float).reshape(2, 3))
files.sort(key=lambda d: np.dot(np.array(d.ImagePositionPatient, float), normal))

vol = np.stack([d.pixel_array for d in files]).astype(np.float32)
vol = vol * float(files[0].RescaleSlope) + float(files[0].RescaleIntercept)  # → HU

dz = abs(np.dot(np.array(files[1].ImagePositionPatient, float) -
                np.array(files[0].ImagePositionPatient, float), normal))
spacing = (dz, *map(float, files[0].PixelSpacing))  # (z, row, col) mm

print(vol.shape, spacing, vol.min(), vol.max())



import matplotlib.pyplot as plt

i = len(vol) // 2
fig, ax = plt.subplots()
im = ax.imshow(vol[i], cmap="gray", vmin=-100, vmax=200); ax.axis("off")
title = ax.set_title(f"slice {i+1}/{len(vol)}")

def show(k):
    global i
    i = max(0, min(len(vol) - 1, k))
    im.set_data(vol[i]); title.set_text(f"slice {i+1}/{len(vol)}"); fig.canvas.draw_idle()

fig.canvas.mpl_connect("scroll_event", lambda e: show(i + (1 if e.button == "up" else -1)))
fig.canvas.mpl_connect("key_press_event",
    lambda e: show(i + {"up": 1, "right": 1, "down": -1, "left": -1}.get(e.key, 0)))




pos = np.array([d.ImagePositionPatient for d in files], float)   # origin of each slice, mm
R = np.array(files[0].ImageOrientationPatient, float).reshape(2, 3)  # row/col direction cosines
dr, dc = map(float, files[0].PixelSpacing)

def fmt(x, y):
    c, r = int(round(x)), int(round(y))
    if not (0 <= r < vol.shape[1] and 0 <= c < vol.shape[2]):
        return ""
    p = pos[i] + R[1] * c * dc + R[0] * r * dr     # DICOM LPS coordinates
    return f"x={p[0]:.1f} y={p[1]:.1f} z={p[2]:.1f} mm  [{vol[i, r, c]:.0f} HU]"

ax.format_coord = fmt




plt.show()