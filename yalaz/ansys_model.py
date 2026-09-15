"""Magnetic field at the MEG sensors from an Ansys Maxwell current-density export.

The export (converted to .npz by fld_to_npz.py) holds the current density J in the saline on a 1 mm
grid, from a DC conduction solve with 1 V between the two active contacts. Its field at a sensor is
the Biot-Savart sum over all current elements, scaled from the solve's current I_sim to the
stimulation current.

Coordinates: the Ansys model's origin is the midpoint between LPA and RPA and its x-z plane contains
the three fiducials, so it maps onto the MEG head frame by a fixed rotation:
    Ansys +x -> toward the left ear  (head -x)
    Ansys +y -> toward the apex      (head +z)   = the lead's shaft direction
    Ansys +z -> toward the nose      (head +y)
Sensors come in the device frame and are moved to the head frame with each recording's own
dev_head_t (its HPI fit), so the model needs no per-recording adjustment.
"""
import numpy as np
from scipy.ndimage import maximum_filter
from scipy.spatial.transform import Rotation
import mne

MU0_OVER_4PI = 1e-7
ANSYS_TO_HEAD = np.array([[-1., 0., 0.],
                          [0., 0., 1.],
                          [0., 1., 0.]])
SHAFT_HEAD = ANSYS_TO_HEAD @ np.array([0., 1., 0.])


def load_export(npz_path, bin_mm=3):
    """Current elements I*dl (A m) on a bin_mm grid, in head coordinates, plus the contact centre.

    Returns (positions (N,3), elements (N,3), contacts (3,)). Binning to 3 mm keeps the field at
    sensors more than a few cm away unchanged while making the sum 27 times cheaper.
    """
    z = np.load(npz_path)
    xs, ys, zs = z['xs'], z['ys'], z['zs']
    J = np.stack([z['Jx'], z['Jy'], z['Jz']], axis=-1).astype(np.float64)
    mag = np.linalg.norm(J, axis=-1)

    # the active contacts are where the current density peaks
    peaks = (mag == maximum_filter(mag, size=7)) & (mag > 0.1 * mag.max())
    ii = np.array(np.nonzero(peaks)).T
    pts = np.stack([xs[ii[:, 0]], ys[ii[:, 1]], zs[ii[:, 2]]], axis=1)
    contacts = (pts * mag[peaks][:, None]).sum(0) / mag[peaks].sum()

    # sum J over bin_mm cubes -> I*dl per cube (1 mm voxels: J * 1e-9 m^3)
    k = bin_mm
    pad = [(-s) % k for s in J.shape[:3]]
    Jp = np.pad(J, ((0, pad[0]), (0, pad[1]), (0, pad[2]), (0, 0)))
    elements = Jp.reshape(Jp.shape[0] // k, k, Jp.shape[1] // k, k, Jp.shape[2] // k, k, 3).sum(axis=(1, 3, 5)) * 1e-9
    centres = lambda a, p: np.pad(a, (0, p), mode='linear_ramp', end_values=a[-1] + p * 1e-3).reshape(-1, k).mean(1)
    grid = np.stack(np.meshgrid(centres(xs, pad[0]), centres(ys, pad[1]), centres(zs, pad[2]), indexing='ij'), -1)
    positions, elements = grid.reshape(-1, 3), elements.reshape(-1, 3)
    keep = np.linalg.norm(elements, axis=1) > 0
    return positions[keep] @ ANSYS_TO_HEAD.T, elements[keep] @ ANSYS_TO_HEAD.T, ANSYS_TO_HEAD @ contacts


def magnetometers_head(info):
    """Magnetometer channel indices, positions and normals in head coordinates."""
    T = info['dev_head_t']['trans']
    R, t = T[:3, :3], T[:3, 3]
    mags = mne.pick_types(info, meg='mag')
    loc = np.array([info['chs'][i]['loc'] for i in mags])
    return mags, loc[:, :3] @ R.T + t, loc[:, 9:12] @ R.T


def rolled(positions, elements, contacts, roll_deg):
    """Current distribution rotated by roll_deg about the shaft axis through the contacts."""
    R = Rotation.from_rotvec(np.radians(roll_deg) * SHAFT_HEAD).as_matrix()
    return (positions - contacts) @ R.T + contacts, elements @ R.T


def field_normal(positions, elements, sensor_pos, sensor_nrm, scale):
    """Normal field component at each sensor (T), Biot-Savart over the current elements."""
    out = np.zeros(len(sensor_pos))
    for i, (p, n) in enumerate(zip(sensor_pos, sensor_nrm)):
        r = p - positions
        out[i] = np.sum(np.cross(elements, r) / np.linalg.norm(r, axis=1)[:, None] ** 3 @ n)
    return MU0_OVER_4PI * out * scale
