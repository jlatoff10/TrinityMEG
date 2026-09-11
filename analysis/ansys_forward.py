"""Biot-Savart field of an Ansys Maxwell current export at the MEG magnetometers.

Usage:  python analysis/ansys_forward.py data/J_2a_3a.npz out_chest.npz

The .npz comes from the corrected .fld parser (C order: z varies fastest in the file). The current
distribution is placed in head coordinates by (1) the electrode tip position TIP_HEAD, (2) the
tip position and shaft axis in the Ansys frame, (3) a rotation taking the Ansys shaft axis to the
shaft direction in head coordinates, and (4) a roll about the shaft. The field is scaled from the
simulated current (1 V drive) to the stimulation current.
"""
import sys, numpy as np
from scipy.spatial.transform import Rotation as Rot
MU0 = 4 * np.pi * 1e-7

TIP_HEAD = np.array([-0.03403, -0.04753, -0.08755])      # electrode tip, head coordinates (m)
CONTACTS_ANSYS = np.array([0.007, -0.057, 0.036])        # centre of the active pair in the Ansys frame (m)
AXIS_ANSYS = np.array([0., 1., 0.])                      # shaft direction in the Ansys frame
TIP_ANSYS = CONTACTS_ANSYS + 0.00375 * AXIS_ANSYS         # tip is 3.75 mm from the 2a-3a pair centre; +y side fits the data
I_SIM, I_STIM = 6.0e-3, 4.0e-3                            # current through the inter-contact plane at 1 V; stimulation current
LOWPASS_ATTEN = 2 * 1650 * 60e-6                          # peak of a 60 us pulse after the 1650 Hz acquisition low-pass (~0.20)

def load_elements(npz, bin_mm=3):
    z = np.load(npz); xs, ys, zs = z['xs'], z['ys'], z['zs']
    J = np.stack([z['Jx'], z['Jy'], z['Jz']], -1).astype(np.float64)
    k = bin_mm; nx, ny, nz = J.shape[:3]; pad = [(-s) % k for s in (nx, ny, nz)]
    Jc = np.pad(J, ((0, pad[0]), (0, pad[1]), (0, pad[2]), (0, 0)))
    Jc = Jc.reshape(-1, k, Jc.shape[1] // k, k, Jc.shape[2] // k, k, 3).sum(axis=(1, 3, 5)) * 1e-9   # I*dl per bin (A m)
    cc = lambda a, p: np.pad(a, (0, p), mode='linear_ramp', end_values=a[-1] + p * 1e-3).reshape(-1, k).mean(1)
    G = np.stack(np.meshgrid(cc(xs, pad[0]), cc(ys, pad[1]), cc(zs, pad[2]), indexing='ij'), -1).reshape(-1, 3)
    Q = Jc.reshape(-1, 3); keep = np.linalg.norm(Q, axis=1) > 0
    return G[keep], Q[keep]

# Rotation Ansys -> head from the fit to the chest recording (rotation + 15 mm tip refinement, corr 0.765).
R_FIT = np.array([[ 0.474464, -0.00193 , -0.880273],
 [ 0.867049, -0.171673,  0.467712],
 [-0.152022, -0.985152, -0.079779]])
TIP_FIT = np.array([-0.025614, -0.049458, -0.099917])

def rotation():
    return R_FIT

def forward(G, Q, pos, nrm, R=None, tip_head=None):
    tip_head = TIP_FIT if tip_head is None else tip_head
    R = rotation() if R is None else R
    P = (G - TIP_ANSYS) @ R.T + tip_head; Qr = Q @ R.T; out = np.zeros(len(pos))
    for i in range(len(pos)):
        r = pos[i] - P
        out[i] = np.sum(np.cross(Qr, r) / np.linalg.norm(r, axis=1)[:, None] ** 3 @ nrm[i])
    return MU0 / (4 * np.pi) * out * I_STIM / I_SIM

if __name__ == '__main__':
    G, Q = load_elements(sys.argv[1])
    d = np.load(sys.argv[2]); T = d['dev_head_t']; m = d['mags']
    pos = d['pos_dev'][m] @ T[:3, :3].T + T[:3, 3]; nrm = d['nrm_dev'][m] @ T[:3, :3].T
    b_model = forward(G, Q, pos, nrm); b_meas = d['pk'][m]
    good = np.array([n != 'MEG0731' for n in d['names'][m]])            # noisy channel in the chest recording
    c = np.corrcoef(b_model[good], b_meas[good])[0, 1]; a = np.dot(b_model[good], b_meas[good]) / np.dot(b_model[good], b_model[good])
    print(f'model peak {np.abs(b_model).max():.2e} T (before low-pass attenuation), measured {np.abs(b_meas).max():.2e} T')
    print(f'correlation {c:.3f}; measured/model amplitude ratio {a:.2f} (expected ~{LOWPASS_ATTEN:.2f} from the 1650 Hz low-pass)')
