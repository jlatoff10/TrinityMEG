"""Directional-lead rotation from MEG recordings of DBS pulses, against an Ansys current export.

Pipeline (per recording, each with its own dev_head_t from the fif header):
  1. period-locked average of the stimulation pulse (exact stimulation period, no peak detection)
  2. signed peak per magnetometer -> measured map
  3. Ansys current export placed in head coordinates by the Ansys->head transform (--exact, recommended:
     the model origin is the fiducial midpoint, so no tip input is needed), or by an electrode tip
     position (--tip), then rotated about the shaft in steps of --roll-step
  4. correlation of the model map with the measured map at each roll; best roll, its noise-limited
     uncertainty, and the same uncertainty extrapolated to --full-duration seconds

Examples
  python analysis/rotation_fit.py --fif r0.fif r11.fif r45.fif --ansys data/r0_ansys.fld_compressed.npz --exact --i-sim <A>
  python analysis/rotation_fit.py --fif r0.fif --ansys data/r*_ansys.fld_compressed.npz --exact --family --i-sim <A>
  python analysis/rotation_fit.py --fif r0.fif --ansys data/J_2a_2b.npz \
        --fiducials-ansys LPAx LPAy LPAz NASx NASy NASz RPAx RPAy RPAz     (metres, Ansys frame)

Outputs (in --out, default analysis/results): <fif>_<ansys>_roll.png, rotation_results.csv, and
the per-sensor maps in <fif>_<ansys>.npz. With several fifs the CSV also lists rotations relative
to the first file, which cancels any constant offset in the roll convention.
"""
import argparse, os, re, sys, csv, numpy as np, mne
from scipy.spatial.transform import Rotation as Rot
from scipy.ndimage import maximum_filter
sys.path.insert(0, os.path.dirname(__file__))
from dbs_period_average import period_average
mne.set_log_level('ERROR')
MU0_4PI = 1e-7
TIP_HEAD = np.array([-0.03403, -0.04753, -0.08755])   # electrode tip, head coordinates (m)
AXIS_ANSYS = np.array([0., 1., 0.])                    # shaft direction in the Ansys frame (tip -> lid, toward the apex)
# Ansys -> head rotation. The Ansys model's origin is the LPA-RPA midpoint and its x-z plane is the fiducial plane;
# its axes are +x toward the left ear (head -x), +y toward the apex (head +z), +z toward the nose (head +y).
# Verified: the r0 export's contacts at Ansys (3, -41, -63) mm map to head (-3, -63, -41) mm, the SolidWorks tip.
R_NOMINAL = np.array([[-1., 0., 0.], [0., 0., 1.], [0., 1., 0.]])
R_FIT = R_NOMINAL
PLANAR_BASELINE = 0.0168                               # Neuromag planar gradiometer baseline (m)


def head_frame_from_fiducials(lpa, nas, rpa):
    """Rotation+translation taking Ansys-frame coordinates to the MNE head frame defined by 3 fiducials."""
    lpa, nas, rpa = map(np.asarray, (lpa, nas, rpa))
    o = (lpa + rpa) / 2; x = rpa - lpa; x /= np.linalg.norm(x)
    y = nas - o; y -= x * np.dot(y, x); y /= np.linalg.norm(y); z = np.cross(x, y)
    R = np.vstack([x, y, z])              # rows: head axes expressed in the Ansys frame
    return R, -R @ o


def load_ansys(npz, bin_mm=3):
    """Current elements (I*dl, A m) on a coarse grid, plus the contact-pair centre found from |J|."""
    z = np.load(npz); xs, ys, zs = z['xs'], z['ys'], z['zs']
    J = np.stack([z['Jx'], z['Jy'], z['Jz']], -1).astype(np.float64); mag = np.linalg.norm(J, axis=-1)
    lm = (mag == maximum_filter(mag, size=7)) & (mag > 0.1 * mag.max()); ii = np.array(np.nonzero(lm)).T
    pts = np.stack([xs[ii[:, 0]], ys[ii[:, 1]], zs[ii[:, 2]]], 1); w = mag[lm]
    contacts = (pts * w[:, None]).sum(0) / w.sum()
    Qnet = J.reshape(-1, 3).sum(0) * 1e-9
    k = bin_mm; nx, ny, nz = J.shape[:3]; pad = [(-s) % k for s in (nx, ny, nz)]
    Jc = np.pad(J, ((0, pad[0]), (0, pad[1]), (0, pad[2]), (0, 0)))
    Jc = Jc.reshape(-1, k, Jc.shape[1] // k, k, Jc.shape[2] // k, k, 3).sum(axis=(1, 3, 5)) * 1e-9
    cc = lambda a, p: np.pad(a, (0, p), mode='linear_ramp', end_values=a[-1] + p * 1e-3).reshape(-1, k).mean(1)
    G = np.stack(np.meshgrid(cc(xs, pad[0]), cc(ys, pad[1]), cc(zs, pad[2]), indexing='ij'), -1).reshape(-1, 3)
    Q = Jc.reshape(-1, 3); keep = np.linalg.norm(Q, axis=1) > 0
    # current between the contacts: flux of J through the mid-plane perpendicular to the net dipole
    qh = Qnet / np.linalg.norm(Qnet); X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    slab = np.abs((X - contacts[0]) * qh[0] + (Y - contacts[1]) * qh[1] + (Z - contacts[2]) * qh[2]) < 0.5e-3
    I_sim = abs(np.sum(J[slab] @ qh) * 1e-6)
    return G[keep], Q[keep], contacts, Qnet, I_sim


def sensors_head(info):
    """Magnetometer positions/normals and planar-gradiometer geometry in head coordinates."""
    T = info['dev_head_t']['trans']; R, t = T[:3, :3], T[:3, 3]
    m = mne.pick_types(info, meg='mag'); g = mne.pick_types(info, meg='grad')
    loc = np.array([info['chs'][i]['loc'] for i in range(len(info['chs']))])
    pos = loc[:, :3] @ R.T + t; ex = loc[:, 3:6] @ R.T; ez = loc[:, 9:12] @ R.T
    return dict(mag=m, grad=g, pos=pos, ex=ex, ez=ez)


def forward(G, Q, sens, chans, R_ah, tip_ans, tip_head, scale):
    """B.n at magnetometers, or (B.n(+d) - B.n(-d))/baseline at planar gradiometers."""
    P = (G - tip_ans) @ R_ah.T + tip_head; Qr = Q @ R_ah.T
    def bn(p, n):
        r = p - P; return MU0_4PI * np.sum(np.cross(Qr, r) / np.linalg.norm(r, axis=1)[:, None] ** 3 @ n) * scale
    out = np.zeros(len(chans))
    for j, i in enumerate(chans):
        if i in sens['mag']:
            out[j] = bn(sens['pos'][i], sens['ez'][i])
        else:
            d = sens['ex'][i] * PLANAR_BASELINE / 2
            out[j] = (bn(sens['pos'][i] + d, sens['ez'][i]) - bn(sens['pos'][i] - d, sens['ez'][i])) / PLANAR_BASELINE
    return out


def roll_matrix(R_ah, shaft_head, roll_deg):
    return (Rot.from_rotvec(np.radians(roll_deg) * shaft_head) * Rot.from_matrix(R_ah)).as_matrix()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fif', nargs='+', required=True); ap.add_argument('--ansys', nargs='+', required=True)
    ap.add_argument('--tip', nargs=3, type=float, default=TIP_HEAD, help='electrode tip, head coords (m)')
    ap.add_argument('--fiducials-ansys', nargs=9, type=float, help='LPA NAS RPA in the Ansys frame (m): exact transform')
    ap.add_argument('--exact', action='store_true', help='Ansys origin = fiducial midpoint and axes as in the convention: place the model by the transform alone, ignore --tip')
    ap.add_argument('--current', type=float, default=7.5e-3, help='stimulation current (A)')
    ap.add_argument('--f0', type=float, default=130.0); ap.add_argument('--roll-step', type=float, default=5.0)
    ap.add_argument('--grads', action='store_true', help='also use the 204 planar gradiometers')
    ap.add_argument('--exclude', nargs='*', default=[], help='channels to drop')
    ap.add_argument('--weighted', action='store_true', help='whiten channels by their split-half noise before the fit')
    ap.add_argument('--stat', choices=['pearson', 'cosine'], default='pearson', help='pearson: mean removed (robust to a uniform offset); cosine: Yalaz-style, no mean removal')
    ap.add_argument('--full-duration', type=float, default=180.0, help='seconds, for the SNR extrapolation')
    ap.add_argument('--i-sim', type=float, default=2.17e-3, help='current between the contacts in the Ansys solve at 1 V (A). Default 2.17 mA = Maxwell surface integral of J.n over the 1 V face of the r0 model')
    ap.add_argument('--family', action='store_true', help='the exports are true rotations of the lead (angle parsed from the file name, e.g. r45_...): compare each recording with every export, refine each with a small numerical roll of +/- --window deg, report the best export and the refined absolute angle')
    ap.add_argument('--window', type=float, default=15.0, help='family mode: numerical roll refinement half-width (deg)')
    ap.add_argument('--out', default='analysis/results')
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    tip_head = np.array(a.tip)
    if a.fiducials_ansys:
        f = np.array(a.fiducials_ansys).reshape(3, 3); R_ah, t_ah = head_frame_from_fiducials(*f); place_by_tip = False
    elif a.exact:
        R_ah, t_ah = R_NOMINAL, np.zeros(3); place_by_tip = False
    else:
        R_ah, t_ah = R_FIT, None; place_by_tip = True
    shaft_head = R_ah @ AXIS_ANSYS; shaft_head /= np.linalg.norm(shaft_head)
    rolls = np.arange(-a.window, a.window + 1e-9, a.roll_step) if a.family else np.arange(0, 360, a.roll_step)
    rows = []
    for npz in a.ansys:
        G, Q, contacts, Qnet, I_sim = load_ansys(npz)
        tip_ans = contacts if place_by_tip else None      # tip offset (<4 mm along the shaft) is irrelevant for the roll
        print(f'[{os.path.basename(npz)}] contacts at {np.round(contacts * 1e3, 1)} mm (Ansys), I_sim = {I_sim * 1e3:.2f} mA at 1 V, '
              f'net dipole dir {np.round(Qnet / np.linalg.norm(Qnet), 2)}')
        if a.i_sim: I_sim = a.i_sim
        scale = a.current / I_sim
        nominal = float(re.search(r'r(\d+)', os.path.basename(npz)).group(1)) if a.family and re.search(r'r(\d+)', os.path.basename(npz)) else 0.0
        for fif in a.fif:
            raw = mne.io.read_raw_fif(fif, preload=False); info = raw.info
            f0, tau, W, pk = period_average(raw, a.f0)
            sens = sensors_head(info)
            chans = list(sens['mag']) + (list(sens['grad']) if a.grads else [])
            chans = [c for c in chans if raw.ch_names[c] not in a.exclude]
            b = pk[chans]
            # split-half noise of the averaged map (same synthesis on each half)
            rawf = raw.copy().load_data().filter(60, None, method='iir', iir_params=dict(order=6, ftype='butter'))
            X = rawf.get_data()[chans]; tt = rawf.times; K = int(info['lowpass'] // f0); h = len(tt) // 2
            def synth(sl):
                Wh = np.zeros((len(chans), len(tau)))
                for k in range(1, K + 1):
                    c = (X[:, sl] * np.exp(-2j * np.pi * k * f0 * tt[sl])).mean(1)
                    Wh += 2 * np.real(c[:, None] * np.exp(2j * np.pi * k * f0 * tau)[None, :])
                return Wh
            dif = (synth(slice(0, h)) - synth(slice(h, None))) / 2
            sigma = np.sqrt((dif ** 2).mean(1))                  # per-channel noise of the averaged map at this duration
            if place_by_tip:
                th = tip_head
            else:
                th = R_ah @ contacts + t_ah                        # exact transform: contacts land where the model puts them
                tip_ans = contacts
            maps = np.array([forward(G, Q, sens, chans, roll_matrix(R_ah, shaft_head, r), tip_ans, th, scale) for r in rolls])
            wgt = 1 / sigma if a.weighted else np.ones_like(sigma)
            # goodness of fit as in Yalaz: cosine similarity without mean subtraction (amplitude-free)
            bw = b * wgt
            cosine = np.array([np.dot(mp * wgt, bw) / (np.linalg.norm(mp * wgt) * np.linalg.norm(bw)) for mp in maps])
            pearson = np.array([np.corrcoef(mp * wgt, bw)[0, 1] for mp in maps])
            corr = pearson if a.stat == 'pearson' else cosine
            other = rolls[int(np.argmax(cosine if a.stat == 'pearson' else pearson))]
            ib = int(np.argmax(corr)); best = rolls[ib]
            # refine with a parabola through the three points around the maximum
            if 0 < ib < len(rolls) - 1:
                y0, y1, y2 = corr[ib - 1], corr[ib], corr[ib + 1]; best = rolls[ib] + a.roll_step * 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
            mb = forward(G, Q, sens, chans, roll_matrix(R_ah, shaft_head, best), tip_ans, th, scale)
            amp = np.dot(mb * wgt, bw) / np.dot(mb * wgt, mb * wgt)
            # noise-limited uncertainty: sigma_theta = 1 / ||d(amp*model)/dtheta / sigma||, amplitude direction projected out
            dth = 1.0
            dm = (forward(G, Q, sens, chans, roll_matrix(R_ah, shaft_head, best + dth), tip_ans, th, scale)
                  - forward(G, Q, sens, chans, roll_matrix(R_ah, shaft_head, best - dth), tip_ans, th, scale)) / (2 * dth) * amp
            u = mb * amp / sigma; v = dm / sigma; v -= u * np.dot(v, u) / np.dot(u, u)
            sig_theta = 1 / np.linalg.norm(v)                     # degrees, at this recording's duration
            dur = raw.times[-1]; sig_full = sig_theta * np.sqrt(dur / a.full_duration)
            resid = b - amp * mb; chi = np.sqrt(np.mean((resid / sigma) ** 2)); pear = np.corrcoef(mb, b)[0, 1]
            tag = f'{os.path.splitext(os.path.basename(fif))[0]}_{os.path.splitext(os.path.basename(npz))[0]}'
            absolute = (nominal + best) % 360 if a.family else best
            print(f'[{tag}] f0 {f0:.4f} Hz, {len(chans)} ch, {dur:.0f} s | best roll {best:.1f} deg, corr {corr[ib]:.3f}, '
                  f'({a.stat}; the other statistic peaks at {other:.0f} deg)' + (f', export nominal {nominal:.0f} deg -> absolute {absolute:.1f} deg, corr at zero roll {corr[np.argmin(np.abs(rolls))]:.3f}' if a.family else '') + f', Pearson r {pear:.3f}, amp ratio {amp:.2f} | noise-limited sigma {sig_theta:.2f} deg now, {sig_full:.2f} deg at {a.full_duration:.0f} s | '
                  f'residual/noise {chi:.0f}x (systematic mismatch dominates when >> 1)')
            rows.append(dict(fif=fif, ansys=npz, nominal_deg=nominal, absolute_deg=absolute, f0=f0, duration_s=dur, best_roll_deg=best, gof=corr[ib], pearson_r=pear, amp_ratio=amp,
                             sigma_now_deg=sig_theta, sigma_full_deg=sig_full, residual_over_noise=chi))
            np.savez(os.path.join(a.out, tag + '.npz'), rolls=rolls, corr=corr, measured=b, model=amp * mb, sigma=sigma,
                     chans=np.array(raw.ch_names)[chans], best_roll=best)
            try:
                import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(5, 3.4)); ax.plot(rolls, corr, color='#1f5fa8', lw=2); ax.axvline(best, color='0.5', ls='--')
                ax.set(xlabel='roll about shaft (deg)', ylabel='goodness of fit', ylim=(-1, 1), title=tag); ax.grid(alpha=.3)
                fig.tight_layout(); fig.savefig(os.path.join(a.out, tag + '_roll.png'), dpi=130); plt.close(fig)
            except Exception as e:
                print('plot skipped:', e)
    if a.family:
        print('--- family summary (best export per recording) ---')
        for fif in a.fif:
            rr = sorted([r for r in rows if r['fif'] == fif], key=lambda r: -r['gof'])
            print(f"{os.path.basename(fif)}: best {os.path.basename(rr[0]['ansys'])} (gof {rr[0]['gof']:.3f}) -> angle {rr[0]['absolute_deg']:.1f} deg; "
                  + 'ranking: ' + ', '.join(f"r{r['nominal_deg']:.0f}:{r['gof']:.3f}" for r in rr))
    with open(os.path.join(a.out, 'rotation_results.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) + ['relative_to_first_deg']); w.writeheader()
        for npz in a.ansys:
            rr = [r for r in rows if r['ansys'] == npz]
            for r in rr:
                r['relative_to_first_deg'] = (r['best_roll_deg'] - rr[0]['best_roll_deg'] + 180) % 360 - 180; w.writerow(r)
    print('written', os.path.join(a.out, 'rotation_results.csv'))


if __name__ == '__main__':
    main()
