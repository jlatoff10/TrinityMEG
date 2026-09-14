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
    ap.add_argument('--method', choices=['harmonic', 'scan'], default='harmonic', help='harmonic: angle from the cos/sin components of the model with the roll-invariant part free (default); scan: best whole-map correlation over roll')
    ap.add_argument('--map', choices=['template', 'peak'], default='template', help='per-channel amplitude: projection on the common waveform (default) or signed peak (Yalaz)')
    ap.add_argument('--chunk', type=float, default=10.0, help='chunk length (s) for the drift-tolerant average')
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
        if a.i_sim: I_sim = a.i_sim
        print(f'[{os.path.basename(npz)}] contacts at {np.round(contacts * 1e3, 1)} mm (Ansys), I_sim = {I_sim * 1e3:.2f} mA at 1 V, '
              f'net dipole dir {np.round(Qnet / np.linalg.norm(Qnet), 2)}')
        scale = a.current / I_sim
        nominal = float(re.search(r'r(\d+)', os.path.basename(npz)).group(1)) if a.family and re.search(r'r(\d+)', os.path.basename(npz)) else 0.0
        for fif in a.fif:
            raw = mne.io.read_raw_fif(fif, preload=False); info = raw.info
            dur = raw.times[-1]
            f0, tau, W, pk, ex = period_average(raw, a.f0, chunk_s=min(a.chunk, dur / 4), return_extras=True)
            print(f'[{os.path.basename(fif)}] {ex["n_chunks"]} chunks of {min(a.chunk, dur / 4):.0f} s: f0 {ex["chunk_f0"].min():.4f}..{ex["chunk_f0"].max():.4f} Hz '
                  f'(spread {np.ptp(ex["chunk_f0"]) * 1e3:.1f} mHz = {np.ptp(ex["chunk_f0"]) / f0 * 1e6:.0f} ppm), phase shifts up to {np.abs(ex["chunk_shift_ms"]).max():.2f} ms; '
                  f'common waveform explains {ex["template_var_fraction"] * 100:.0f}% of the magnetometer variance')
            sens = sensors_head(info); loc_dev = np.array([ch['loc'][:3] for ch in info['chs']]); loc_nrm = np.array([ch['loc'][9:12] for ch in info['chs']])
            chans = list(sens['mag']) + (list(sens['grad']) if a.grads else [])
            chans = [c for c in chans if raw.ch_names[c] not in a.exclude]
            b = (ex['amp'] if a.map == 'template' else pk)[chans]
            sigma = ex['noise'][chans]                             # per-channel standard error from the chunk scatter
            if np.any(~np.isfinite(sigma)): sigma = np.full(len(chans), np.nanmedian(sigma) if np.any(np.isfinite(sigma)) else np.abs(b).max() * 0.01)
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
            # --- harmonic method: model(theta) = M0 + Mc cos(theta) + Ms sin(theta) + ...; fit the measured map with free
            #     coefficients for the roll-invariant part M0 (absorbs its model error) and the rotating part; theta = atan2(c, b)
            if a.method == 'harmonic' and not a.family:
                th_r = np.radians(rolls)
                M0 = maps.mean(0); Mc = 2 * np.mean(maps * np.cos(th_r)[:, None], 0); Ms = 2 * np.mean(maps * np.sin(th_r)[:, None], 0)
                M2c = 2 * np.mean(maps * np.cos(2 * th_r)[:, None], 0); M2s = 2 * np.mean(maps * np.sin(2 * th_r)[:, None], 0)
                Wd = np.diag(wgt)
                A = np.column_stack([M0, Mc, Ms, np.ones_like(M0)])
                coef, *_ = np.linalg.lstsq(Wd @ A, wgt * b, rcond=None)
                theta = np.degrees(np.arctan2(coef[2], coef[1])) % 360
                fit = A @ coef; r2 = 1 - np.sum((wgt * (b - fit)) ** 2) / np.sum((wgt * (b - np.average(b, weights=wgt ** 2))) ** 2)
                rot_power = np.sum((wgt * (coef[1] * Mc + coef[2] * Ms)) ** 2) / np.sum((wgt * fit) ** 2)
                # uncertainty from the linear fit covariance (noise-limited), propagated to the angle
                Aw = A / sigma[:, None]; cov = np.linalg.inv(Aw.T @ Aw)            # covariance of the coefficients given the per-channel noise
                Jt = np.array([0, -coef[2], coef[1], 0]) / (coef[1] ** 2 + coef[2] ** 2)
                sig_h = np.degrees(np.sqrt(Jt @ cov @ Jt))
                second = np.sum(M2c ** 2 + M2s ** 2) / np.sum(Mc ** 2 + Ms ** 2)
                print(f'  harmonic fit: theta = {theta:.1f} deg, noise-limited sigma {sig_h:.2f} deg, fit R2 {r2:.3f}, rotating part = {rot_power * 100:.0f}% of the fitted power, '
                      f'invariant-part amplitude ratio {coef[0]:.2f}, rotating-part amplitude {np.hypot(coef[1], coef[2]):.2f}, 2nd-harmonic/1st power in the model {second:.2f}')
                best = theta
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
            sig_full = sig_theta * np.sqrt(dur / a.full_duration)
            resid = b - amp * mb; chi = np.sqrt(np.mean((resid / sigma) ** 2)); pear = np.corrcoef(mb, b)[0, 1]
            tag = f'{os.path.splitext(os.path.basename(fif))[0]}_{os.path.splitext(os.path.basename(npz))[0]}'
            absolute = (nominal + best) % 360 if a.family else best
            print(f'[{tag}] f0 {f0:.4f} Hz, {len(chans)} ch, {dur:.0f} s | best roll {best:.1f} deg, corr {corr[ib]:.3f}, '
                  f'({a.stat}; the other statistic peaks at {other:.0f} deg)' + (f', export nominal {nominal:.0f} deg -> absolute {absolute:.1f} deg, corr at zero roll {corr[np.argmin(np.abs(rolls))]:.3f}' if a.family else '') + f', Pearson r {pear:.3f}, amp ratio {amp:.2f} | noise-limited sigma {sig_theta:.2f} deg now, {sig_full:.2f} deg at {a.full_duration:.0f} s | '
                  f'residual/noise {chi:.0f}x (systematic mismatch dominates when >> 1)')
            rows.append(dict(fif=fif, ansys=npz, nominal_deg=nominal, absolute_deg=absolute, f0=f0, duration_s=dur, best_roll_deg=best, gof=corr[ib], pearson_r=pear, amp_ratio=amp,
                             sigma_now_deg=sig_theta, sigma_full_deg=sig_full, residual_over_noise=chi))
            np.savez(os.path.join(a.out, tag + '.npz'), rolls=rolls, corr=corr, measured=b, model=amp * mb, sigma=sigma,
                     chans=np.array(raw.ch_names)[chans], best_roll=best, chunk_amps=ex['chunk_amps'][:, chans], chunk_starts_s=ex['chunk_starts_s'],
                     chunk_f0=ex['chunk_f0'], chunk_shift_ms=ex['chunk_shift_ms'], pos_dev=sens['pos_dev'][chans] if 'pos_dev' in sens else loc_dev[chans], nrm_dev=loc_nrm[chans])
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
    def set_angle(fif):
        m = re.search(r'_r(\d+)_', os.path.basename(fif)); return float(m.group(1)) if m else float('nan')
    with open(os.path.join(a.out, 'rotation_results.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) + ['relative_to_first_deg', 'set_angle_deg', 'set_relative_to_first_deg', 'error_deg']); w.writeheader()
        for npz in a.ansys:
            rr = [r for r in rows if r['ansys'] == npz]
            for r in rr:
                r['relative_to_first_deg'] = (r['best_roll_deg'] - rr[0]['best_roll_deg'] + 180) % 360 - 180
                r['set_angle_deg'] = set_angle(r['fif']); r['set_relative_to_first_deg'] = set_angle(r['fif']) - set_angle(rr[0]['fif'])
                # measured relative rotation vs the set one, checked for either rotation sense (the sign convention is not calibrated)
                e_plus = (r['relative_to_first_deg'] - r['set_relative_to_first_deg'] + 180) % 360 - 180
                e_minus = (r['relative_to_first_deg'] + r['set_relative_to_first_deg'] + 180) % 360 - 180
                r['error_deg'] = e_plus if abs(e_plus) <= abs(e_minus) else e_minus
                w.writerow(r)
        print('summary (set angle -> measured rotation relative to the first file):')
        for npz in a.ansys:
            for r in [x for x in rows if x['ansys'] == npz]:
                print(f"  {os.path.basename(r['fif'])}: set {r['set_angle_deg']:.0f} deg, measured relative {r['relative_to_first_deg']:+.1f} deg, error {r['error_deg']:+.1f} deg, corr {r['gof']:.3f}")
    print('written', os.path.join(a.out, 'rotation_results.csv'))


if __name__ == '__main__':
    main()
