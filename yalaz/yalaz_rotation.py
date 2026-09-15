"""Lead rotation from MEG recordings of DBS pulses, by Yalaz's method.

For each recording:
  1. average the stimulation pulse per magnetometer (dbs_average.py)          -> measured map
  2. drop magnetometers close to a known artifact (--exclude-near)
  3. rotate the Ansys current distribution about the lead's shaft in steps and, at each roll,
     correlate its field pattern at the same magnetometers with the measured map
  4. the roll with the highest correlation is the lead orientation (refined by a parabola through
     the three points around the maximum)
Rotations are reported relative to the first recording, and compared with the set angle parsed from
the file name (..._r<angle>_...), for either rotation sense.

Example
  python yalaz/yalaz_rotation.py --ansys data/r0_ansys.fld_compressed.npz \
      --exclude-near 0.095 -0.04 0.085 0.12 --out yalaz/results \
      --fif 260911/260911_r0_2a_2b_7p5_60_130_raw.fif 260911/260911_r2_2a_2b_7p5_60_130_raw.fif ...
"""
import argparse, csv, os, re
import numpy as np
import mne
from dbs_average import pulse_average
from ansys_model import load_export, magnetometers_head, rolled, field_normal

mne.set_log_level('ERROR')


def set_angle(path):
    m = re.search(r'_r(\d+)_', os.path.basename(path))
    return float(m.group(1)) if m else np.nan


def wrap(deg):
    return (deg + 180) % 360 - 180


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fif', nargs='+', required=True, help='recordings; the first is the reference')
    ap.add_argument('--ansys', required=True, help='Ansys export (.npz from fld_to_npz.py)')
    ap.add_argument('--current', type=float, default=7.5e-3, help='stimulation current (A)')
    ap.add_argument('--i-sim', type=float, default=2.17e-3, help='current in the 1 V Ansys solve (A)')
    ap.add_argument('--f0', type=float, default=130.0, help='nominal stimulation frequency (Hz)')
    ap.add_argument('--roll-step', type=float, default=1.0, help='roll scan step (deg)')
    ap.add_argument('--exclude-near', nargs=4, type=float, metavar=('X', 'Y', 'Z', 'R'),
                    help='drop magnetometers within R m of (X,Y,Z), device coordinates')
    ap.add_argument('--out', default='yalaz/results')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    positions, elements, contacts = load_export(a.ansys)
    scale = a.current / a.i_sim
    rolls = np.arange(0, 360, a.roll_step)
    rows = []
    for fif in a.fif:
        raw = mne.io.read_raw_fif(fif, preload=False)
        duration = raw.times[-1]
        f0, tau, W, amp, info = pulse_average(raw, a.f0, chunk_s=min(10.0, duration / 4))
        mags, pos, nrm = magnetometers_head(raw.info)
        keep = np.ones(len(mags), bool)
        if a.exclude_near:
            dev = np.array([raw.info['chs'][i]['loc'][:3] for i in mags])
            keep = np.linalg.norm(dev - np.array(a.exclude_near[:3]), axis=1) > a.exclude_near[3]
        measured = amp[mags][keep]
        pos, nrm = pos[keep], nrm[keep]

        # roll scan: correlation of the model pattern with the measured map
        corr = np.array([np.corrcoef(field_normal(*rolled(positions, elements, contacts, r), pos, nrm, scale), measured)[0, 1]
                         for r in rolls])
        i = int(np.argmax(corr)); best = rolls[i]
        if 0 < i < len(rolls) - 1:
            y0, y1, y2 = corr[i - 1], corr[i], corr[i + 1]
            best = rolls[i] + a.roll_step * 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
        model = field_normal(*rolled(positions, elements, contacts, best), pos, nrm, scale)
        amp_ratio = np.dot(model, measured) / np.dot(model, model)

        name = os.path.basename(fif)
        print(f'{name}: {duration:.0f} s, {info["n_chunks"]} chunks, f0 {f0:.4f} Hz '
              f'(chunk spread {np.ptp(info["chunk_f0"]) / f0 * 1e6:.0f} ppm, phase shifts <= {np.abs(info["chunk_shift_ms"]).max():.2f} ms), '
              f'{keep.sum()} magnetometers | best roll {best:.1f} deg, correlation {corr[i]:.3f}, measured/model amplitude {amp_ratio:.3f}')
        rows.append(dict(fif=name, set_angle_deg=set_angle(fif), duration_s=duration, f0_hz=f0, n_sensors=int(keep.sum()),
                         best_roll_deg=best, correlation=corr[i], amp_ratio=amp_ratio))
        np.savez(os.path.join(a.out, os.path.splitext(name)[0] + '.npz'), rolls=rolls, corr=corr, measured=measured,
                 model=amp_ratio * model, channels=np.array(raw.ch_names)[mags][keep], dev_head_t=raw.info['dev_head_t']['trans'])

    ref = rows[0]
    for r in rows:
        r['relative_deg'] = wrap(r['best_roll_deg'] - ref['best_roll_deg'])
        r['set_relative_deg'] = r['set_angle_deg'] - ref['set_angle_deg']
        e_plus, e_minus = wrap(r['relative_deg'] - r['set_relative_deg']), wrap(r['relative_deg'] + r['set_relative_deg'])
        r['error_deg'] = e_plus if abs(e_plus) <= abs(e_minus) else e_minus     # either rotation sense
    with open(os.path.join(a.out, 'rotation_results.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print('\nset angle -> rotation relative to the first recording:')
    for r in rows:
        print(f"  {r['fif']}: set {r['set_angle_deg']:.0f} deg, measured {r['relative_deg']:+.1f} deg, error {r['error_deg']:+.1f} deg, correlation {r['correlation']:.3f}")
    print('written', os.path.join(a.out, 'rotation_results.csv'))


if __name__ == '__main__':
    main()
