"""Model-free rotation calibration from the measured maps of several recordings at known set angles.

For a segment pair rotating about the shaft, the field at the sensors is linear in the dipole vector, so the
measured map is m(theta) = C + A cos(theta) + B sin(theta) to first order. With recordings at known angles, C, A, B
are fitted by least squares from the data alone; each recording's angle is then estimated leave-one-out
(calibrating on the others) and compared with its set angle. No Ansys model involved.

Usage: python analysis/rotation_empirical.py analysis/results/260911_r*_r0_ansys.fld_compressed.npz
       (the per-recording .npz files written by rotation_fit.py; the set angle is parsed from the file name)
"""
import sys, re, os, glob, numpy as np

def load(paths):
    recs = []
    for p in paths:
        m = re.search(r'_r(\d+)_', os.path.basename(p))
        if not m: continue
        d = np.load(p); recs.append((float(m.group(1)), d['measured'], d['sigma'] if 'sigma' in d else None, os.path.basename(p)))
    return sorted(recs, key=lambda r: r[0])

def fit_cab(angles, maps, sense=1.0):
    th = np.radians(angles) * sense
    X = np.column_stack([np.ones_like(th), np.cos(th), np.sin(th)])       # (n_rec, 3)
    coef, *_ = np.linalg.lstsq(X, np.array(maps), rcond=None)             # (3, n_ch): C, A, B maps
    return coef

def estimate(coef, m):
    C, A, B = coef
    X = np.column_stack([C, A, B, np.ones_like(C)])
    c, *_ = np.linalg.lstsq(X, m, rcond=None)                              # free scale on C absorbs amplitude drift
    return np.degrees(np.arctan2(c[2], c[1])) % 360, c

def main(paths, exclude=None):
    recs = load(paths)
    if exclude is not None:
        d0 = np.load([p for p in paths if p.endswith('.npz')][0])
        if 'pos_dev' in d0:
            keep = np.linalg.norm(d0['pos_dev'] - np.array(exclude[:3]), axis=1) > exclude[3]
            recs = [(a, m[keep], (s[keep] if s is not None else None), n) for a, m, s, n in recs]; print(f'using {keep.sum()} sensors farther than {exclude[3] * 100:.0f} cm from the excluded point')
    if len(recs) < 4: print('need at least 4 recordings at different set angles'); return
    angles = np.array([r[0] for r in recs]); maps = [r[1] for r in recs]
    print('recordings:', ', '.join(f'r{int(a)}' for a in angles))
    print('pairwise correlation of the measured maps (should be ~1 for small angle differences):')
    for i, r in enumerate(recs): print(f'  r{int(r[0]):<3d} ' + ' '.join(f'{np.corrcoef(maps[i], maps[j])[0, 1]:5.2f}' for j in range(len(recs))))
    for sense in (1.0, -1.0):
        errs = []
        print(f'\nleave-one-out, rotation sense {"+" if sense > 0 else "-"}:')
        for i, r in enumerate(recs):
            idx = [j for j in range(len(recs)) if j != i]
            coef = fit_cab(angles[idx], [maps[j] for j in idx], sense)
            est, c = estimate(coef, maps[i]); est = (est / sense) % 360
            err = (est - angles[i] + 180) % 360 - 180; errs.append(err)
            print(f'  r{int(angles[i]):<3d}: estimated {est:6.1f} deg, error {err:+6.1f} deg, scale on invariant part {c[0]:.2f}, rotating amplitude {np.hypot(c[1], c[2]):.2f}')
        print(f'  rms error {np.sqrt(np.mean(np.square(errs))):.1f} deg, max {np.max(np.abs(errs)):.1f} deg')

if __name__ == '__main__':
    args = sys.argv[1:]; exclude = None
    if '--exclude-near' in args:
        i = args.index('--exclude-near'); exclude = [float(v) for v in args[i + 1:i + 5]]; args = args[:i] + args[i + 5:]
    main(sum([glob.glob(p) for p in args], []), exclude)
