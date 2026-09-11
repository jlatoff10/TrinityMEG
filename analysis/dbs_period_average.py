"""Period-locked average of DBS pulses in a raw MEG recording (no peak detection needed).

Usage:  python analysis/dbs_period_average.py data/chest_20s_raw.fif out_chest.npz [f_nominal]

Method (after Yalaz et al.): 60 Hz high-pass, estimate the exact stimulation frequency f0 from the
harmonic comb, then synthesise the average period from the complex Fourier coefficients at k*f0
(k = 1 .. f_lowpass/f0). That is equivalent to cutting the record into one-period segments and
averaging, but exact for a non-integer number of samples per period. The per-sensor signed
peak of the averaged period is the quantity Yalaz compares with the FEM model.
"""
import sys, numpy as np, mne
mne.set_log_level('ERROR')

def period_average(raw, f_nominal=130.0, f_lp=None):
    raw = raw.copy().load_data().filter(60, None, method='iir', iir_params=dict(order=6, ftype='butter'))
    X, t, fs = raw.get_data(), raw.times, raw.info['sfreq']
    f_lp = f_lp or raw.info['lowpass']
    mags = mne.pick_types(raw.info, meg='mag')
    F = np.fft.rfftfreq(X.shape[1], 1 / fs); S = np.abs(np.fft.rfft(X[mags], axis=1))
    k = np.argmin(np.abs(F - f_nominal)); strong = mags[np.argsort(-S[:, k - 2:k + 3].max(1))[:10]]

    def hpow(f0):
        return sum(np.sum(np.abs((X[strong] * np.exp(-2j * np.pi * h * f0 * t)).mean(1)) ** 2) for h in range(1, 13))
    g = np.arange(f_nominal - 1, f_nominal + 1, 0.01); f0 = g[np.argmax([hpow(f) for f in g])]
    g = np.arange(f0 - 0.02, f0 + 0.02, 0.0005); f0 = g[np.argmax([hpow(f) for f in g])]

    K = int(f_lp // f0); tau = np.arange(0, 1 / f0, 1 / fs); W = np.zeros((X.shape[0], len(tau)))
    for h in range(1, K + 1):
        c = (X * np.exp(-2j * np.pi * h * f0 * t)).mean(1)
        W += 2 * np.real(c[:, None] * np.exp(2j * np.pi * h * f0 * tau)[None, :])
    pk = W[np.arange(len(W)), np.argmax(np.abs(W), 1)]          # signed peak per channel
    return f0, tau, W, pk

if __name__ == '__main__':
    raw = mne.io.read_raw_fif(sys.argv[1], preload=False)
    f0, tau, W, pk = period_average(raw, float(sys.argv[3]) if len(sys.argv) > 3 else 130.0)
    T = raw.info['dev_head_t']['trans']
    loc = np.array([ch['loc'] for ch in raw.info['chs']])
    np.savez(sys.argv[2], f0=f0, tau=tau, W=W, pk=pk, names=np.array(raw.ch_names),
             mags=mne.pick_types(raw.info, meg='mag'), grads=mne.pick_types(raw.info, meg='grad'),
             pos_dev=loc[:, :3], nrm_dev=loc[:, 9:12], dev_head_t=T)
    m = mne.pick_types(raw.info, meg='mag')
    print(f'f0 = {f0:.4f} Hz; magnetometer signed peaks: min {pk[m].min():.2e} max {pk[m].max():.2e} T, '
          f'{(pk[m] > 0).sum()} positive / {(pk[m] < 0).sum()} negative')
