"""Period-locked average of DBS pulses in a raw MEG recording (no peak detection needed).

Usage:  python analysis/dbs_period_average.py data/chest_20s_raw.fif out_chest.npz [f_nominal]

Method (after Yalaz et al.): 60 Hz high-pass, then the average pulse is synthesised from the complex
Fourier coefficients at the harmonics k*f0 (k = 1 .. f_lowpass/f0). That equals cutting the record
into one-period segments and averaging, but exact for a non-integer number of samples per period.

The IPG clock is not locked to the MEG clock and drifts by parts per million, which smears a single
f0 average over minutes. So the record is processed in chunks: f0 is re-estimated per chunk, each
chunk's average period is phase-aligned to the first chunk on the strongest channels, and the aligned
chunk averages are combined. The per-channel amplitude is the projection onto the common waveform
(first principal component), which is robust where a per-channel signed peak can jump between lobes.
"""
import sys, numpy as np, mne
mne.set_log_level('ERROR')


def _coefs(X, t, f0, K):
    """Complex Fourier coefficients at k*f0, k=1..K: shape (n_ch, K)."""
    return np.stack([(X * np.exp(-2j * np.pi * k * f0 * t)).mean(1) for k in range(1, K + 1)], 1)


def _synth(C, f0, tau):
    """Waveform (n_ch, len(tau)) from coefficients."""
    k = np.arange(1, C.shape[1] + 1)
    return 2 * np.real(C @ np.exp(2j * np.pi * np.outer(k, f0 * tau)))


def _hpow(X, t, f0, K):
    return np.sum(np.abs(_coefs(X, t, f0, K)) ** 2)


def period_average(raw, f_nominal=130.0, f_lp=None, chunk_s=10.0, return_extras=False):
    raw = raw.copy().load_data().filter(60, None, method='iir', iir_params=dict(order=6, ftype='butter'))
    X, t, fs = raw.get_data(), raw.times, raw.info['sfreq']
    f_lp = f_lp or raw.info['lowpass']
    mags = mne.pick_types(raw.info, meg='mag')
    # strongest magnetometers at the nominal frequency, used for the frequency/phase tracking
    F = np.fft.rfftfreq(X.shape[1], 1 / fs); S = np.abs(np.fft.rfft(X[mags], axis=1))
    k0 = np.argmin(np.abs(F - f_nominal)); strong = mags[np.argsort(-S[:, k0 - 2:k0 + 3].max(1))[:10]]
    K = int(f_lp // f_nominal)
    # global f0: coarse then fine search of the harmonic power on the whole record
    g = np.arange(f_nominal - 1, f_nominal + 1, 0.01); f0 = g[np.argmax([_hpow(X[strong], t, f, 12) for f in g])]
    g = np.arange(f0 - 0.02, f0 + 0.02, 0.0005); f0 = g[np.argmax([_hpow(X[strong], t, f, 12) for f in g])]
    K = int(f_lp // f0); tau = np.arange(0, 1 / f0, 1 / fs)
    # chunks: per-chunk f0, per-chunk coefficients, phase-aligned to the first chunk
    n = X.shape[1]; L = int(chunk_s * fs); starts = list(range(0, n - L + 1, L)) or [0]
    if starts[-1] + L < n - L // 2: starts.append(n - L)
    Cs, f0s, shifts, wts = [], [], [], []
    ref = None; k = np.arange(1, K + 1)
    for s0 in starts:
        sl = slice(s0, min(s0 + L, n)); ts = t[sl]                      # absolute time: chunk phases comparable without drift
        g = np.arange(f0 - 0.1, f0 + 0.1, 0.0005)                         # +/- 770 ppm search per chunk
        fc = g[np.argmax([_hpow(X[strong][:, sl], ts, f, 12) for f in g])]
        C = _coefs(X[:, sl], ts, fc, K)
        if ref is None:
            ref = C[strong]; d = 0.0
        else:
            # cyclic time shift d that maximises the match of this chunk to the reference on the strong channels
            dd = np.linspace(-0.5 / f0, 0.5 / f0, 2001)
            ph = np.exp(-2j * np.pi * np.outer(dd, k * f0))                       # (n_d, K)
            score = np.real(np.einsum('ck,dk,ck->d', C[strong], ph, np.conj(ref)))
            d = dd[int(np.argmax(score))]
            C = C * np.exp(-2j * np.pi * k * f0 * d)[None, :]
        Cs.append(C); f0s.append(fc); shifts.append(d); wts.append(sl.stop - sl.start)
    wts = np.array(wts, float) / np.sum(wts)
    Cmean = np.sum([w * C for w, C in zip(wts, Cs)], 0)
    W = _synth(Cmean, f0, tau)
    pk = W[np.arange(len(W)), np.argmax(np.abs(W), 1)]                     # signed peak per channel (Yalaz)
    # common waveform: first principal component over the magnetometers; amplitude = projection onto it
    u, s, vt = np.linalg.svd(W[mags], full_matrices=False); templ = vt[0]
    if templ[np.argmax(np.abs(templ))] < 0: templ = -templ                 # sign: the template's main lobe is positive, so amp has the sign of the signed peak
    amp = W @ templ
    # per-channel noise of the amplitude from the scatter between chunk averages (standard error of the mean)
    amps_c = np.array([_synth(C, f0, tau) @ templ for C in Cs])
    noise = amps_c.std(0, ddof=1) / np.sqrt(len(Cs)) if len(Cs) > 1 else np.full(len(amp), np.nan)
    extras = dict(amp=amp, noise=noise, template=templ, chunk_f0=np.array(f0s), chunk_shift_ms=np.array(shifts) * 1e3,
                  chunk_starts_s=np.array([t[s] for s in starts]), n_chunks=len(Cs), template_var_fraction=s[0] ** 2 / np.sum(s ** 2))
    return (f0, tau, W, pk, extras) if return_extras else (f0, tau, W, pk)


if __name__ == '__main__':
    raw = mne.io.read_raw_fif(sys.argv[1], preload=False)
    f0, tau, W, pk, ex = period_average(raw, float(sys.argv[3]) if len(sys.argv) > 3 else 130.0, return_extras=True)
    T = raw.info['dev_head_t']['trans']
    loc = np.array([ch['loc'] for ch in raw.info['chs']])
    np.savez(sys.argv[2], f0=f0, tau=tau, W=W, pk=pk, amp=ex['amp'], noise=ex['noise'], names=np.array(raw.ch_names),
             mags=mne.pick_types(raw.info, meg='mag'), grads=mne.pick_types(raw.info, meg='grad'),
             pos_dev=loc[:, :3], nrm_dev=loc[:, 9:12], dev_head_t=T, chunk_f0=ex['chunk_f0'], chunk_shift_ms=ex['chunk_shift_ms'])
    m = mne.pick_types(raw.info, meg='mag')
    print(f'f0 = {f0:.4f} Hz over {ex["n_chunks"]} chunks; per-chunk f0 spread {np.ptp(ex["chunk_f0"]) * 1e3:.1f} mHz, '
          f'phase shifts {np.round(ex["chunk_shift_ms"], 3)} ms')
    print(f'template explains {ex["template_var_fraction"] * 100:.0f}% of the magnetometer waveform variance; '
          f'signed peaks {pk[m].min():.2e}..{pk[m].max():.2e} T, projected amplitudes {ex["amp"][m].min():.2e}..{ex["amp"][m].max():.2e} T, '
          f'median noise {np.nanmedian(ex["noise"][m]):.1e} T')
