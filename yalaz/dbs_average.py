"""Average DBS pulse per MEG sensor, locked to the stimulation period.

The stimulator runs at f0 (about 130 Hz) on its own clock. Instead of detecting pulses, we take the
complex Fourier coefficients of every channel at the harmonics k*f0 (k = 1 .. lowpass/f0) and
synthesise the average period from them. This is identical to cutting the recording into one-period
segments and averaging, but works for a non-integer number of samples per period.

Because the stimulator clock drifts by a few ppm against the MEG clock, the recording is processed
in 10 s chunks: f0 is estimated per chunk, each chunk's coefficients are phase-aligned to the first
chunk, and the aligned chunks are averaged. The value used for each channel is the projection of its
average period onto the waveform shared by all channels (first principal component), which is more
robust than the largest sample of a possibly smeared waveform.
"""
import numpy as np
import mne

mne.set_log_level('ERROR')


def fourier_coefficients(X, t, f0, K):
    """Complex coefficients of X (n_channels, n_samples) at k*f0 for k = 1..K -> (n_channels, K)."""
    return np.stack([(X * np.exp(-2j * np.pi * k * f0 * t)).mean(axis=1) for k in range(1, K + 1)], axis=1)


def synthesise(C, f0, tau):
    """Average-period waveform (n_channels, len(tau)) from coefficients C."""
    k = np.arange(1, C.shape[1] + 1)
    return 2 * np.real(C @ np.exp(2j * np.pi * np.outer(k, f0 * tau)))


def harmonic_power(X, t, f0, K=12):
    return np.sum(np.abs(fourier_coefficients(X, t, f0, K)) ** 2)


def best_frequency(X, t, centre, halfwidth, step):
    grid = np.arange(centre - halfwidth, centre + halfwidth, step)
    return grid[np.argmax([harmonic_power(X, t, f) for f in grid])]


def pulse_average(raw, f_nominal=130.0, chunk_s=10.0):
    """Return (f0, tau, W, amp, info) for a raw recording.

    W    : average period per channel (n_channels, len(tau)), all channels in the fif
    amp  : per-channel amplitude = projection of W on the common waveform (the measured map)
    info : dict with the per-chunk frequencies and phase shifts (drift diagnostics)
    """
    raw = raw.copy().load_data().filter(60, None, method='iir', iir_params=dict(order=6, ftype='butter'))
    X, t, fs = raw.get_data(), raw.times, raw.info['sfreq']
    mags = mne.pick_types(raw.info, meg='mag')

    # the 10 magnetometers with the largest spectral peak at the nominal frequency drive the tracking
    F = np.fft.rfftfreq(X.shape[1], 1 / fs)
    S = np.abs(np.fft.rfft(X[mags], axis=1))
    k0 = np.argmin(np.abs(F - f_nominal))
    strong = mags[np.argsort(-S[:, k0 - 2:k0 + 3].max(axis=1))[:10]]

    # stimulation frequency over the whole recording: coarse then fine search
    f0 = best_frequency(X[strong], t, f_nominal, 1.0, 0.01)
    f0 = best_frequency(X[strong], t, f0, 0.02, 0.0005)
    K = int(raw.info['lowpass'] // f0)
    tau = np.arange(0, 1 / f0, 1 / fs)
    k = np.arange(1, K + 1)

    # chunks of chunk_s seconds (the last one shifted back to end at the recording's end)
    n, L = X.shape[1], int(chunk_s * fs)
    starts = list(range(0, n - L + 1, L)) or [0]
    if starts[-1] + L < n - L // 2:
        starts.append(n - L)

    coefs, weights, chunk_f0, chunk_shift, ref = [], [], [], [], None
    for s0 in starts:
        sl = slice(s0, min(s0 + L, n))
        fc = best_frequency(X[strong][:, sl], t[sl], f0, 0.1, 0.0005)          # +/- 770 ppm
        C = fourier_coefficients(X[:, sl], t[sl], fc, K)                          # absolute time keeps phases comparable
        if ref is None:
            ref, shift = C[strong], 0.0
        else:
            # cyclic time shift that best aligns this chunk with the first one on the strong channels
            shifts = np.linspace(-0.5 / f0, 0.5 / f0, 2001)
            phase = np.exp(-2j * np.pi * np.outer(shifts, k * f0))
            score = np.real(np.einsum('ck,dk,ck->d', C[strong], phase, np.conj(ref)))
            shift = shifts[int(np.argmax(score))]
            C = C * np.exp(-2j * np.pi * k * f0 * shift)[None, :]
        coefs.append(C); weights.append(sl.stop - sl.start); chunk_f0.append(fc); chunk_shift.append(shift)

    weights = np.array(weights, float) / np.sum(weights)
    W = synthesise(sum(w * C for w, C in zip(weights, coefs)), f0, tau)

    # common waveform across magnetometers; sign so that its main lobe is positive
    template = np.linalg.svd(W[mags], full_matrices=False)[2][0]
    if template[np.argmax(np.abs(template))] < 0:
        template = -template
    amp = W @ template

    info = dict(chunk_f0=np.array(chunk_f0), chunk_shift_ms=np.array(chunk_shift) * 1e3, n_chunks=len(coefs))
    return f0, tau, W, amp, info
