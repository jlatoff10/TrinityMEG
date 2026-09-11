# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 12:28:38 2026

@author: trini
"""

import mne
raw = mne.io.read_raw_fif(
    'C:/Users/nebula/Desktop/260904_r0_a2a_c2b_4_60_130_raw.fif', 
    preload=True)

# # Crop the data from tmin (start) to tmax (end) in seconds
# raw.crop(tmin=0.0, tmax=100.0)

# # Save the shortened file
# raw.save('C:/Users/trini/Documents/noah_noah/260818/260818/260818_a2a_c2b_1.fif', overwrite=True)



raw.compute_psd(fmax=150).plot(picks="data", amplitude=False)





import numpy as np
from mne.preprocessing import peak_finder
import matplotlib.pyplot as plt


preload=True
iir_params = dict(order=6, ftype="butter")
raw_notched = raw.copy().notch_filter(
    freqs=[60, 120, 180, 240], 
    method='fir', 
    notch_widths=1.0) 
    #iir_params=iir_params)
raw_notched.compute_psd(fmax=150).plot()
#plt.show()

iir_params = dict(order=6, ftype='butter', )
raw_filtered = raw_notched.filter(
    l_freq=60,
    h_freq=None,
    method='iir',
    iir_params=iir_params)


data = raw_filtered.get_data()
times = raw_filtered.times

results = []
for ch_idx, ch_name in enumerate(raw_filtered.ch_names):
    channel_data = data[ch_idx, :]
    peak_locs, peak_mags = peak_finder(channel_data)

    if len(peak_locs) == 0:
        continue
    peak_times = times[peak_locs]

    for loc, t, mag in zip(peak_locs, peak_times, peak_mags):
        results.append({
            "channel": ch_name,
            "time_sec": t,
            "amplitude": mag})



from scipy.signal import butter, filtfilt
b, a = butter(6, 60/(5000/2), 'highpass')
xf = filtfilt(b, a, raw_filtered.get_data(picks='grad'))








sfreq = raw_filtered.info['sfreq']
peak_samples = raw.time_as_index(peak_times) + raw_filtered.first_samp

peaks = {}
for i in range(xf.shape[0]):
    ch = xf[i]
    sign = 1 if abs(ch.max()) > abs(ch.min()) else -1
    thresh = 4 * np.median(np.abs(ch - np.median(ch)))
    locs, mags = peak_finder(ch, thresh=thresh, extrema=sign)
    keep = np.abs(mags) > 4e-10
    print('rejected =', (~keep).sum())
    print('rejected = ',len(mags) - len(keep))
    locs, mags = locs[keep], mags[keep]
    peaks[i] = (locs, mags)



plt.show()


#----------------------------------------------------------------------------

print(len(raw))


grad_names = raw_filtered.copy().pick('grad').ch_names
n_peaks = np.array([len(peaks[i][0]) for i in range(xf.shape[0])])

for name, n in zip(grad_names, n_peaks):
    print(f"{name}: {n} peaks")

for i in range(len(peaks)):
    locs, mags = peaks[i]
    if len(mags) == 23554:
        print(i)
        
#meg2543: 23554, 195

#----------------------------------------------------------------------------    



locs, mags = peaks[195]
sfreq = raw_filtered.info['sfreq']
peak_samples = locs + raw_filtered.first_samp


events = np.column_stack((
    peak_samples,
    np.zeros(len(peak_samples), dtype=int),
    np.ones(len(peak_samples), dtype=int)))

epoched = mne.Epochs(
    raw_filtered,
    events=events,
    event_id=1,
    tmin=-0.003,
    tmax=0.003,
    baseline=(-0.002, 0),
    proj=False,
    reject=None,
    preload=True)

evoked=epoched.average(method='mean', by_event_type=False)






ch_name = 'MEG2543'
grad_names = raw_filtered.copy().pick('grad').ch_names
i = grad_names.index(ch_name)

ch = xf[i]
locs, mags = peaks[i]

t = np.arange(len(ch)) / sfreq

fig, ax = plt.subplots(1, 1, figsize=(13, 9))

w = (t >= 10.0) & (t < 10.05)
ax.plot(t[w] * 1e3, ch[w] * 1e12, lw=0.8, label=ch_name)

sel = locs[(t[locs] >= 10.0) & (t[locs] < 10.05)]
ax.plot(t[sel] * 1e3, ch[sel] * 1e12, 'rx', ms=9, label='peaks')

ax.set_xlabel('Time (ms)')
ax.set_ylabel('Amplitude (pT/cm)')  # adjust units label if needed for grad
ax.legend()
plt.show()





epoched.plot(picks='grad', n_epochs=1, butterfly=True)
evoked.plot(spatial_colors=True)
evoked.plot_topo()
evoked.save('act_r0_2a_2b.fif')





data = evoked.picks('mags').data  # shape: (n_channels, n_times)
times = evoked.times

# Find peak (max absolute amplitude) per sensor
peak_indices = np.argmax(np.abs(data), axis=1)
peak_amplitudes = data[np.arange(data.shape[0]), peak_indices]
peak_latencies = times[peak_indices]

for ch_name, amp, lat in zip(evoked.ch_names, peak_amplitudes, peak_latencies):
    print(f"channel_name: {ch_name},{amp:.5e}")
    print(len(peak_amplitudes))




