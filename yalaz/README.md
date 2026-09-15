# Directional DBS lead rotation from MEG, Yalaz's method

Self-contained code and results for determining the rotation of a directional DBS lead in the
bucket phantom from MEG recordings of the stimulation pulses, following Yalaz et al. (measured
magnetometer pattern correlated with a simulated pattern at every rotation of the lead). Three short
modules, no dependency on the rest of the repository:

| File | What it does |
|---|---|
| `dbs_average.py` | average stimulation pulse per sensor, locked to the stimulation period, drift-tolerant |
| `ansys_model.py` | field of the Ansys current export at the magnetometers, rotation of the lead about its shaft |
| `yalaz_rotation.py` | the analysis: roll scan per recording, rotation relative to the first recording, error vs set angle |
| `fld_to_npz.py` | converts an Ansys `.fld` grid export to the `.npz` the model reads |

Requirements: Python 3, numpy, scipy, pandas, mne.

## Setup that was analysed

Bucket phantom filled with 0.9% saline, Boston Scientific Vercise PC and Cartesia directional lead,
bipolar stimulation between adjacent segments of one level (2a anode, 2b cathode), 7.5 mA, 60 µs,
130 Hz. The IPG was at the chest, outside the helmet. Recordings of 2 to 3 minutes each at set
rotations of the lead of 0, 2, 5, 11, 22, 45 and 90 degrees (session 260911, files
`260911_r<angle>_2a_2b_7p5_60_130_raw.fif`), MEGIN TRIUX, 5 kHz, 0.1 to 1650 Hz, 102 magnetometers used.

The Ansys Maxwell model is a DC conduction solve of the saline with 1 V between the two contacts,
exported as current density on a 1 mm grid covering the whole saline volume (`data/r0_ansys.fld_compressed.npz`).
The current in that solve is 2.17 mA (Maxwell surface integral of J·n over the 1 V face). The model's
coordinate origin is the LPA–RPA midpoint and its x–z plane is the fiducial plane, so it sits in the MEG
head frame by a fixed axis permutation (Ansys +x left ear, +y apex, +z nose); each recording's HPI fit
(`dev_head_t`) moves the sensors into that frame. The contacts are at head (−3, −63, −41) mm, which is
41 mm below the fiducial plane, about 2 cm below the helmet rim, 65 mm from the nearest magnetometer.

## Method

1. **Pulse average per sensor.** The recording is high-passed at 60 Hz. The stimulation frequency f0 is
   found from the harmonic comb (130.1425 Hz). The average period is synthesised from the complex
   Fourier coefficients at k·f0, k = 1..12 (up to the 1650 Hz low-pass), which equals segmenting the
   recording into single periods and averaging. Because the IPG clock drifts a few ppm against the MEG
   clock, this is done in 10 s chunks, each with its own f0 and phase-aligned to the first chunk, then
   averaged. The value per sensor is the projection of its average period onto the waveform common to
   all sensors (first principal component). The result is one signed number per magnetometer: the
   measured map. Noise after averaging is about 20 fT per sensor for a 10 s chunk against signals of 5 to 13 pT.
2. **Sensor selection.** Magnetometers within 12 cm of device point (0.095, −0.04, 0.085) m are dropped
   (72 of 102 remain). See "What was found" for why.
3. **Model pattern.** The current elements of the export (binned to 3 mm) are rotated about the lead's
   shaft in 1° steps; at each roll the normal field component at the selected magnetometers is computed
   by Biot–Savart, scaled by 7.5 mA / 2.17 mA.
4. **Fit.** Pearson correlation between model pattern and measured map at every roll; the maximum,
   refined by a parabola through its three neighbouring points, is the lead orientation. Rotations
   are reported relative to the first recording (the model's zero is arbitrary), and compared with
   the set angle for either rotation sense.

Run:

```
python yalaz/yalaz_rotation.py --ansys data/r0_ansys.fld_compressed.npz --exclude-near 0.095 -0.04 0.085 0.12 \
    --out yalaz/results --fif 260911/260911_r0_2a_2b_7p5_60_130_raw.fif 260911/260911_r2_2a_2b_7p5_60_130_raw.fif \
    260911/260911_r5_2a_2b_7p5_60_130_raw.fif 260911/260911_r11_7p5_60_130_raw.fif 260911/260911_r22_2a_2b_7p5_60_130_raw.fif \
    260911/260911_r45_2a_2b_7p5_60_130_raw.fif 260911/260911_r90_2a_2b_7p5_60_130_raw.fif
```

About two minutes per recording. Output: `rotation_results.csv` and one `.npz` per recording with the
measured map, the model map at the best roll and the roll curve.

## Results (session 260911, 72 magnetometers)

| Set angle | Best roll | Correlation | Rotation relative to r0 | Error |
|---|---|---|---|---|
| 0° | 202.2° | 0.81 | 0 | 0 |
| 2° | 200.5° | 0.74 | −1.7° | +0.3° |
| 5° | 193.2° | 0.74 | −9.0° | −4.0° |
| 11° | 185.0° | 0.75 | −17.2° | −6.2° |
| 22° | 198.0° | 0.56 | −4.2° | +17.8° |
| 45° | 34.6° | 0.44 | −167.6° | −122.6° |
| 90° | 255.0° | 0.94 | +52.8° | −37.2° |

- The 2°, 5° and 11° steps are recovered within 6° with a consistent rotation sense (4° rms).
  This matches the accuracy Yalaz reports for a single horizontal-pair recording in his cylinder phantom.
- Amplitude: the measured map is 0.06 to 0.07 of the model. About 0.19 is expected from the 1650 Hz
  low-pass acting on a 60 µs pulse, so the recorded field is roughly 3 times weaker than the model at
  7.5 mA. Possible reasons are compliance-limited current at the high impedance of a segment pair,
  or a saline volume in the phantom that differs from the model. Rotation is unaffected by amplitude.
- r22 fits the model worse than its neighbours (0.56) and its rotation is not recovered.
- r45 and r90 are not rotations of the r0 pattern by 45° and 90°: r45 fits nothing, r90 is the cleanest
  recording of the set (0.94) but sits 53° from r0 in the opposite sense to the small-angle series.
  The lead was not at those set angles; the mechanism or its zero moved.

## What was found along the way

- **Signed-peak maps fail on long recordings.** A single-frequency average over 2 minutes smears
  because of the IPG clock drift, and picking each sensor's largest sample of a smeared waveform
  scrambles the map. The chunked, phase-aligned average with template projection fixes this.
- **A moving artifact inside the helmet.** Differences between recordings that should be near-identical
  (r0, r2, r5) localise to a compact source near the upper-right sensors, device (90 to 100, −35 to
  −47, 80 to 88) mm, with the moment of a millimetre-sized residual loop of the stimulation current:
  the lead connector or excess lead, handled at each angle change and relaxing during recordings.
  It contributes 300 fT of scatter within a recording and tens of pT between recordings. Its field falls
  as 1/r³, so magnetometers more than 12 cm from it are clean; hence the sensor selection. With all 102
  sensors the same method gives errors above 100°.
- **Noise is irrelevant.** The noise-limited rotation uncertainty is 0.2° per recording; all error is
  systematic (artifact, model, mechanism).
- **The earlier Ansys export of the same configuration had the lead in the wrong place** (68 mm off) and
  was discarded; the r-series exports are consistent with the digitised bucket (same taper, tilt and
  axis within 5 mm).

## Recommendations for the next session

Route the lead's proximal end, connector and any excess lead straight out of the helmet and fix them,
at least 15 cm below the rim; do not touch them between angle settings. Verify with two recordings at
the same angle before the series (maps should correlate above 0.98). Use angles spread over the full
range (0, 30, 60, 90, 120, 150, and 0 again at the end) so the rotation sense is determined and drift
is measured. Getting the contacts inside the sensor array, as in Yalaz's phantom, would raise the signal
by an order of magnitude.
