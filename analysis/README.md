# analysis

- `dbs_period_average.py` – period-locked average of the stimulation pulse from a raw fif (exact period from the
  harmonic comb, no peak detection); signed peak per sensor.
- `rotation_fit.py` – lead rotation per recording against an Ansys current export. Uses each fif's own dev_head_t,
  the Ansys axis convention (+x posterior, +y apex, +z left ear) or exact fiducials, a roll scan about the shaft,
  and reports noise-limited and relative rotations. This is the script to run on the full recordings.
- `fit_2a2b_roll_scan.py`, `ansys_forward.py` – earlier one-off comparisons, superseded.

Data notes
- `data/r0_ansys.fld_compressed.npz` ... `r90_...` are the valid exports (corrected geometry, true lead rotations; contacts at Ansys (3, -41, -63) mm = head (-3, -63, -41) mm). `J_2a_2b.npz` is the earlier geometry.
- `data/J_2a_3a.npz` has the lead in the wrong place (68 mm off) and must not be used for comparisons.
- `data/helmet_20s_raw.fif` (IPG on the bucket inside the helmet) is dominated by a 1.4 cm² current loop below the
  front rim; `data/chest_20s_raw.fif` (4 mA) and `data/chest_12s_new_raw.fif` (7.5 mA, 2a-2b) are clean.
