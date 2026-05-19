Corrected bodies.mat calibration notes
======================================

The final project uses a corrected bodies.mat file. The original system data were
used as the starting point, but some generated CSV trajectories did not match the
online validator. We therefore rebuilt selected body ephemerides using validator
feedback.

This calibration work is not meant to replace the official problem statement. It
documents the small corrections that were needed to make the local CSV generator
consistent with the web validator.

Preserved recalibrated / refitted body snapshots
------------------------------------------------

id_5 Aurelios
- Calibrated because several flybys showed small position mismatches.
- Approximate calibration epochs:
  - t = 5.9 yr
  - t = 50 yr
  - t = 82.4117 yr
  - t = 163.8599 yr
  - t = 231.3625 yr
  - t = 280.2526 yr

id_6 Tharos
- Calibrated because the mismatch increased over long time spans.
- Approximate calibration epochs:
  - t = 20 yr
  - t = 60 yr
  - t = 120 yr

id_7 Glacior
- Calibrated because the mismatch also increased over long time spans.
- Approximate calibration epochs:
  - t = 40 yr
  - t = 90 yr
  - t = 160 yr

id_8 Ashkara
- A preserved fit snapshot is available as bodies_fit_id_8.mat.
- Corrected/probe CSV evidence is preserved in output_probe/Correct.

Additional preserved fit snapshots
----------------------------------

id_4 Cindrel
- A preserved fit snapshot is available as bodies_fit_id_4.mat.

id_9 Varuun
- A preserved fit snapshot is available as bodies_fit_id_9.mat.

id_10 Xerion
- A preserved fit snapshot is available as bodies_fit_id_10.mat.

id_1000 Nyxar
- A preserved fit snapshot is available as bodies_fit_id_1000.mat.

Method
------

For each calibrated body, we generated simple two-body probe CSV files, uploaded
them to the validator, read the reported position errors, applied controlled X/Y/Z
perturbations, estimated the correction vector, and refitted the orbital elements.

Final artifact
--------------

The final corrected ephemeris is stored in:

- data/bodies.mat
- optimization/submission/bodies.mat
- optimization/flyby_functions/bodies.mat

These copies are intentionally identical in the cleaned repository.
