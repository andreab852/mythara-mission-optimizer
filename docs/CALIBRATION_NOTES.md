# Calibration Notes

The final project uses a corrected `bodies.mat` file.

During development, CSV trajectories generated from the original system data were
not consistent with the online validator. The working hypothesis was that some of
the initially provided ephemeris data did not match the validator's internal model.
To recover a validator-consistent model, the team rebuilt the affected bodies one at
a time.

## Calibration Workflow

The calibration work is preserved in:

```text
optimization/orbit_calibration/
optimization/README_bodies_calibration.txt
```

The workflow was:

1. Generate simple probe CSV files for one target body.
2. Upload those probe CSV files to the validator.
3. Read the validator position errors.
4. Apply controlled X/Y/Z perturbations.
5. Estimate a correction vector from the observed error response.
6. Fit corrected orbital elements from multiple validated correction points.
7. Install the corrected body into the final `bodies.mat`.

## Calibrated Bodies

The preserved fit snapshots are:

```text
optimization/orbit_calibration/bodies_fit/bodies_fit_id_4.mat
optimization/orbit_calibration/bodies_fit/bodies_fit_id_5.mat
optimization/orbit_calibration/bodies_fit/bodies_fit_id_6.mat
optimization/orbit_calibration/bodies_fit/bodies_fit_id_7.mat
optimization/orbit_calibration/bodies_fit/bodies_fit_id_8.mat
optimization/orbit_calibration/bodies_fit/bodies_fit_id_9.mat
optimization/orbit_calibration/bodies_fit/bodies_fit_id_10.mat
optimization/orbit_calibration/bodies_fit/bodies_fit_id_1000.mat
```

These correspond to Cindrel, Aurelios, Tharos, Glacior, Ashkara, Varuun, Xerion, and Nyxar. Each `.mat` file stores a full `bodies_fit` system struct after fitting the body identified in the filename.

The corresponding probe and corrected CSV files are preserved in:

```text
optimization/orbit_calibration/output_probe/Correct/
```

## Final Artifacts

The corrected final ephemeris is:

```text
data/bodies.mat
optimization/submission/bodies.mat
optimization/flyby_functions/bodies.mat
```

These three copies are intentionally identical in the cleaned repository. The copy inside `optimization/submission` is the one used by the final search and CSV exporter. The copy inside `data` is provided as the canonical data artifact for the cleaned repository.
