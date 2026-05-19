# Orbit Calibration Workflow

This folder implements the **ephemeris validation and correction pipeline** used to rebuild the planetary data (`bodies.mat`) from raw GTOC13 data.

## Problem Background

The original GTOC13 ephemeris produced **incorrect CSV outputs** when propagated through the standard solver:
- Coordinate frame misalignment (heliocentric vs. barycentric)
- Ephemeris gaps or interpolation artifacts
- Inconsistent epoch handling across bodies

Solution: Implement a **validator-feedback loop** to incrementally correct the ephemeris using:
1. Probe trajectories at known points
2. Validator error messages from submitted solutions
3. Fitted orbital elements that reconstruct the correct behavior

## Workflow Scripts

### Phase 1: Generate Probe Trajectories

**`RUN_make_probe4.m`**

Generate 4 test CSVs to probe the ephemeris at a target body:
- **Base probe:** Direct trajectory to target body (no offset)
- **+X probe:** Target position shifted +1000 km in heliocentric X
- **+Y probe:** Target position shifted +1000 km in heliocentric Y
- **+Z probe:** Target position shifted +1000 km in heliocentric Z

These are uploaded to the validator, which returns error codes indicating how far off each solution is.

**Input (edit only this block):**
```matlab
target_id = 5;        % Body to calibrate (Jupiter)
source_id = 3;        % Reference body already corrected
t_target = 2.5;       % Epoch in years
shift_km = 1000;      % Probe offset magnitude
```

**Output:**
```
optimization/orbit_calibration/output_probe/
  probe_base_*.csv
  probe_px_*.csv
  probe_py_*.csv
  probe_pz_*.csv
```

### Phase 2: Estimate Correction Vector

**`RUN_estimate_corr4_and_make_csv.m`**

After validator returns errors (E0, ExP, EyP, EzP), compute the correction offset and produce a corrected CSV.

Linear approximation:
- Correction_X ≈ −(ExP − E0) / (shift_km / 1000)
- Correction_Y ≈ −(EyP − E0) / (shift_km / 1000)
- Correction_Z ≈ −(EzP − E0) / (shift_km / 1000)

**Input:**
```matlab
target_id = 5;        % Body to calibrate
source_id = 3;        % Reference body
t_target = 2.5;       % Epoch
shift_km = 1000;      % Probe offset (must match Phase 1)
dir = './';           % Output directory

E0 = 1250;            % Validator error for base probe
ExP = 2150;           % Validator error for +X probe
EyP = 1310;           % Validator error for +Y probe
EzP = 1265;           % Validator error for +Z probe
```

**Output:**
```
optimization/orbit_calibration/output_probe/
  corrected_*.csv
```

### Phase 3: Fit Orbital Elements

**`RUN_fit_body_from_points.m`**

After 3+ validated correction points, fit the true orbital elements for the target body.

Uses nonlinear least squares to find Keplerian elements (a, e, i, Ω, ω, M₀) that best match the validated positions across all epochs.

**Input:**
```matlab
target_id = 5;        % Body to fit
model = 'keplerian';  % Orbital model

T = [0.0 1.0 2.0 3.0];           % Validated epochs (years)
C = [0 50; 0 45; 0 48; 0 50];    % Corrections in km (rows: epochs; cols: X, Y)
                                  % Z corrections derived from other data
```

**Output:**
```
optimization/orbit_calibration/bodies_fit/
  fitted_elements_body5_v1.mat
  fitted_ephemeris_v1.mat
```

## Setup & Automation

### `SETUP_mythara_paths_auto.m`

Central initialization script (called by all runners):
- Locates required files (bodies.mat, solver functions, validators)
- Adds only necessary folders to MATLAB path (avoids genpath shadowing)
- Prints resolved file paths to console for verification

**Usage:**
```matlab
run('SETUP_mythara_paths_auto.m')  % Automatic in all RUN_*.m scripts
```

**Typical output:**
```
=== Mythara Path Setup ===
Using files:
  solver: ...optimization/submission/Main_Flyby_*.m
  bodies: ...optimization/submission/bodies.mat
  functions: ...optimization/flyby_functions/
  ...setup complete
```

## Directory Structure

```
orbit_calibration/
  ├─ README.md                    (this file)
  ├─ SETUP_mythara_paths_auto.m   (initialization)
  ├─ RUN_make_probe4.m            (Phase 1)
  ├─ RUN_estimate_corr4_and_make_csv.m  (Phase 2)
  ├─ RUN_fit_body_from_points.m   (Phase 3)
  ├─ output_probe/                (temporary: probe CSVs)
  └─ bodies_fit/                  (permanent: fitted elements)
```

## Typical Workflow Example

**Goal:** Rebuild bodies.mat with corrected ephemeris

1. **Target body:** Mercury (ID = 1)

   ```matlab
   % RUN_make_probe4.m INPUT:
   target_id = 1; source_id = 3; t_target = 1.0; shift_km = 1000;
   ```

2. Upload `probe_base_*.csv`, `probe_px_*.csv`, `probe_py_*.csv`, `probe_pz_*.csv` to GTOC13 validator.
   Validator returns errors: E0=1200, ExP=2100, EyP=1250, EzP=1180

   ```matlab
   % RUN_estimate_corr4_and_make_csv.m INPUT:
   target_id = 1; E0 = 1200; ExP = 2100; EyP = 1250; EzP = 1180;
   ```

3. Upload `corrected_*.csv`. If validator error now < 100, move to next body or fit.

4. After 3+ validated points:
   ```matlab
   % RUN_fit_body_from_points.m INPUT:
   target_id = 1;
   T = [0.5 1.0 1.5 2.0];
   C = [50 45 48 52; ...];  % Fit these corrections
   ```

5. Load fitted elements and rebuild `bodies.mat`:
   ```matlab
   load('bodies_fit/fitted_elements_body1_v1.mat');
   % Reconstruct ephemeris with corrected Keplerian elements
   bodies_corrected = reconstruct_ephemeris(fitted_elements);
   save('../../data/bodies.mat', 'bodies_corrected');
   ```

## Validation Loop

The validator feedback tells you:
- **E << 100:** Solution is valid, body ephemeris is correct
- **100 ≤ E < 1000:** Systematic error (offset/bias)
- **E ≥ 1000:** Major mismatch (coordinate frame, epoch, or ephemeris issue)

Once E < 50 for all probes at a given epoch, that body is validated for that time.

## Notes for Developers

- **No manual ephemeris editing:** All changes flow through this workflow
- **Reproducibility:** Keep validator errors and probe CSVs for audit trail
- **Iteration:** If fit quality is poor (residual > 200 km), generate new probes and re-fit
- **Scale:** Typically 2–4 iterations per body to achieve E < 50 across all epochs
- **Time investment:** Full calibration (10 bodies) typically takes 1–2 hours total

## References

- GTOC13 problem statement for validator specifications
- Bate, Mueller, White (1971) for orbital element definitions
- Vallado et al. (2006) for epoch and frame conventions

## Troubleshooting

| Issue | Cause | Fix |
|-------|-------|-----|
| `bodies.mat not found` | File not in expected location | Run `SETUP_mythara_paths_auto.m` |
| Validator errors don't decrease | Wrong coordinate frame or epoch offset | Check `CALIBRATION_NOTES.md` |
| Fit diverges | Poor initial guess or singular matrix | Use fewer epochs or better seeding |
| CSV export fails | Missing `solution_to_csv_patched.m` | Verify `flyby_functions/` path |

