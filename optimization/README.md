# Optimization & Solver

This folder contains the complete trajectory optimization pipeline: the final multi-strategy parallel search solver, core flyby mechanics functions, and calibration utilities.

## Folder Structure

### `submission/`

**Entry point for running the final solver.**

- `Main_Flyby_multistrategy.m` — **Primary script** for the parallel multi-body search with 16 strategies and adaptive tuning. Generates all candidate solutions.
- `Main_Flyby_multistrategy_LIGHT.m` — Lightweight variant for quick prototyping.
- `bodies.mat` — Canonical planetary ephemeris (copy of `../data/bodies.mat`).
- `utils/` — Utility functions supporting the main solver.
- `Correct/` — Corrected trajectory validation and CSV export routines.

**How to run:**
```matlab
cd('optimization/submission');
Main_Flyby_multistrategy
```

Or from the project root:
```matlab
run_final_search
```

### `flyby_functions/`

Core orbital mechanics library:
- **`lambert_universal.m`** — Universal Lambert solver for Keplerian transfers
- **`KM_solver.m`** — Kerbal Motion propagator for high-precision orbit integration
- **`flyby_opt.m`** — Multi-impulse flyby optimization (Δv minimization)
- **`solution_to_csv_patched.m`** — Trajectory-to-CSV converter with epoch handling
- Supporting files for orbital element extraction, frame transformations, and constraint checking

These are used by the main solver but can be called independently for analysis or custom trajectory design.

### `orbit_calibration/`

Ephemeris validation and calibration workflow:

- `README.md` — Detailed guide to the calibration process
- `RUN_make_probe4.m` — Generate 4-point probe for target body correction
- `RUN_estimate_corr4_and_make_csv.m` — Compute correction vector from validator feedback
- `RUN_fit_body_from_points.m` — Fit orbital elements using validated correction points
- `SETUP_mythara_paths_auto.m` — Automatic path setup (no manual genpath needed)
- `output_probe/` — Temporary output for probe CSVs during calibration
- `bodies_fit/` — Orbit fit snapshots and corrected ephemeris versions

**When to use:**
- If you need to rebuild or refine the planetary ephemeris
- If validation suggests coordinate system issues
- Normally not required for running the final solver

### Root Files

- `SETUP_mythara_paths_auto.m` — Initializes paths for the optimization module (called by all runners)
- `RUN_*.m` — Calibration workflow scripts (see `orbit_calibration/README.md`)
- `README_bodies_calibration.txt` — Legacy notes on the ephemeris fix

## Quick Start

To generate new candidate trajectories:

1. **From MATLAB project root:**
   ```matlab
   run_final_search
   ```

2. **Or directly in optimization/submission:**
   ```matlab
   cd('optimization/submission');
   Main_Flyby_multistrategy;
   ```

The script will:
- Load bodies.mat and validate ephemeris
- Execute 16 parallel search strategies
- Output candidate trajectories to CSV
- Display performance statistics

## Algorithm Overview

**Multi-Strategy Parallel Search:**
- 16 concurrent search strategies (Lambert, impulsive, flyby-optimized, timing-adaptive, etc.)
- Each strategy explores a different region of the solution space
- Results are merged, deduplicated, and ranked by objective function
- Adaptive tuning adjusts parameters based on intermediate results

**Optimization Objective:**
Minimize: Σ(Δv_i) + time_penalty + constraint_violations

Subject to: Solar distance limits, flyby safety margins, epoch windows for each body.

## Troubleshooting

- **Missing bodies.mat:** Run `SETUP_mythara_paths_auto.m` or copy from `../data/bodies.mat`
- **Path errors:** Ensure you're running from the project root or `optimization/submission/`
- **Performance too slow:** Use `Main_Flyby_multistrategy_LIGHT.m` for testing; reduce time windows or search depth
- **Unexpected validator errors:** Check `orbit_calibration/CALIBRATION_NOTES.md` and re-validate ephemeris

## Results

Final candidate solutions are written to:
```
../csv/topscore_J*.csv
../csv/runnerup_J*.csv
../csv/third_J*.csv
../csv/grandtour_J*.csv
```

Intermediate results and statistics appear in the MATLAB console and (optionally) log files.

## References

- Lambert solver: Universal formulation (Bate, Mueller, White)
- KM propagator: Keplerian motion with numerical integration
- Flyby optimization: Sequential delta-v minimization with timing constraints
- Project report: See `../docs/PROJECT_REPORT.md`
