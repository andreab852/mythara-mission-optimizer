# Submission: Final Solver & Results

This folder contains the **ready-to-run multi-body trajectory optimizer** with all dependencies included.

## Main Scripts

### `Main_Flyby_multistrategy.m`
**Primary entry point for the final search.**

Implements a 16-strategy parallel evolutionary search for multi-body trajectories:
- **16 parallel search strategies** exploring different regions of the solution space simultaneously
- **Adaptive tuning** adjusts strategy parameters based on intermediate results
- **Multi-body feasibility** enforces solar distance, flyby safety, and time window constraints
- **CSV output** exports top candidates in GTOC13 format

**Run time:** ~1–4 hours depending on hardware and search depth.

**Example:**
```matlab
Main_Flyby_multistrategy
% Outputs: CSV files with candidates ranked by objective J
```

### `Main_Flyby_multistrategy_LIGHT.m`
**Lightweight prototype** with reduced search scope. Completes in ~5–10 minutes. Use for algorithm validation or quick feasibility checks.

## Dependencies

- `bodies.mat` — Planetary ephemeris (canonical source: `../data/bodies.mat`)
- `utils/` — Utility functions for trajectory generation and validation
- `Correct/` — CSV export and post-processing routines
- Core functions from `../flyby_functions/`

## Quick Start

```matlab
% Option 1: From project root
cd('D:\path\to\mythara-matlab-project');
run_final_search  % Executes Main_Flyby_multistrategy

% Option 2: Direct run
cd('optimization/submission');
Main_Flyby_multistrategy
```

## Outputs

CSV files in the working directory:
- `topscore_J*.csv` — Best candidate
- `runnerup_J*.csv` — Second best
- `third_J*.csv` — Third best
- Plus all other validated candidates ranked by objective

See `../../csv/README.md` for interpretation of score (J) and file format.

## Algorithm Parameters

Typical tuning parameters (can be modified in the script):
- **Search depth:** Number of generations or iterations per strategy
- **Population size:** Trajectories per generation
- **Time windows:** Min/max epochs for body visits
- **Δv budget:** Maximum mission delta-v
- **Flyby altitude:** Minimum safe distance from body surface

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `bodies.mat not found` | Copy from `../data/bodies.mat` or run `SETUP_mythara_paths_auto.m` |
| Missing functions in `utils/` or `Correct/` | Verify paths are set correctly; check `../orbit_calibration/SETUP_mythara_paths_auto.m` |
| Very slow performance | Use `Main_Flyby_multistrategy_LIGHT.m` or reduce search depth parameter |
| Validator errors in output CSV | Check `../../docs/CALIBRATION_NOTES.md` and bodies.mat ephemeris |

## Performance Notes

- **Parallel execution:** Strategies run concurrently; speedup depends on available CPU cores
- **Memory:** Typical run uses 2–8 GB RAM for large solution libraries
- **GPU acceleration:** Available in specialized variants (see `../flyby_functions/`)

For reproducibility or benchmarking, record the random seed or strategy order used.
