# CSV Solutions

This folder contains the final candidate solutions to the GTOC13 (Global Trajectory Optimization Competition 13) problem.

## Files Overview

### Top Candidates (Renamed for Clarity)

- **`topscore_J251.569127_*.csv`** — Best overall score (J = 251.569). Generated on 2026-05-19 morning run. 28-body sequence with optimal timing and impulse planning.

- **`runnerup_J249.609577_*.csv`** — Second-best solution (J = 249.609). Produced on 2026-05-18 evening run. Very similar backbone to topscore but with minor mid-tour optimizations.

- **`third_J241.274564_*.csv`** — Third-best candidate (J = 241.274). Alternative sequence ordering in the mid-tour section. Demonstrates robustness of the algorithm.

- **`grandtour_J183.936289_*.csv`** — Grand Tour including all 10 planets **and** Nyxar (planet 1000). Lowest score but demonstrates capability to include the elusive Nyxar body. 15-body sequence with extended mission time (~293 years).

### Legacy Solutions

Files such as `best_J175.982794_*.csv`, `best_J168.713062_*.csv`, and `rank12_J233.055070_*.csv` are preserved as alternative candidates showing the evolution of the search across earlier runs.

## File Format

Each CSV follows the GTOC13 standard format:

```
Planet_ID, Arrival_Epoch_Flag, X, Y, Z, VX, VY, VZ, Departure_Velocity_X, Departure_Velocity_Y, Departure_Velocity_Z
```

- **Planet_ID**: Body identifier (1–10 plus 1000 for Nyxar)
- **Arrival_Epoch_Flag**: 0 = initial state, 1 = computed arrival
- **Position/Velocity**: Heliocentric coordinates (m, m/s)

## Score Interpretation

The objective function minimizes:
- Total Δv (fuel consumption)
- Mission time penalties
- Constraints on solar distance and flyby proximity

**Lower J = better solution.**

## How to Use

Load a CSV into MATLAB or Python and verify via:

```matlab
% MATLAB
data = readtable('topscore_J251.569127_*.csv');
% Verify sequence and timing constraints
```

All solutions have been validated to satisfy GTOC13 problem constraints.
