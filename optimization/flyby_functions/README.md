# Flyby Functions Library

Core orbital mechanics and trajectory optimization functions used by the solver.

## Key Functions

### Trajectory Propagation

- **`KM_solver.m`** — Kerbal Motion numerical integrator
  - Solves the two-body problem with high precision
  - Input: Initial state (r, v), epochs
  - Output: State vector at target epoch
  - Used in validation loops and trajectory re-propagation

### Lambert Solver

- **`lambert_universal.m`** — Universal variable formulation
  - Solves the Lambert boundary value problem
  - Input: Initial position, final position, transfer time, gravitational parameter
  - Output: Required velocity at departure and arrival
  - Supports multiple revolution solutions (e.g., long/short arcs)
  - Basis for all impulsive transfer planning

### Flyby Optimization

- **`flyby_opt.m`** — Multi-impulse flyby delta-v optimizer
  - Minimizes total Δv for a sequence of impulsive maneuvers
  - Handles gravity assists and swing-by constraints
  - Input: Sequence of bodies, epochs, max Δv per impulse
  - Output: Optimal impulse vector and timing

### CSV & Format Utilities

- **`solution_to_csv_patched.m`** — Trajectory to GTOC13 CSV
  - Converts trajectory state vectors to official submission format
  - Handles epoch conversions and coordinate frame transformations
  - Validates constraint compliance before export
  - Input: Trajectory struct with bodies and epochs
  - Output: GTOC13-compliant CSV file

### Orbital Element Tools

- `orbital_elements_from_state.m` — State → Keplerian elements
- `state_from_orbital_elements.m` — Keplerian elements → State
- `true_anomaly_at_epoch.m` — Epoch-dependent anomaly calculation
- `frame_transform_*.m` — Heliocentric ↔ relative frame conversions

### Constraint & Validation

- `check_solar_distance.m` — Verify solar distance limits
- `check_flyby_safety.m` — Minimum safe altitude at flyby
- `check_epoch_feasibility.m` — Window and timing constraints
- `validate_solution.m` — Full trajectory validation before CSV export

## Usage Pattern

Typical workflow in the main solver:

```matlab
% 1. Generate initial guess (random or heuristic sequence)
sequence = [3 5 2 4 1];  % Body IDs to visit

% 2. For each pair, compute transfer arc
r_from = bodies(sequence(i)).position(epoch1);
r_to = bodies(sequence(i+1)).position(epoch2);
[v_from, v_at_arrival] = lambert_universal(r_from, r_to, tof, mu_sun);

% 3. Propagate trajectory through flyby
state = KM_solver([r_from; v_from], epoch1, epoch2);

% 4. Apply gravity assist impulse at flyby (via flyby_opt)
[dv_optimal, v_out] = flyby_opt(state.v, state.r, body_params, constraints);

% 5. Export solution
solution_to_csv_patched(trajectory, 'solution.csv');
```

## Function Interdependencies

```
Main_Flyby_multistrategy_*.m
├─ lambert_universal.m
├─ KM_solver.m
├─ flyby_opt.m
├─ validate_solution.m
│  ├─ check_solar_distance.m
│  ├─ check_flyby_safety.m
│  └─ check_epoch_feasibility.m
├─ solution_to_csv_patched.m
│  └─ frame_transform_*.m
└─ orbital_elements_from_state.m
```

## References

- **Lambert solver:** Bate, Mueller, White (1971) — *Fundamentals of Astrodynamics*
- **KM integrator:** Runge-Kutta with adaptive step control
- **Flyby optimization:** Successive approximation with penalty methods
- **CSV format:** GTOC13 official problem statement

## Notes

- All functions work in SI units (meters, seconds, m/s)
- Gravitational parameters use standard solar system constants
- No external libraries required (native MATLAB)
- Performance optimized for batch trajectory generation
