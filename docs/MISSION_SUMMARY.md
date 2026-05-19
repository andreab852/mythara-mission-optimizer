# Mission Summary

This document summarizes the official Mythara project material included in
`docs/mission/`.

## Mission Objective

Design a long-duration spacecraft tour of the Mythara system using only gravity
assists. The spacecraft must arrive from an interstellar asymptote and explore the
system without propulsive maneuvers during the planetary tour.

The objective is to maximize the scientific score `J`, which depends on:

- the scientific value of each visited body;
- the number of scientific flybys per body, up to 13 counted flybys per body;
- the hyperbolic excess velocity at each encounter;
- seasonal diversity, measured from the heliocentric direction of repeated flybys;
- a Grand Tour bonus if every relevant body is visited scientifically.

## Body Weights

| ID | Body | Weight | Notes |
| ---: | --- | ---: | --- |
| 4 | Cindrel | 3 | Inclined terrestrial planet |
| 1000 | Nyxar | 5 | Massless dwarf planet in the trajectory model |
| 5 | Aurelios | 7 | Ringed gas giant |
| 6 | Tharos | 10 | Massive super-Jovian gravity-assist body |
| 7 | Glacior | 15 | Ice giant |
| 8 | Ashkara | 20 | Stripped dense planetary core |
| 9 | Varuun | 35 | Captured retrograde Jovian planet |
| 10 | Xerion | 50 | Highly eccentric, high-value outer giant |

## Grand Tour Bonus

The score multiplier is:

- `b = 1.2` if the trajectory includes scientific flybys of all planets and Nyxar;
- `b = 1.0` otherwise.

## Initial State

The spacecraft starts near an incoming interstellar asymptote:

- `X = -200 AU`;
- `Y` and `Z` are free;
- `Vx` is free;
- `Vy = 0`;
- `Vz = 0`;
- the start epoch must lie within the 300-year mission window.

The positive X direction is aligned with the incoming velocity direction toward
Mythara. The coordinate frame is the same frame used by the provided ephemerides.

## Constants

| Constant | Value |
| --- | ---: |
| AU | `149597870.691 km` |
| Mythara gravitational parameter | `139348062043.343 km^3/s^2` |
| Day | `86400 s` |
| Year | `365.25 days` |

## Trajectory Constraints

- All events must occur within 300 years from the reference epoch.
- Flyby positions must match the target-body positions within validator tolerance.
- Planetary flybys use a patched-conics gravity-assist model.
- Incoming and outgoing `V_infinity` magnitudes must match at massive bodies.
- Flyby altitude must be between `0.1` and `100` body radii above the surface.
- Nyxar is massless, so the incoming and outgoing `V_infinity` vectors must be
  continuous in both magnitude and direction.
- If two successive flybys are of the same body, their time separation must be at
  least one third of that body's orbital period.
- Heliocentric perihelion must stay above `0.05 AU`, except for one allowed thermal
  protection passage down to `0.01 AU`.

## Flyby Model

Flybys are modeled as zero-time events:

- incoming and outgoing flyby positions are identical;
- incoming and outgoing flyby epochs are identical;
- inertial velocity is the planet velocity plus the relative `V_infinity`;
- the turn angle is limited by the patched-conics hyperbolic passage and the allowed
  periapsis altitude range.

## CSV Submission Format

The converter writes a 12-column CSV:

```text
body_id, flag, epoch_sec, px, py, pz, vx, vy, vz, cx, cy, cz
```

where:

- `body_id` is the body identifier for flyby rows and `0` for conic-arc rows;
- `flag` is `1` for flyby rows and `0` for conic-arc rows;
- epoch is expressed in seconds from `t = 0`;
- positions are in kilometers;
- velocities are in kilometers per second;
- the control vector stores `V_infinity` for flyby rows and zeros for conic arcs.

The generated CSV is uploaded to the online validator, which checks continuity,
flyby feasibility, perihelion constraints, and the final score.
