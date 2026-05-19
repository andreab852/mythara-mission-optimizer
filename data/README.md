# Data Files

This folder contains the canonical planetary ephemeris and supporting datasets required by the Mythara solver.

## Files

### `bodies.mat`

Binary MATLAB data file containing the Mythara system planetary and moon ephemeris.

**History:**
- Original ephemeris (from GTOC13 data) produced inconsistent CSV outputs due to coordinate system transformations
- **Fixed via dedicated validator** that reconstructs the ephemeris with corrected orbital elements and heliocentric frames
- This corrected version is the source of truth for all trajectory computations in the optimization folder

**Content:**
- 10 planets (Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto-equivalent)
- 1 special body: Nyxar (ID 1000) — target of exploration in extended tour scenarios
- Orbital elements: Semi-major axis, eccentricity, inclination, RAAN, argument of perigee, mean anomaly
- Ephemeris interpolation tables for rapid position/velocity lookup

**Usage in MATLAB:**
```matlab
load('data/bodies.mat', 'bodies_struct');
% Access ephemeris for body 5 (Jupiter) at Julian Date jd
planet_pos = bodies_struct.get_position(5, jd);
```

**Locations:**
- Canonical: `data/bodies.mat`
- Copy: `optimization/submission/bodies.mat` (used by Main_Flyby scripts)

### `MytharaSystemData.csv`

Reference CSV containing static properties of the Mythara system bodies:
- Body ID, Name, Orbital period, Semi-major axis, Eccentricity, etc.

Used for validation and documentation purposes. The working ephemeris is in `bodies.mat`.

## Validation

To verify bodies.mat integrity:

```matlab
load('data/bodies.mat');
% Check that all 11 bodies (1–10 plus 1000) are present
assert(length(bodies_struct.bodies) == 11, 'Incomplete ephemeris');
```

If you suspect ephemeris corruption or need to rebuild from raw GTOC13 data, see `optimization/orbit_calibration/` for the validator and calibration scripts.
