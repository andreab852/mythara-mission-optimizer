# Project Report

## Overview

The Mythara Grand Tour project is an astrodynamics optimization problem. The task is
to construct a feasible multi-flyby trajectory through an artificial exoplanetary
system, using gravity assists only, and to maximize the scientific score returned by
an online validator.

The final repository contains:

- the corrected planetary model used by the final scripts;
- the MATLAB trajectory search code;
- the CSV exporter used for validator submissions;
- preserved final CSV candidates;
- documentation of the mission rules and calibration process.

## Problem Interpretation

The spacecraft begins from an incoming interstellar state near `X = -200 AU`. Its
trajectory must remain inside a 300-year mission window. Each transfer leg is modeled
as a heliocentric conic arc solved through Lambert's problem, and each encounter is
modeled as a zero-duration patched-conics flyby.

The search therefore has to satisfy three layers of constraints:

1. Lambert transfers between body positions at selected epochs.
2. Flyby feasibility through `V_infinity` matching and turn-angle limits.
3. Validator-level constraints such as perihelion limits, repeated-flyby timing,
   Nyxar vector continuity, and numerical tolerances.

## Scoring Strategy

The score favors repeated scientifically valuable flybys, especially of high-weight
bodies such as Varuun and Xerion. However, repeated flybys of the same body are
penalized if they occur at similar heliocentric seasonal geometry. High flyby speed
is also penalized.

The search therefore balances:

- visiting high-weight bodies;
- preserving seasonal diversity;
- keeping `V_infinity` in favorable ranges;
- avoiding infeasible flyby turns;
- staying inside the 300-year mission window.

## Search Implementation

The final operational script is:

```text
optimization/submission/Main_Flyby_multistrategy.m
```

It performs a multi-strategy beam search over body sequences, transfer times, and
Lambert branches. The version exposes the main exploration parameters and
parallel-worker settings directly in the `default_opts` section.

The supporting astrodynamics code is in:

```text
optimization/flyby_functions/
optimization/submission/utils/
```

## Calibration Work

The original data model did not reproduce validator-consistent CSV trajectories for
all tested cases. To solve this, the team used validator feedback to infer corrected
body states and rebuild a final `bodies.mat`.

The calibration workflow is preserved in:

```text
optimization/orbit_calibration/
docs/CALIBRATION_NOTES.md
optimization/README_bodies_calibration.txt
```

## Final Outputs

The final candidate/submitted CSV trajectories are stored in:

```text
csv/
```

These files are kept intentionally. Heavy intermediate archives and personal
development folders are not included in this cleaned repository.
