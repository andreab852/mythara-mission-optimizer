# Mythara Gravity-Assist Grand Tour

This repository contains the final cleaned MATLAB project for the 2026 Astrodynamics
group assignment on the Mythara exoplanetary system.

The goal is to design propellant-free multi-flyby trajectories that maximize the
scientific score assigned by the online validator. The final workflow is centered on
a corrected `bodies.mat` ephemeris model, a multi-strategy flyby search, and CSV
exports ready for upload to the validator.

## Final Entry Point

Open MATLAB in this repository root and run:

```matlab
run_final_search
```

This launches:

```text
optimization/submission/Main_Flyby_multistrategy.m
```

The script contains the parallel version of the final search workflow. Its
main parameters are documented in the `default_opts` section of the script.

## Repository Layout

Each folder contains a dedicated `README.md` with detailed explanations:

```text
csv/                         Final candidate/submitted CSV trajectories
                             ├─ README.md (explains 4 top solutions with story tags)
data/                        Corrected final bodies.mat and original system CSV
                             ├─ README.md (explains ephemeris calibration history)
docs/                        Project documentation and official mission material
                             ├─ README.md (navigation guide to all documents)
                             ├─ MISSION_SUMMARY.md
                             ├─ PROJECT_REPORT.md
                             ├─ CALIBRATION_NOTES.md
                             └─ RESULTS.md
optimization/                MATLAB search, scoring, export, and calibration code
                             ├─ README.md (overview of solver architecture)
                             ├─ submission/ (entry point)
                             │  ├─ README.md (how to run the solver)
                             │  └─ Main_Flyby_multistrategy.m
                             ├─ flyby_functions/ (core astrodynamics)
                             │  └─ README.md (function reference and usage)
                             └─ orbit_calibration/ (ephemeris validation)
                                └─ README.md (calibration workflow guide)
```

**Key links:**
- Start here: [`csv/README.md`](csv/README.md) for top 4 solutions
- For portfolio reviewers: [`docs/PROJECT_REPORT.md`](docs/PROJECT_REPORT.md)
- To run the solver: [`optimization/submission/README.md`](optimization/submission/README.md)
- Understand the science: [`docs/MISSION_SUMMARY.md`](docs/MISSION_SUMMARY.md)

Historical development folders, personal sandboxes, backups, and heavy run archives
are intentionally not included here.

## Final Data Model

The final corrected ephemeris file is stored in three intentionally identical copies:

```text
data/bodies.mat
optimization/submission/bodies.mat
optimization/flyby_functions/bodies.mat
```

The search and CSV exporter use the copy inside `optimization/submission`. The copy
inside `data` is kept as the canonical project-level data artifact. The copy inside
`optimization/flyby_functions` keeps legacy helper scripts consistent with the final
model.

## Final CSV Files

The `csv/` folder contains the final candidate/submitted trajectories preserved from
the project root. These files are intentionally versioned and should not be treated as
temporary outputs.

## Mission Summary

See:

```text
docs/MISSION_SUMMARY.md
```

for a concise English summary of the official mission objective, scoring terms,
constraints, flyby assumptions, and CSV submission format.

## Body Calibration

During development, the original system data did not reproduce validator-consistent
CSV results. The final workflow therefore includes a calibration track used to infer
and rebuild a corrected `bodies.mat` one body at a time from validator feedback.

See:

```text
docs/CALIBRATION_NOTES.md
optimization/README_bodies_calibration.txt
optimization/orbit_calibration/
```

## Main MATLAB Files

```text
run_final_search.m
optimization/submission/Main_Flyby_multistrategy.m
optimization/submission/Main_Flyby_multistrategy_LIGHT.m
optimization/submission/utils/solution_to_csv_patched.m
optimization/flyby_functions/KM_solver.m
optimization/flyby_functions/lambert_universal.m
optimization/flyby_functions/compute_gtoc13_score.m
```

## Requirements

- MATLAB
- Optimization Toolbox recommended
- Parallel Computing Toolbox recommended for the tunable parallel script

No neural-network model or external online access is required to generate CSV files.
The online validator is only needed to score uploaded CSV submissions.

## Notes for Reviewers

This cleaned repository is meant to present the final project, not the full historical
development process. Experimental branches, personal sandboxes, old backups, and
large generated archives were deliberately removed to keep the project readable.
