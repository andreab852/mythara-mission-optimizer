# Documentation

This folder contains comprehensive documentation of the Mythara project, including mission background, technical reports, calibration history, and results analysis.

## Files

### `MISSION_SUMMARY.md`

High-level overview of the GTOC13 problem and the Mythara mission concept:
- Problem statement and constraints
- Mission objectives (visit bodies, minimize fuel, meet time windows)
- System topology and body properties
- Performance benchmarks

Start here for a quick understanding of the project scope.

### `PROJECT_REPORT.md`

Complete technical report covering:
- Algorithm design (Lambert universal solver, KM method for orbital propagation)
- Multi-body trajectory optimization strategy
- Parallel search architecture
- Validation and verification procedures
- Results summary and key findings

Suitable for academic or portfolio review.

### `CALIBRATION_NOTES.md`

Detailed explanation of the planetary ephemeris calibration process:
- Original GTOC13 data issues and anomalies
- Validator algorithm for reconstructing correct ephemeris
- Fix methodology: orbital element extraction, frame transformation, validation checks
- Bodies.mat generation and integrity verification

Essential reference for understanding why `data/bodies.mat` was rebuilt.

### `RESULTS.md`

Summary of top candidate solutions with story-tagged descriptions:
- **Topscore** (J=251.57): Best overall solution with optimal sequence and timing
- **Runnerup** (J=249.61): Alternative backbone with similar performance
- **Third** (J=241.27): Mid-tour optimization variant
- **Grand Tour** (J=183.94): Extended mission including Nyxar (special body 1000)

Each entry includes mission duration, sequence, and key characterization.

### Official PDFs

- `GTOC13_Problem_Statement.pdf` — Official competition problem description
- `Mythara_Final_Presentation.pdf` — Slide deck summarizing research and findings

## Navigation

- **For portfolio/hiring:** Start with `PROJECT_REPORT.md` and `RESULTS.md`
- **For technical deep-dive:** Read `CALIBRATION_NOTES.md` and `PROJECT_REPORT.md`
- **For mission background:** Start with `MISSION_SUMMARY.md`

## Cross-References

- Top solutions: See `../csv/README.md`
- Source code: See `../optimization/README.md`
- Run instructions: See `../README.md` and `../run_final_search.m`
