# Progressive MDOF Building Degradation Simulator

A MATLAB-based research/demo code that simulates nonlinear, progressive stiffness degradation in multi-degree-of-freedom (MDOF) shear-building models under earthquake ground motions. The code assembles per-floor mass and stiffness from simple geometry/material inputs, runs a Newmark time‑integration solver with progressive damage (compiled C MEX), and writes detailed per-building results (MAT files, PNG plots, and human-readable summaries). A project-level summary is produced after all simulations.

## Highlights
- Multistory shear-building models with per-floor mass and stiffness computed from simple geometry.
- Nonlinear progressive stiffness degradation controlled by yield/ultimate drift, residual strength and degradation rate.
- High-performance time integration implemented in C as a MEX file (`mdof_degrade_mex.c`) for speed.
- Robust AT2 (PEER) ground-motion reader (`readAT2.m`) and project-summary writer.
- Outputs: MAT files with full histories, PNG plots, and per-building summary text. Consolidated `results/project_summary.txt`.

## Project structure (important files/folders)
- `src/` — MATLAB source + C MEX (`run_simulations.m`, `run_mdof_mex.m`, `mdof_degrade_mex.c`, helpers).
- `data/` — Place input ground motions (PEER `.AT2`) here. Preferred filename: `earthquake_data.AT2`.
- `results/` — Generated automatically; contains per-country and per-building folders with MAT, PNG, and `_summary.txt` files.
- `docs/` — Report templates and documentation.
- `bin/` — (optional) place precompiled MEX binaries here if you cannot compile locally.

Key MATLAB scripts:
- `run_simulations.m` — orchestrates locating the project root, compiling MEX (if needed), running all countries/buildings, and creating the consolidated summary.
- `countries.m` — defines Countries and building inventories (geometry, materials, nonlinear params); edit this to add or change buildings.
- `create_building.m`, `get_m.m`, `get_k.m` — helpers that compute floor mass and stiffness.
- `run_mdof_mex.m` — runs the MEX solver for a single building, postprocesses, plots and saves results.
- `create_project_summary.m` — generates `results/project_summary.txt` (compact table across buildings).

## Inputs
- Ground motion: PEER `.AT2` file(s). Place at least one `.AT2` file in `data/` (preferred name `earthquake_data.AT2`). `readAT2.m` reads `DT` and `NPTS` robustly and returns acceleration in units of g.
- Building definitions: edit `countries.m` to add or modify building inventories and geometric/material parameters.

## Outputs
For each building run you will find under `results/<Country>/<Building>/`:
- `<Building>_mdof_results_building<idx>.mat` — full MATLAB workspace variables (displacements, drifts, stiffness history, final_stiffness, perFloorStiffRed_pct, etc.).
- `<Building>_figure_*.png` — plots (stacked displacements, etc.).
- `<Building>_building<idx>_summary.txt` — human-readable summary & damage classification.
- `results/project_summary.txt` — consolidated table summarizing max displacement, average stiffness reduction, and damage state for all buildings.

### Report

The final report is in the docs folder saved as report.pdf

## Quick start — run everything
1. Open MATLAB and make sure your current folder is the project root (the folder that contains `src/`).
2. Ensure a supported C compiler is configured:
   - Run `mex -setup` in MATLAB and follow prompts.
3. Place at least one `.AT2` file in `data/` (or in project root as fallback).
4. From MATLAB command prompt run:
   - `run_simulations`
   This will attempt to compile `mdof_degrade_mex.c` (if not found), run simulations for every country/building defined in `countries.m`, and then create `results/project_summary.txt`.

Alternatively, from a system terminal (non-interactive):
- cd dir/to/projectfolder/src
- matlab -batch "run_simulations"

Notes:
- If automatic compilation fails, manually compile in MATLAB:
  - cd to `src/` and run `mex mdof_degrade_mex.c`
- If you cannot compile MEX (e.g., missing compiler), you can place a precompiled MEX binary for your platform into `src/` (name must match `mdof_degrade_mex` + platform extension, e.g. `mdof_degrade_mex.mexa64` on Linux).


## Runtime options and debugging
- Plots can be suppressed by `run_mdof_mex(building, country, idx, false)` (default is off in the batch runner).
- If the MEX solver reports singular Keff or unexpected results, check inputs:
  - Positive definite mass and initial stiffness (no zero/negative values).
  - Consistent vector sizes for `story_heights`, `wall_thickness`, `density`, `E`.
- Warnings and summary files are written into each building folder — consult them for diagnostics.

## Dependencies & compatibility
- MATLAB (recommended R2019b or newer). Basic functionality may work on older releases but MEX compilation behavior can vary.
- A working C compiler compatible with MATLAB MEX (e.g., GCC on Linux, Visual Studio on Windows, Xcode on macOS).
- Optional: exportgraphics (R2020a+) used for high-quality PNG export; fallback to `saveas` is present.
- Tested on Linux/Windows/macOS with appropriate MEX builds.

## Troubleshooting
- "Could not find any .AT2 files" — ensure `data/` contains `.AT2` or place an `.AT2` in the current folder.
- MEX compile errors — run `mex -setup` and ensure a native compiler is installed. See MATLAB docs for platform-specific setup.
- If `run_simulations` cannot find project root: run it from inside the project tree or ensure your project has a `src/` folder.

## Development notes / extension points
- The solver returns and saves `final_stiffness`, `perFloorStiffRed_pct`, and `max_drift_per_story` to make downstream summary robust.
- You can replace the ground motion reader or add other formats by modifying `run_mdof_mex.m`.
- The damage classification logic is in `run_mdof_mex.m` (summary writer) and `create_project_summary.m` (project-level classification); customize thresholds to match your fragility/damage model.

## Contributors
- Mathis Olivier Horst
- Munongo Kazembe Luc-Faustin

## License & citation
- Add a LICENSE file to define reuse and citation policy. If this is for publication, include a citation string and DOI here.

## Acknowledgments

Data sources:
- PEER ground motion archive (earthquake_data.AT2). PEER Center, University of California, Berkeley. https://peer.berkeley.edu/


Code and libraries:
- mdof_degrade_mex.c adapted from Structural Dynamics of Earthquake Engineering_ Theory and -- Sekaran Rajasekaran -- Woodhead Publishing Series in Civil and Structural 

Tooling / assistance:
- the general structure of mdof_degrade_mex.c was generated with the help of Open AI GPT5-mini
