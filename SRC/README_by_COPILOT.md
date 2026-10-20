# pkg_cleaned

Cleaned modular layout for the `NetworkA_single_v2` simulation.

## Run

1. Open and edit `run_pkg_cleaned.m` for your chosen inputs (`xA`, `Ntotal`, `fA`, plotting/video flags, etc.).
2. Run `run_pkg_cleaned.m`.

The script calls `runNetworkA_cleaned(cfg)` and returns a `results` struct in the workspace.
Each invocation creates a fresh run folder under `runs_pkg_cleaned/`, so outputs are never overwritten.

## Batch run

1. Open and edit `run_pkg_cleaned_batch.m`.
2. Set `runOverrides` with per-run parameter changes (for example different `xA` values).
3. Run `run_pkg_cleaned_batch.m`.

Batch outputs are saved under a unique folder in `batch_runs_pkg_cleaned/`:
- one unique subfolder per run,
- per-run data/movie/results files,
- `batch_summary.csv` and `batch_summary.mat` at batch root.

## File roles

- `run_pkg_cleaned.m`: user-facing main driver (edit inputs here).
- `run_pkg_cleaned_batch.m`: batch driver for multiple runs with per-run overrides.
- `defaultNetworkAConfig.m`: defaults for all tunable parameters.
- `runNetworkA_cleaned.m`: SSA loop + event execution + output reporting.
- `createUniqueRunFolder.m`: unique folder creator for non-overwriting run outputs.
- `equilibrateNetwork.m`, `networkForces.m`, `networkEnergy.m`: mechanics.
- `initVisualization.m`, `updateLivePlot.m`, `drawMotion.m`: live plotting.
- `buildHeatmap.m`, `buildSkeleton.m`, `wrapPoint.m`, `bondVector.m`: geometry/visual helpers.
- `computeDegrees.m`, `makeAdjacency.m`, `validateNetworkState.m`: topology utilities.

## Notes

- Physics/chemistry behavior is preserved from the current v2 model:
  - `b*n^0.6` new-junction placement.
  - tensile-only bonded spring.
  - harmonic repulsion with feed-based constant `rEq` (including bonded pairs by default).
  - periodic bond-image tracking.

- Equilibration now supports two modes via config:
  - `eqIntegrator = 'damped2nd'` (default, mass+damping),
  - `eqIntegrator = 'overdamped'` (legacy first-order).

- Performance note:
  - `useRepulsionBinning` is available, but for current MATLAB-scale network sizes
    the vectorized all-pairs path is typically faster, so default is `false`.
