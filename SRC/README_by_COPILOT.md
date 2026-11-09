# pkg_cleaned

Shared source for the growing-chain (Case A) and preformed-chain (Case B)
simulations. See the project-root `README.md` for the full driver map and Case B
settings. Legacy root scripts are archived in `tmp/legacy_root_drivers_20261109`.

## Run

1. Open and edit `drivers/growing_chains/run_single_growing.m` for your chosen inputs (`xA`, `Ntotal`, `fA`, plotting/video flags, etc.).
2. Run that script. For Case B, use `drivers/preformed_chains/run_single_preformed.m`.

The script calls `runNetworkA_cleaned(cfg)` and returns a `results` struct in the workspace.
Each invocation creates a fresh run folder under `runs_pkg_cleaned/`, so outputs are never overwritten.

## Batch run

1. Open and edit `drivers/growing_chains/run_sweep_growing.m` or `drivers/preformed_chains/run_sweep_preformed.m`.
2. Set the sweep mode, compositions, and seeds (plus optional Case B parameter sets).
3. Run the selected script.

Batch outputs are saved under a unique folder in `batch_runs_pkg_cleaned_sweeps/`
(Case A) or `batch_runs_preformed/` (Case B):
- one unique subfolder per run,
- per-run data/movie/results files,
- `batch_summary.csv` and `batch_summary.mat` at batch root.

## File roles

- `drivers/growing_chains/`, `drivers/preformed_chains/`: editable run and postprocessing drivers.
- `defaultNetworkAConfig.m`, `defaultNetworkBConfig.m`: defaults for each route.
- `runNetworkA_cleaned.m`, `runNetworkB_cleaned.m`: route entry points.
- `runNetwork_cleaned.m`: shared SSA loop, event execution, and output reporting.
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
