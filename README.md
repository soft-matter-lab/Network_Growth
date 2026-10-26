# Network Growth

MATLAB Case-A stochastic network growth with periodic boundaries. Add `SRC` to
the MATLAB path, or use the drivers which do this automatically.

## Drivers

- `run_pkg_cleaned.m`: one realization; progress dumps are enabled.
- `run_pkg_cleaned_sweeps.m`: new batch driver with explicit seed/composition modes.
- `run_pkg_cleaned_batch.m`: older custom-override driver, retained separately.
- `visualize_pkg_cleaned.m`: offline driver; point it at saved dumps to export a
  network movie and distribution figures without running the simulation.

Edit the new driver's settings before running it. Its default plan is 15 runs:
five compositions, repeated for three seeds, at fixed `fA = 4`.

```matlab
sweep.mode = 'replicates';             % fixed baseCfg.xA, every seed
sweep.mode = 'composition';            % sweep.xA, only sweep.seeds(1)
sweep.mode = 'replicated_composition'; % sweep.xA for every seed
sweep.xA = [.01 .025 .05 .075 .10];
sweep.seeds = [1 2 3];
```

Choose **one** mode assignment. `xA` is the feed mole fraction A/(A+B);
the A/B ratio is `xA/(1-xA)`. Replicate blocks use distinct seeds; within a
block the same seed is reused at each composition. Different compositions still
have different event sequences. Seed values are saved in the plan and run config.

`sweep.NtotalScaleByXA = []` uses `baseCfg.NtotalScale` for every composition.
Optionally set one scale per `sweep.xA`, e.g. `[1.175 1.25 1.275 1.375 1.55]`
for the supplied fA=4 calibration. This option requires computed feed
(`baseCfg.Ntotal = []`) and is disabled in seed-only mode. Leave the vector empty
to compare compositions at a common scale. The fill-mode repulsion radius still
uses the **nominal, unscaled** feed count.

Preview without running:

```matlab
plan = makeNetworkSweepPlan(baseCfg,sweep);
disp(plan)
[summary,batchDir] = runNetworkSweep(baseCfg,sweep,outputRoot);
```

Every batch and run gets a fresh folder. `run_plan.csv` and `batch_config.mat`
record the plan; each run has `run_config.mat`, optional `results_struct.mat`, and
its dumps. `batch_summary.csv/.mat` update before and after each run, with status,
seed, composition, counts, final conversion, termination reason and output path.
`completed` means the simulation returned normally, including a chemically
jammed state; inspect `stopReason` and `usagePercent`. Failed runs retain their
completed snapshots and an error report. `sweep.continueOnError` controls whether
later runs proceed after a MATLAB error. Results are not accumulated in memory.

## Progress dumps

The library default is `saveDumps = false`; both the single-run and new sweep
drivers explicitly enable it. These options apply to `cfg` (single run) or
`baseCfg` (batch):

| Option | Meaning |
| --- | --- |
| `saveDumps = true` | Enable incremental MAT snapshots in `outputDir/dumps` |
| `dumpUsagePercent = 0:5:100` | Crosslinker-consumption milestones; `[]` disables this trigger |
| `dumpEveryEvents = 0` | 0 disables; 1 captures every reaction; N captures every N reactions |
| `dumpTimeInterval = 0` | 0 disables; positive value enables SSA-time interval sampling |
| `dumpSaveDistributions = true` | Raw completed and dangling chain lengths |
| `dumpSaveNetwork = true` | Node positions/degrees, bonds with periodic images, active arms |
| `dumpSaveFinal = true` | Separate snapshot after final mechanical relaxation |

At least one payload must be enabled. An initial snapshot is always saved when
dumps are enabled. Usage is `100*(initialFreeA - nA)/initialFreeA`, where
`initialFreeA = round(xA*Ntotal)-1`: the initial seed is excluded, so the start is
0%. Consumption includes crosslinker attachment and, when enabled, nucleation.
If there are no initially free crosslinkers, usage is undefined (`NaN`); only the
initial 0% marker is recorded, and no 100% milestone is invented.

States are captured **after a reaction and its mechanical relaxation**. A finite
reservoir may jump over several targets in one reaction; those targets share one
snapshot. Each snapshot stores `usageTargets` and the actual `usagePercent`.
There is no interpolated network. Event, time and usage triggers are combined,
so coincident triggers do not duplicate a state. Time sampling saves the first
post-event state crossing an interval; `timeBinRange` records the inclusive
indices of all intervals crossed in that time jump. This is SSA time, not wall
time or mechanical integration steps.

A jammed run can stop before 100%; unreached milestones stay absent. The 100%
checkpoint is the first state with no free A remaining. Chain growth/joining may
continue afterward, so the separately marked `phase = 'final'` snapshot can have
different distributions at the same conversion. Snapshot files use a temporary
`.part` file before being renamed; completed snapshots can be read after an
interruption. A folder containing existing dumps is never silently overwritten.

`metadata.mat` stores the requested and resolved configuration, initial counts,
RNG state at initialization, repulsion radius and column definitions. Each
`snapshot_000001.mat`, etc. contains one `snapshot` struct:

- `completed.kuhnCount`: n for each completed bond; `contourLength = n*b`;
  `endToEnd = |x_j + image*L - x_i|`, using the stored bond image.
- `dangling.kuhnCount` and `dangling.contourLength`: growing chains with n>0;
  `bareArmCount` counts n=0 arms separately. The model represents dangling ends
  as clouds anchored at their parent, so it has no explicit dangling tip position
  or measured end-to-end distance to save.
- Optional `nodes = [x_unwrapped,y_unwrapped,degree]`,
  `bonds = [i,j,n,imageX,imageY]`, and
  `activeEnds = [x_wrapped,y_wrapped,parent,n]`. These are junction and connectivity
  data, not explicit monomer coordinates. Rows give the raw distributions without
  committing to histogram bins.
- Event number, SSA time, reservoir counts, requested/actual usage, trigger
  information, phase and final termination reason.

Every-event full-network dumps can create many files; milestone-only capture is
the driver default. Set `saveFinalDistributions = false` to skip the legacy final
PNG/CSV export; it is independent of dumps, `saveData`, and `showFinalHistogram`.

## Offline visualization and analysis

**Start with `visualize_pkg_cleaned.m`, next to the run drivers.** Set
`inputFolder` to one run folder or its `dumps` subfolder, then run the script.
Leaving it empty opens a folder picker. Choose an individual run inside a batch,
not the batch's top-level directory.

```matlab
inputFolder = 'E:\...\runs_pkg_cleaned\single_xA_...';
```

The script's `vizOptions` controls these exports:

- `makeVideo`: `network_evolution.mp4`, replaying clouds, periodic bonds and
  junctions alongside reservoir histories. `videoFrameEvery` selects every Nth
  **saved snapshot** (the last is always included); `videoFPS` and
  `videoHoldLastSeconds` control playback. `videoPreview` shows rendering live.
- `saveCheckpointPlots`: one PNG per distinct saved state matching
  `checkpointUsagePercent` (default `0:5:100`), with network, Kuhn-count and
  end-to-end histograms. Several crossed milestones can share
  one state/PNG. `checkpoint_index.csv` maps the files to requested/actual usage.
- `saveFinalPlot`: a separate `final_state.png` for the final relaxed network.
- `saveComparison`: `distribution_comparison.png`, overlaying Kuhn-count,
  physical contour-length and completed-chain end-to-end distributions using
  shared bins. `comparisonUsagePercent` defaults to `[0 25 50 75 100]`; set it to
  `0:5:100` for all checkpoints. `population` chooses the contour population and
  `showComparisonFigure` controls whether the comparison stays open.

Checkpoint and final-state figures contain three panels: the network on the
left, contour length (Kuhn count) at upper right, and end-to-end distance at lower
right. Histogram **bins and x/y limits stay fixed across all snapshots**. The
visualizer scans every saved state, including final relaxation, once for ranges
and once for peak bin probabilities. Each histogram has its own fixed probability-axis limits.
Each bin height is its count divided by the total sample count for that plot in
that snapshot; the title reports N. With the full range, heights sum to one.
This is probability per bin, not density (which additionally divides by bin
width). Normalization only rescales heights and preserves within-snapshot shape
for unchanged bins. Empty populations are marked "No samples". Explicitly
cropped ranges retain the original denominator, so visible probabilities can
sum to less than one. Raw dump samples are unchanged.
Set `vizOptions.kuhnXLim`, `vizOptions.endToEndXLim`, and
`vizOptions.kuhnYLim` / `vizOptions.endToEndYLim` to explicit two-element ranges to compare different
runs; `[]` chooses common ranges automatically within the selected run.
`endToEndBins` sets the shared end-to-end bin count. Resolved bins/limits are
saved in `visualization_report.mat` as `report.histogramStyle`.

Outputs go to a new `postprocessing/visualization_<timestamp>` folder under the
run. Set `vizOptions.outputRoot` to use another disk or directory. Saved dumps
are unchanged; rerunning the visualizer creates another output folder.
`video_frames.csv` records each snapshot's video-frame range, SSA event/time and
conversion. Settings and an output report are also saved as MAT files.

The movie advances one stored state per frame at equal playback intervals; it
does not interpolate SSA time or reconstruct unsaved events or mechanical
substeps. With only 5% checkpoints the movie is correspondingly coarse. Capture
more frequently with `dumpEveryEvents` or `dumpTimeInterval` for future runs.
Video requires the network payload; distribution-only dumps still produce
distribution figures, with the video skipped. Missing checkpoints/final states
are reported and skipped. The movie renderer streams selected snapshots instead
of holding every network in memory.

For programmatic use, call `postprocessNetworkDumps(inputFolder,vizOptions)`.
The lower-level functions remain available:

```matlab
addpath('SRC')
[snapshots,metadata,index] = loadNetworkDumps(runDir); % run or dumps folder
disp(index)
n = snapshots{end}.completed.kuhnCount;
r = snapshots{end}.completed.endToEnd;

% Shared histogram bins, normalized to probability per bin.
[fig,selection] = plotNetworkDumps(runDir, ...
    struct('usagePercent',[0 25 50 75 100]));
% Use 0:5:100 to compare all recorded milestones.
plotNetworkDumps({runDir1,runDir2},struct('usagePercent',[25 50 75 100], ...
    'labels',{{'seed 1','seed 2'}},'outputFile','comparison.png'));

% Inspect a single saved state with shared histogram axes across the run.
style = networkDumpHistogramStyle(runDir);
fig = plotNetworkDump(snapshots{end},struct('histogramStyle',style));
% A filename inside a dumps folder automatically uses that run's shared axes.
plotNetworkDump(fullfile(runDir,'dumps','snapshot_000010.mat'));
```

The comparison plots Kuhn counts, physical contour lengths, and completed-chain
end-to-end distances. `population` selects `'completed'` (default), `'dangling'`,
or `'all'` for the contour panels; the end-to-end panel always uses completed
chains. Empty populations are labeled, not replaced with zero-length chains.
`selection` reports files, events, actual conversions and sample counts.
Only saved milestone markers are selected; missing targets produce a warning.
Multiple targets can refer to the same file in small systems. Replicates are
shown separately, without automatic pooling or averaging.

Optional `nEdges`, `contourEdges`, `rEdges` fix common bins across comparisons;
`visible = 'off'` and `outputFile = 'plot.png'` support noninteractive export.
`loadNetworkDumps(runDir)` loads all snapshots. The comparison plot retains only
the requested milestones; it scans other files one at a time. You can also use
`loadNetworkDumps(runDir,[25 50 75 100])` to retain selected milestone states, or
`plotNetworkDump(filename)` / MATLAB `load(filename,'snapshot')` for one frame.

`networkRdf(nodes,L,edges)` computes the 2D minimum-image **junction-junction**
RDF, including bonded pairs, using unordered pair counts and annulus areas.
Edges must lie in `[0,L/2]`; normalization is
`N*(N-1)/2 * pi*diff(edges.^2)/L^2` (ideal uniform points give g(r) near 1).
Fewer than two nodes gives undefined (`NaN`) g(r). RDF is computed offline from
positions; it is not a monomer RDF or a relaxation/convergence diagnostic.

## Replicate-averaged evolution plots

Run `average_pkg_cleaned.m` beside the other drivers. Enter the replicate run
folders in `runFolders` (e.g. five runs with distinct seeds), or leave that list
empty and set `batchFolder` and `selectedXA` to select a composition from a sweep.
The averager checks that all non-output settings match and seeds are distinct.

`averageOptions.usagePercent = [25 50 100]` selects evolution points. Completed
chain histograms use common bins across all runs and points. Each replicate is
normalized to probability per bin **before** averaging, giving each run equal
weight rather than giving runs with more chains more weight. The three panels
show Kuhn count, physical contour length, and end-to-end distance. At fixed b,
the first two describe the same distribution in different units.

Shading (`showSpread`) is plus/minus one sample standard deviation across
replicates, clipped to [0,1] for display; it is not a confidence interval. Raw SDs
are saved without clipping. Averaging reduces sampling fluctuations but does not
guarantee a smooth curve. No smoothing kernel is applied. Increase `nBinWidth`
or reduce `distanceBins` for coarser, explicitly controlled histograms.

By default all replicates must have reached each checkpoint and have a nonempty
completed-chain population. This also means the empty 0% distribution cannot be
averaged. Missing checkpoints are never replaced with a later or final state.
`requireAllReplicates = false` explicitly permits available-run averages; the
plot reports the contributing count at each point, and missing/empty populations
are excluded rather than treated as zero probability. Such averages can have
selection bias when some runs stop early. A single contributing run has undefined
SD. Actual conversion (which can overshoot a target in a small system), seed,
sample count and snapshot filename are recorded for every replicate/target.

Outputs are written to a unique folder under `replicate_averages`: a PNG, a MAT
file containing means, SDs, bin edges and per-replicate histograms, and a CSV
selection table. Edit `outputRoot` in the driver to choose another destination.
Programmatic use: `[fig,data] = averageNetworkDumps(runFolders,averageOptions)`.

## Verification

```matlab
results = runtests('tests');
assertSuccess(results)
```

Tests check seeded trajectory invariance with dumping on/off, raw segment
conservation, conversion and time/event sampling, periodic distances, RDF
normalization, sweep plans and nominal-feed radius preservation, and a small
four-run batch with offline figure exports. Visualization tests also export and
decode an MP4, check frame selection and final-frame holds, and verify that the
source dumps are unchanged. Artifacts are written to a unique
MATLAB temporary folder; the test prints its path. No production sweep is run.
