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
  `checkpointUsagePercent` (default `0:5:100`), with network, chain-number and
  weight-fraction histograms, and 2D displacement. Several crossed milestones can share
  one state/PNG. `checkpoint_index.csv` maps the files to requested/actual usage.
- `saveFinalPlot`: a separate `final_state.png` for the final relaxed network.
- `saveComparison`: `distribution_comparison.png`, overlaying Kuhn-count,
  physical contour-length and completed-chain end-to-end distributions using
  shared bins. `comparisonUsagePercent` defaults to `[0 25 50 75 100]`; set it to
  `0:5:100` for all checkpoints. `population` chooses the contour population and
  `showComparisonFigure` controls whether the comparison stays open.

Checkpoint and final-state figures contain four panels: network (upper left),
chain-number fraction versus Kuhn count (upper right), chain weight fraction
versus Kuhn count (lower left), and
a 2D end-to-end displacement histogram (lower right).
Histogram **bins and x/y limits stay fixed across all snapshots**. The
visualizer scans the saved states, including final relaxation, to determine ranges
and peak bin probabilities. Each histogram has its own fixed probability-axis limits.
The number-fraction panel divides bin chain count by the number of completed
chains; its title reports N. The weight-fraction panel uses
`w(n) = n * numberOfChainsWithLength(n) / sum(completed.kuhnCount)`.
It answers: what fraction of completed-chain Kuhn segments belong to chains of
this length? For a bin spanning several lengths, the numerator is the sum of
the actual n values of all chains in that bin, not count times bin center.
The denominator includes **only completed chains**, excluding dangling chains,
free monomers and junction mass. This is segment-weighted fraction, with the
total segment count shown in the title. Longer chains contribute more weight.
Both distributions sum to one over the full range when their denominators are
positive; these are fractions per bin, not densities. Empty/zero-segment
populations are marked as unavailable. Explicitly cropped ranges retain the
full denominator, so visible fractions can sum to less than one.
Raw dump samples and the separate single-run overlay plots are unchanged.
Replicate-average figures are described below.
Set `vizOptions.kuhnXLim`, `vizOptions.weightXLim`, and
`vizOptions.kuhnYLim` / `vizOptions.weightYLim` to explicit two-element ranges to compare different
runs; `[]` chooses common ranges automatically within the selected run.
Both length panels use integer-centered bins. Resolved bins/limits are
saved in `visualization_report.mat` as `report.histogramStyle`.

The 2D panel translates the first endpoint of each completed bond to (0,0) and
bins the other endpoint at `dv = x_j + image*L - x_i`. Translation preserves
the bond's direction: vectors are neither rotated nor wrapped to a new minimum
image. There is one point per bond; reversing the stored endpoint order reverses
its sign. No mirrored points are added. Color is bin count divided by that
snapshot's total bond count (probability per 2D bin, not probability density).
The origin is marked, both axes have equal spatial scale, and bins, symmetric
coordinate ranges and color limits are held fixed across all saved states.
`displacementBins` sets bins per axis; `displacementLimit = R` fixes both axes to
`[-R,R]`; `displacementCLim` fixes the probability color range, e.g. `[0 .02]`.
Empty bins are white, with a sequential white/gold/orange/deep-red color map.
When `displacementCLim=[]`, the upper color limit is the
`displacementColorPercentile` percentile (default 98) of positive bin probabilities
across all saved states. This shared cap improves contrast when a few concentrated
early snapshots dominate the maximum. Lower the percentile (e.g. 95) for more
contrast, or use 100 for the full run-wide maximum; an explicit `displacementCLim`
overrides the percentile. Bins above the cap use the highest color, and a
greater-than-or-equal top colorbar label identifies saturated panels. The raw
probabilities are unchanged, with no smoothing or per-frame color rescaling.
Empty coordinate limits use run-wide automatic values. Cropping retains the
full sample denominator. The chosen color limits are saved in `report.histogramStyle`.

New distribution payloads include `completed.endToEndVector` (columns dx, dy).
Older snapshots can reconstruct these vectors from nodes, bonds and periodic
images. Older distribution-only dumps contain magnitudes without directions;
the panel reports unavailable data rather than inventing directions.

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
One run is also accepted, e.g. `runFolders = {'E:\path\to\one_run'}` (or a
single path string). It produces the same three figures, snapshot panels,
Flory fits and number-predicted weight curves, labeled **Single run**. With
only one run, replicate SDs and whole-replicate bootstrap confidence intervals
are undefined (`NaN` in saved results), so no SD shading is drawn even when
`showSpread=true`. Selecting just one completed job from a batch also works.

`averageOptions.usagePercent = [25 50 100]` selects evolution points. Completed
chain histograms use common bins across all runs and points. Each replicate is
normalized to probability per bin **before** averaging, giving each run equal
weight rather than giving runs with more chains more weight. Three separate
figures show Kuhn-count number fraction, Kuhn-count weight fraction, and
end-to-end distance; contour length is omitted. Each figure places the selected
usage points left to right (three subplots for three selected points), without
overlaying different usage points. All subplots within a figure share its x/y
limits, including room for displayed SD bands. Limits are computed independently
for each distribution, not forced to match across figures.
`averageOptions.kuhnYLim = [0 .5]` fixes the number-fraction figure's y range
across all snapshots. Use `[]` for automatic limits; other figures remain
independently scaled. Curves outside a manual range are clipped only for display.

For weight fraction, each replicate uses `sum(n in bin)/sum(n over completed
chains)` before the equal-replicate mean and SD are calculated. Grouped bins sum
actual Kuhn counts, not bin centers. This weights chains by their Kuhn count
within a run, while retaining equal weight for each run. A zero total Kuhn count
is undefined and causes an error by default; with `requireAllReplicates=false`,
that run is excluded from the weight distribution and its contributing count
is reported separately.

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

Outputs are written to a unique folder under `replicate_averages`: three PNGs, a MAT
file containing means, SDs, bin edges and per-replicate histograms, and a CSV
selection table. Edit `outputRoot` in the driver to choose another destination.
Programmatic use: `[fig,data] = averageNetworkDumps(runFolders,averageOptions)`.
`fig` contains three figure handles in the order given by `data.distributionNames`:
`kuhn_number`, `kuhn_weight`, `end_to_end`. The cells in `edges`, `mean`, `sd`, and
`perReplicate` use that order (the second entry is now weight fraction, not contour
length). `distributionReplicateCount` records contributing runs per distribution
and target. `outputFile` is a base filename; the three exports append those
distribution names, and their paths are returned in `data.outputFiles`.

### Flory number-fraction fits

The main number-fraction plot has one model switch in `average_pkg_cleaned.m`:

```matlab
averageOptions.floryModel = 'shifted'; % 'shifted' or 'original'
averageOptions.floryMinCount = 4;
averageOptions.floryOriginalAnchorWeight = 1; % original only; 10 restores soft anchor
```

Both choices give smooth dashed curves and properly normalized discrete
probabilities. The driver produces the usual three horizontal-panel figures;
the selected model appears on the number-fraction figure. `fitFlory` controls
number fitting; `fitFloryWeight` controls weight fitting (see below). Distance
panels are not fitted. Set both switches false to hide fitting.

- **shifted** (driver default): `P(n)=(1-p)*p^(n-n_min)` for integer
  `n>=n_min`, with `n_min=floryMinCount`. Probability is zero below
  that minimum. This is a minimum-length-constrained geometric (Flory-type)
  distribution, not a backward extension of a conditional fit. A minimum
  exceeding an observed count raises an error rather than discarding chains
  or silently rescaling the curve. Choose the minimum appropriate to the run.
- **original**: `P(n)=(1-p)*p^(n-1)` for integer `n>=1`, zero at n=0.
  Here `floryMinCount` sets only the fitting cutoff (at least 1); it does
  not shift the model. The full curve begins at n=1, even when fitting n>=4.
  See Flory (1936), [Equation (4), p. 1879](https://pubs.acs.org/doi/10.1021/ja01301a016).

For the shifted model, equal-replicate raw-count log likelihood gives
`p=(mu-n_min)/(1+mu-n_min)`, where mu is the mean of the replicate mean
Kuhn counts. It has no free amplitude or forced peak anchor. With minimum 4
the estimates are unchanged from the earlier conditional-tail fit when every
chain has n>=4. The invalid upward extrapolation below four has been removed.
For p=0 the model is a point mass at n_min.

The original model minimizes squared error against the equal-replicate mean
**absolute** number fractions over n>=the fitting cutoff. Each replicate is
normalized over its entire observed population, not just the fitted tail.
`floryOriginalAnchorWeight=1` gives ordinary least squares. Higher values
increase the squared-error weight of the bin at the fitting cutoff; this is
a soft anchor. An exact n=4 anchor above `27/256 = 0.10546875` is impossible
for the original formula. The optimizer searches both parameter branches and
compares boundaries. The unobserved infinite tail contributes analytically.
Unit-count bins are used internally, so fitting is independent of display
bin width or cropping. Original-fit details, including anchor feasibility
and unweighted/weighted errors, are in each fit's `originalDiagnostics`.

Neither mode extrapolates outside its support or normalizes over the displayed
range. Saved bin probabilities sum the discrete model within each bin, retaining
the probability outside the plot separately. With unit bins the curve directly
evaluates the model; with grouped bins a continuous envelope passes through
the grouped probabilities. A grouped boundary bin can straddle the minimum;
its center represents the entire bin, not a forbidden single chain length.

`floryBootstrapSamples=2000` and `floryBootstrapSeed=1729` give reproducible
95% percentile intervals by resampling whole replicates without changing the
global random stream. Intervals measure replicate uncertainty, not model
adequacy. Fewer than two fitting replicates, or zero bootstrap draws, gives
no interval. Original-mode replicates without any count at/above the cutoff
cause an error by default; `requireAllReplicates=false` excludes them from
fitting and reports the smaller cohort.

The figure legend identifies the model; the console, `flory_fits.csv`, and
MAT file report p and its interval. The CSV also records `model`,
`modelMinCount` (support), `fitMinCount` (fitting cutoff), and
`originalAnchorWeight` (NaN when inapplicable). `fullCdfDistance` and
`fullBinRMSE` now compare the **selected normalized model** with the full
empirical distribution. `tailCdfDistance` compares conditional tail shapes;
it is not the original mode's absolute-fraction least-squares objective.
These are discrepancy measures, not hypothesis-test p-values. Model
probabilities, normalization, and diagnostics now consistently describe the
chosen mode rather than a separate zero-based reference distribution.

Fitted p is an effective continuation/shape parameter, not automatically
crosslinker usage or chemical conversion. The simulation's length-dependent
joining rates and changing concentrations need not satisfy Flory's assumptions.
[NIST's geometric model](https://www.itl.nist.gov/div898/software/dataplot/refman2/ch8/geopdf.pdf)
uses stopping probability 1-p in this notation.

The earlier separate comparison utility remains available for saved averages:
```matlab
loaded = load('path/to/replicate_average.mat');
options = struct('minCount',4,'anchorCount',4,'anchorWeight',10);
[fig,comparison] = compareOriginalFloryFits(loaded.averageData,options);
```
It overlays the original fit and, when present, the saved shifted fit.
It is no longer run automatically by the averaging driver.

### Flory weight fractions and consistency of p

Weight fraction weights the number distribution by **actual Kuhn count n**:
`W(n)=n*P(n)/E[n]`. For a shifted geometric number model with minimum c,
`E[n]=c+p/(1-p)`, hence

```text
W(n) = n*(1-p)^2*p^(n-c) / (c*(1-p)+p),   n=c,c+1,...
c=4: W(n) = n*(1-p)^2*p^(n-4) / (4-3*p)
c=1: W(n) = n*(1-p)^2*p^(n-1)             (original Flory)
```

These are normalized over the full integer support. The prefactor is n, not
n-c: all segments, including the compulsory minimum, contribute to Kuhn weight.
The original result is Flory's Equation (2), p. 1878 of the paper linked above.
The shifted expression follows by size-weighting the shifted number model;
simply shifting the original weight formula without its new denominator is wrong.

The driver enables:

```matlab
averageOptions.fitFlory = true;              % fit number fractions and supply p
averageOptions.fitFloryWeight = false;        % skip independent weight fitting
averageOptions.showWeightFromNumberFit = true; % predict weight with number-fit p
averageOptions.weightYLim = [];              % own shared y range; [] = automatic
```

Both distributions use the existing `floryModel` and `floryMinCount` choices.
With these settings the weight plot shows a black dashed prediction using the
number-fit p, labeled p_N per snapshot. No weight fitting is performed. The
prediction has no additional fitted parameter and works independently of
`fitFloryWeight`; it requires number fitting (`fitFlory=true`). Set
`showWeightFromNumberFit=false` to hide it. If `fitFloryWeight=true`, the plot
adds a black dashed independent weight fit and shows the number-based prediction
in gray dash-dot, with p_W and p_N labeled per snapshot. Weight fitting can be
used alone. In programmatic calls, an omitted
`fitFloryWeight` follows `fitFlory`; the driver sets both explicitly.

For an exact common geometric population the number and weight laws share p.
Separately fitted p values need not coincide in finite data: number likelihood
and weight least squares use different objectives, and weight fractions emphasize
long chains. Each replicate is weight-normalized separately before averaging;
reweighting the averaged number histogram is generally a different operation.
The estimates use the same observations and are not statistically independent.

The weight fit minimizes the squared difference from the equal-replicate mean
weight PMF, on individual integer counts at/above the fitting cutoff, including
the infinite zero-observation tail analytically. Fractions retain the FULL
per-replicate Kuhn total as denominator, even in original mode when fitting only
n>=4. No peak anchor or extra amplitude is fitted. Cropping and display binning
cannot change p. A shifted minimum exceeding observed counts raises an error.
Zero-Kuhn-mass populations are excluded only when `requireAllReplicates=false`,
as for the empirical weight plots. Missing original-model tails are also checked.

`flory_weight_fits.csv` and the MAT fields `floryWeightFitTable` / `floryWeightFit`
save the weight estimate and its 95% whole-replicate bootstrap interval, number
estimate, p_W-p_N, and a paired bootstrap interval for that difference. Pairing
is used only with the same fitting cohort, bootstrap seed and number of draws.
Otherwise the difference interval is NaN. With five replicates these intervals
are approximate, not proof of model validity or a calibrated hypothesis test.
The table also records full CDF/RMSE discrepancies, the weight-fit tail SSE,
the number-prediction tail SSE, fit counts, support and normalization details.

At p=0 with c>=1 the weight model is a point mass at c. If a zero minimum is
chosen, zero-length chains carry no weight; for p>0 its weight law equals the
c=1 law. At c=0,p=0 the number population has zero Kuhn mass, so the fitter
explicitly labels and uses the p->0 weight-distribution limit at n=1.

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
