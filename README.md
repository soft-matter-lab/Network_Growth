# Network Growth

MATLAB stochastic network growth with periodic boundaries: growing chains
(Case A) or end-linking preformed chains (Case B). Add the single shared `SRC` to
the MATLAB path, or use the drivers which do this automatically.

## Driver folders and shared source

There is **one package and one `SRC` folder**. The two driver folders contain
editable settings and entry scripts only:

| Task | `drivers/growing_chains/` | `drivers/preformed_chains/` |
| --- | --- | --- |
| Single run | `run_single_growing.m` | `run_single_preformed.m` |
| Sweep / replicates | `run_sweep_growing.m` | `run_sweep_preformed.m` |
| Saved-state visualizer / video / prestretch | `visualize_growing.m` | `visualize_preformed.m` |
| Replicate distributions | `average_growing.m` | `average_preformed.m` |
| Offline LAMMPS builder | `build_lammps_growing.m` | `build_lammps_preformed.m` |

All ten drivers find `SRC` relative to their own location. New driver output
paths are anchored to the project root, independent of MATLAB's current folder.
Case A uses `runs_pkg_cleaned` / `batch_runs_pkg_cleaned_sweeps`; Case B uses
`runs_preformed` / `batch_runs_preformed`. Each sample has its own folder.
The five legacy root scripts are archived unchanged under
`tmp/legacy_root_drivers_20261109/` for reference to their old settings.
Use the drivers above; the archived scripts retain their old root-relative paths.

`runNetwork_cleaned` shares scheduling, mechanics, bookkeeping and output.
`buildGrowingChainEvents` and `buildPreformedChainEvents` implement separate
chemistries. `runNetworkA_cleaned` and `runNetworkB_cleaned` are route wrappers;
`runNetworkSweep` dispatches from `cfg.growthMode`.

## Preformed chains (Case B)

```matlab
cfg = defaultNetworkBConfig();
cfg.chainLengths = [20 60 120];
cfg.chainNumberFractions = [.5 .3 .2]; % fractions of chains, not Kuhn segments
cfg.chainCounts = [];                % optional exact count per species
cfg.fA = 4;
cfg.kAddition = 1;                   % area/time, includes both chain ends
cfg.kChem = 1;                       % fresh-crosslinker rate constant
cfg.kEndEnd = 1;                     % existing-junction multiplier on kChem
cfg.allowSelfLoops = false;
cfg.kSelfLoop = .01;
cfg.selfLoopDistanceFactor = 1;
```

Supply one to three distinct positive integer lengths. Number fractions must
be nonnegative and sum to one. The meaning of `xA` is unchanged:
`NA = round(xA*Ntotal)` and the Kuhn budget is `Ntotal-NA`. With number fractions
`f`, ideal counts are `budget*f/sum(f.*chainLengths)`. Counts are rounded down,
then at most one chain per species is added in descending fractional-remainder
order when it fits the remaining budget. This does not guarantee exact budget
filling or exact requested fractions. Both actual counts/fractions and unused
budget are reported. The unused rounding budget is not a free-monomer reservoir
and does not enter mechanics or LAMMPS material counts.

If `chainCounts` is supplied, it overrides number fractions and must fit the
Kuhn budget. `resolvePreformedFeed(cfg)` previews counts and the stoichiometric
ratio `2*sum(chainCounts)/(fA*NA)`, including the seed in NA. Ratio one means
balanced initial functional groups, not guaranteed complete conversion.

Only three event types occur:

1. A bare arm captures a whole chain of species s with rate
   `kAddition*freeChainCounts(s)/L^2`. Summing the q bare arms at a junction gives
   `kAddition*q*C_P,s`. There is no extra factor of two. The chain becomes a
   dangling chain of length N_s; its reservoir count decreases by one.
2. A dangling chain captures a fresh crosslinker with rate `fA*kChem*nA/L^2`.
   Placement samples the same 2D Gaussian cloud as the existing-node kernel:
   `r = b*sqrt(N_s)*sqrt(-log(U))`, uniform polar angle. Thus RMS radius is
   `b*sqrt(N_s)` (not the mean radius, nor a fixed sampled bond distance).
3. A dangling chain terminates on a bare arm of an existing junction with rate
   `kEndEnd*kChem/(pi*N_s*b^2)*exp(-r^2/(N_s*b^2))` per arm, using the periodic
   minimum-image anchor distance. A self-loop uses `kSelfLoop` in place of
   `kEndEnd` and `r_eff = selfLoopDistanceFactor*b*sqrt(N_s)` in this rate only.

No segment growth, two-dangling joining, or independent nucleation occurs.
Case A's minimum-growth thresholds do not apply to these already formed
chains. Every dangling/completed chain retains one prescribed length. Geometry
can relax while each species' cloud width remains constant. A self-loop needs
a spare bare arm on its own anchor. Both routes permit at most one self-loop
per junction and one strand per distinct junction pair. Junction functionality
counts completed bonds (self-loops twice), dangling attachments and bare arms.

`NtotalScale` increases computed feed. Fill-mode repulsion continues to use
nominal unscaled crosslinker count; an explicit `Ntotal` still ignores scale.
Feed-mode repulsion uses the actual allocated Kuhn count. LAMMPS export uses
incorporated Kuhn segments (completed, self-loop, dangling), divided by formed
crosslinkers. Free chains and unused budget are excluded. Only distinct-node
completed strands become LAMMPS physical bonds, with their own N_s.

Case B saves the existing geometry arrays plus `preformed` reservoir and
conversion records. In `results.nB` and history column 6, B denotes free Kuhn
segments in whole supplied chains. Version-2 dumps explicitly store
`freeKuhnSegments`, `freeMonomers=0`, and per-species free/dangling/completed
counts; metadata records the event and history definitions. `chainCountHistory`
contains the initial reservoir followed by one row per reaction. Old Case A
dumps remain readable. All default usage milestones still measure consumption
of initially free crosslinkers, excluding the seed; they do not mean chain
capture or chain completion. Event-interval dumps also capture subsequent
changes after crosslinker depletion.

The Case B visualizer adds `preformed_progress.png/.csv`: remaining chains and
captured/completed conversion by species. Captured means one or both ends are
attached; completed includes self-loops. Number/weight distributions and
prestretch use the shared plotting functions. Flory fits default off for Case B.

Sweeps support the same seed/composition modes plus explicit parameter sets:

```matlab
sweep.parameterSets = {
    struct('chainLengths',20,'chainNumberFractions',1,'chainCounts',[],'fA',4)
    struct('chainLengths',[20 60],'chainNumberFractions',[.7 .3],'chainCounts',[],'fA',4)
};
```

Each set is crossed with the selected compositions and seeds. An empty cell
array uses `baseCfg` alone. Plans record per-species lengths/counts, budget
remainder and end/arm ratio; summaries record remaining species and conversion.
Select one `parameterSetIndex` and xA when averaging replicates; differing
physical settings are rejected by the shared averaging function.

## Common sweep options

Use `drivers/growing_chains/run_sweep_growing.m` or
`drivers/preformed_chains/run_sweep_preformed.m` for the selected synthesis route.

Edit a driver's settings before running it. Preview the generated plan to see
the exact number of runs and resolved feed for your selected settings.

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

## Self-loop kinetics (Case A)

The growing-chain single and sweep drivers expose `allowSelfLoops`, `kSelfLoop`, and
`selfLoopDistanceFactor` (`cfg` for single runs, `baseCfg` for sweeps):

```matlab
cfg.allowSelfLoops = true;       % false disables both mechanisms
cfg.kSelfLoop = 0.01;            % independent multiplier on kChem
cfg.selfLoopDistanceFactor = 1;  % alpha in d_eff = alpha*b*sqrt(n1+n2)
```

The library and drivers keep `allowSelfLoops = false` until enabled.
Two same-anchor active ends can close through either mechanism:

| Mechanism | Reactants | Completed loop Kuhn count |
| --- | --- | --- |
| Dangling chain returns to its anchor | One growing end with n>0 plus a spare bare arm with n=0 | n |
| Two dangling chains join on their common anchor | Two growing ends with n1>0 and n2>0 | n1+n2 |

A return to an anchor requires a **spare reactive arm**; it cannot react with
an already occupied attachment site. Each event consumes two active arms,
adds one `[p,p,n1+n2,0,0]` bond, and adds two to the anchor's degree. It consumes
no free crosslinker or background monomer. The existing eligibility threshold
is `n1+n2 >= 2*minKuhnEndEnd`, with positive total length. Bare-to-bare closure
is excluded. The existing no-duplicate-bond rule permits at most one self-loop
per anchor. Distinct eligible arm pairs each have their own propensity and
are enumerated once; there is no extra functionality multiplier.

For a self-loop the reaction kernel uses a **rate-only effective distance**:

```text
nTotal = n1+n2
d_eff = selfLoopDistanceFactor * b * sqrt(nTotal)
w_loop = kSelfLoop * kChem / (pi*nTotal*b^2)
         * exp(-d_eff^2/(nTotal*b^2))
       = kSelfLoop * kChem / (pi*nTotal*b^2)
         * exp(-selfLoopDistanceFactor^2)
```

For one dangling chain this gives `d_eff = b*sqrt(n)` at the default factor 1;
for two dangling chains it gives `b*sqrt(n1+n2)`. This is a deliberate kinetic
suppression, not a measured tip separation or a new physical bond length.
The completed self-loop still has coincident endpoints and zero end-to-end
length. With factor 1, the added suppression is `exp(-1)` for both mechanisms;
it introduces no further chain-length dependence beyond the existing `1/nTotal`.

`kSelfLoop` replaces `kEndEnd` for same-anchor events, so changing `kEndEnd`
affects joining between different anchors only. Setting `kSelfLoop = 0`
disables loops. To reproduce the old zero-distance loop rate, enable loops,
set `selfLoopDistanceFactor = 0`, and set `kSelfLoop` equal to the old `kEndEnd`.
Both mechanisms remain history event type 3. Completed loops still count
toward chain-bound material and the exported density radius, and are omitted
from physical LAMMPS bonds.

## LAMMPS export

The single-run and sweep drivers in both mode folders enable
`saveLammps = true`. After final relaxation, each sample folder receives:

- `network.data`: LAMMPS data file with `Atoms # bpm/sphere` and `Bonds`.
- `bond.table`: NetworkGen's `KEY` section, with rows
  `bondID atom1 atom2 N_Kuhn b`, for your custom `bpm/poly` bond style.
- `network.localdensity.table`: one type-1 local-density potential.

This is a final network export, independent of MAT progress dumps, `saveData`,
and `sweep.saveResults`. The library default is off; enable `cfg.saveLammps`
and set `cfg.outputDir` for programmatic runs. Existing export files are never
overwritten. Paths, counts, omitted loop rows, and potential parameters are
returned in `results.lammps` and included in the saved final MAT data when enabled.
Sweep summaries include `lammpsStatus` and `lammpsRepulsiveRadius`.

The radius uses the **actual formed crosslinker count**, including the seed and
any isolated nuclei, rather than the feed count:

```text
N_kuhn_tot = sum(bonds(:,3)) + sum(activeEnds(:,4))
N_formed  = size(nodes,1)
A_node    = (N_kuhn_tot / N_formed) * pi*b^2/4
r_repulsive = sqrt(A_node/pi) = (b/2)*sqrt(N_kuhn_tot/N_formed)
```

By default, `lammpsIncludeDanglingKuhn = true`: count all Kuhn segments in
completed inter-junction chains, self-loops, and dangling chains. Unused
background Kuhn segments and unused background crosslinkers are excluded from
the numerator and denominator, respectively. Bare arms contain zero segments.
Dangling tips have no explicit coordinates in this growth model and are not
exported as extra atoms or bonds. Setting `lammpsIncludeDanglingKuhn = false`
restores the reference generator's completed-chain-only counting convention.

Following the reference potential, `R1 = 0.8*r_repulsive`,
`R2 = 4*r_repulsive`, `rho0 = 0.8*(R2/r_repulsive)^2 = 12.8`, and
`F(rho) = lammpsKLD*(rho-rho0)^2`. The associated pair repulsion cutoff is
`rc = 2*r_repulsive`, available as `results.lammps.pairRepulsionCutoff` and
in the data-file comment. The atom diameter/density remain 1/1, matching
NetworkGen; the effective repulsion radius is carried by the potential cutoffs.
This exported potential does not change the growth simulation's `rEq`,
pairwise relaxation, or chain model.

| Setting | Default |
| --- | --- |
| `lammpsKLD` | `0.414` (NetworkGen class default; independently adjustable) |
| `lammpsNRho` | `100000` density samples |
| `lammpsRhoRange` | `[0 500]` |
| `lammpsZThickness` | `10`, with atoms at z=0 |
| `lammpsIncludeDanglingKuhn` | `true` (completed chains, loops, and dangling chains) |
| `lammpsPeriodicPolicy` | `'warn'`; use `'error'` to reject image changes |

The density table has exactly `lammpsNRho` samples with spacing
`(rho_max-rho_min)/(lammpsNRho-1)`, consistent with the
[LAMMPS local-density reader](https://github.com/lammps/lammps/blob/stable/src/MANYBODY/pair_local_density.cpp).
This corrects the reference generator's extra sample and inconsistent spacing.
The bond table is a per-chain parameter table for your custom style, not the
energy/force grid used by LAMMPS `bond_style table`.

Use `dimension 2`, `atom_style bpm/sphere`, and periodic x/y boundaries in your
LAMMPS input. Coordinates are wrapped into `[0,L)`, preserving node IDs.
For your existing hybrid pair setup, point `read_data` at `network.data`,
the `bpm/poly` bond coefficients at `bond.table KEY`, and
`pair_coeff * * local/density` at `network.localdensity.table`.
Choose the other force, damping, and fracture coefficients in your LAMMPS input.

Self-loops are omitted from both exported bond lists because LAMMPS rejects
atom-to-itself bonds; their segments still contribute to the radius. Original
bond-row mappings and omitted rows are recorded in `results.lammps`. They remain
unchanged in the growth results and snapshots. A state with no inter-junction
bonds is reported as `skipped_no_interjunction_bonds`, with a warning and no
three-file export, since the custom bond reader rejects an empty table.

LAMMPS `bpm/poly` uses minimum-image bond vectors. If a growth strand's stored
periodic image selects a different vector, export emits a warning, records its
original row in `nonMinimumImageBondRows`, and reports
`written_with_image_changes`. Such a state will start with different strand
extensions in LAMMPS. Use `lammpsPeriodicPolicy = 'error'` to prevent that export;
a standard LAMMPS data file cannot encode independent per-bond image offsets.

### Build ingredients from an older saved simulation

Run **`build_lammps_growing.m`** or **`build_lammps_preformed.m`** in the
corresponding driver folder. Set `inputFolder`
to one simulation folder (or its `dumps` subfolder), or leave it empty for a
folder picker. The driver writes the three LAMMPS files directly into that
simulation folder without rerunning growth or changing any saved input data.

Automatic source selection prefers `results_struct.mat`, then the raw final
`Network_PartA_v2_Cleaned_Data.mat` (or a custom-named MAT containing `nodes`,
`bonds`, `activeEnds`, and `cfg`), then a final network snapshot. It reads the
original box size `L`, Kuhn length `b`, coordinates, Kuhn counts, and periodic
images. Older files do not need any of the new LAMMPS configuration fields.
The builder always includes connected, loop, and dangling segments in the
radius calculation and uses only formed crosslinkers, excluding both free
reservoirs. Saved output paths are ignored in favor of the selected folder.

The additional `lammps_build_report.mat` records source filename/phase, counts,
radius/cutoffs, output paths, and export settings. Existing LAMMPS files or a
build report are never overwritten; set `builderOptions.outputDir` to a fresh
folder to make another export. `builderOptions.sourceFile` can select a
specific MAT file (absolute path or relative to the simulation folder).
If no final state exists, the builder stops. Set
`builderOptions.allowLatestSnapshot = true` explicitly to use the latest
available network checkpoint; it warns and records that this state is non-final.

Distribution-only dumps, histogram CSVs, or an atom-only text dump cannot
reconstruct the complete export. The builder requires the saved node positions,
bond endpoints/Kuhn counts/periodic images, dangling-chain data, and `L`/`b`.
It reports missing data instead of inventing geometry or treating all feed
segments as chain-bound material.

Programmatic use:

```matlab
addpath('SRC')
info = buildLammpsIngredients(runDir);
% Optional export settings (never overrides saved L or b):
options.exportConfig = struct('lammpsKLD',0.414,'lammpsNRho',100000);
options.outputDir = fullfile(runDir,'lammps_rebuild');
info = buildLammpsIngredients(runDir,options);
```

The lower-level exporter also accepts a previously loaded final state:

```matlab
addpath('SRC')
s = load(fullfile(runDir,'results_struct.mat'),'results');
r = s.results;
cfg = r.config;
cfg.outputDir = runDir; % or a fresh folder if already exported
info = writeLammpsNetwork(r.nodes,r.bonds,r.activeEnds,cfg);
```

Run the export and saved-state builder regression tests with
`runtests({'tests/testLammpsExport.m','tests/testLammpsIngredientsBuilder.m'})`.

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

**Start with `visualize_growing.m` or `visualize_preformed.m` in the corresponding driver folder.** Set
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
  weight-fraction histograms, 2D displacement, and bond prestretch. Several crossed milestones can share
  one state/PNG. `checkpoint_index.csv` maps the files to requested/actual usage.
- `saveFinalPlot`: a separate `final_state.png` for the final relaxed network.
- `saveComparison`: `distribution_comparison.png`, overlaying Kuhn-count,
  physical contour-length, completed-chain end-to-end, and prestretch distributions using
  shared bins. `comparisonUsagePercent` defaults to `[0 25 50 75 100]`; set it to
  `0:5:100` for all checkpoints. `population` chooses the contour population and
  `showComparisonFigure` controls whether the comparison stays open.

`vizOptions.showDanglingClouds = true` overlays Gaussian clouds for growing
dangling chains on the network panels in both the replay movie and checkpoint /
final PNGs. Set it to `false` for bonds and junctions only. Clouds have dark-blue
centers fading outward through light blue to white; `danglingCloudOpacity`
(default 0.85, range 0 to 1) controls their visibility and `heatGrid` controls
resolution. Clouds use each dangling chain's parent position and Kuhn count,
with the existing periodic envelope width proportional to `b*sqrt(n)`.
Bare arms (`n=0`) and completed chains, including self-loops, add no clouds here.
These are relative visualization envelopes (overlaps use the maximum), not a
normalized density or a pair RDF g(r). The simulation's reaction rates and
live rendering are unchanged.

Checkpoint and final-state figures contain five panels: network, chain-number
fraction, and chain weight fraction across the top; 2D end-to-end displacement
at lower left; and bond prestretch across the remaining lower panels.
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

The prestretch distribution uses the dimensionless quantity
`lambda_i = r_i/(n_i*b)`, the end-to-end length divided by the contour length
of each completed bond. This is computed directly from existing saved data;
no simulation rerun or new dump fields are needed. Magnitude-only distribution
dumps use `completed.endToEnd` and `completed.kuhnCount` with the saved `b`
(or saved `contourLength` if `b` is absent). Network-only snapshots reconstruct
lengths using `nodes`, `bonds`, and the stored periodic images. The bond vectors
are not changed to a new minimum image.

Self-loops contribute at lambda=0; dangling chains are excluded because no
explicit tip position/end-to-end length was saved. Invalid or zero-contour
samples are excluded. Values above lambda=1 are retained; the dotted line at
1 marks the contour length. This is the requested contour-normalized extension,
not normalization by the Gaussian coil size `sqrt(n)*b`.

`vizOptions.prestretchBins` (default 30), `prestretchXLim`, and `prestretchYLim`
control the new panel. Empty limits select common ranges across every saved
state, including final relaxation. Probability per bin is normalized by the
number of valid completed bonds; manually cropped ranges retain that denominator.
The comparison figure also includes prestretch and uses the same bins and y
limits in the visualizer. `comparison_index.csv` records its sample count and
availability. `networkDumpPrestretch(snapshot)` returns the raw lambda values.

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

The comparison plots Kuhn counts, physical contour lengths, completed-chain
end-to-end distances, and prestretch. `population` selects `'completed'` (default), `'dangling'`,
or `'all'` for the contour panels; the end-to-end and prestretch panels always use completed
chains. Empty populations are labeled, not replaced with zero-length chains.
`selection` reports files, events, actual conversions and sample counts.
Only saved milestone markers are selected; missing targets produce a warning.
Multiple targets can refer to the same file in small systems. Replicates are
shown separately, without automatic pooling or averaging.

Optional `nEdges`, `contourEdges`, `rEdges`, `prestretchEdges` fix common bins across comparisons;
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

Run `average_growing.m` or `average_preformed.m` in the corresponding driver
folder. Enter the replicate run
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

The growing-chain number-fraction plot has one model switch in `average_growing.m`:

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
