% New batch driver: edit base settings and sweep mode, then run this script.
addpath(fullfile(fileparts(mfilename('fullpath')),'SRC'));
baseCfg = defaultNetworkAConfig();

% Shared physics (fixed across all runs except xA and optional NtotalScale).
baseCfg.L = 100;
baseCfg.Ntotal = [];             % computed feed; preserves nominal fill-mode rEq
baseCfg.NtotalScale = 1.35;
baseCfg.areaFraction = 1;
baseCfg.xA = .05;                % used in 'replicates' mode
baseCfg.fA = 4;
baseCfg.fB = 2;                  % metadata; the current model uses B2
baseCfg.b = 1.6;
baseCfg.kT = 4.14;
baseCfg.chainModel = 'gaussian';
baseCfg.kEndEnd = 1;
baseCfg.minKuhnCrosslink = 4;
baseCfg.minKuhnEndEnd = 2;        % resulting joined chain has n >= 4
baseCfg.nucleation = false;
baseCfg.kNuc = .001;
baseCfg.nucAvoidOccupied = true;
baseCfg.kRep = .828;
baseCfg.kRepAnh = 1.656;
baseCfg.rEqMode = 'fill';
baseCfg.nodeAreaFraction = 1;
baseCfg.eqIntegrator = 'damped2nd';
baseCfg.massEq = 1;
baseCfg.gammaEq = .5;

% Batch output: capture data now, plot it later.
baseCfg.showLivePlot = false;
baseCfg.makeVideo = false;
baseCfg.showFinalHistogram = false;
baseCfg.saveFinalDistributions = false;
baseCfg.saveData = false;        % full final state is saved by sweep.saveResults
baseCfg.saveDumps = true;
baseCfg.dumpUsagePercent = 0:5:100;
baseCfg.dumpEveryEvents = 0;     % 1 = every SSA reaction, N = every N reactions
baseCfg.dumpTimeInterval = 0;    % positive SSA-time interval, 0 = off
baseCfg.dumpSaveDistributions = true;
baseCfg.dumpSaveNetwork = true;  % positions, bonds/images, active arms; offline RDF
baseCfg.dumpSaveFinal = true;

% Modes:
% 'replicates'             : baseCfg.xA at every listed seed
% 'composition'            : every sweep.xA at sweep.seeds(1)
% 'replicated_composition' : every sweep.xA at every listed seed (15 runs below)
sweep = struct();
sweep.mode = 'replicated_composition';
sweep.xA = [.01 .025 .05 .075 .10]; % mole fractions A/(A+B), not A/B
sweep.seeds = [11231 2015 3215464];          % one distinct seed per replicate block
sweep.NtotalScaleByXA = [];      % [] = same baseCfg.NtotalScale for all compositions
% Optional fA=4 calibration, in the same order as sweep.xA:
% xA (%)       1       2.5     5          7.5        10
% scale guide  1.15-1.2 1.25   1.25-1.3   1.35-1.4   1.5-1.6
% sweep.NtotalScaleByXA = [1.175 1.25 1.275 1.375 1.55];
% Clear NtotalScaleByXA when using 'replicates'.
sweep.continueOnError = false;
sweep.saveResults = true;

plan = makeNetworkSweepPlan(baseCfg,sweep);
disp(plan);
outputRoot = fullfile(pwd,'batch_runs_pkg_cleaned_sweeps');
[summary,batchDir] = runNetworkSweep(baseCfg,sweep,outputRoot);
disp(summary);
