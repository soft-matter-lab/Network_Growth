% Case B sweeps: mixtures, functionality, composition, and independent seeds.
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(projectRoot,'SRC'));
baseCfg = defaultNetworkBConfig();
baseCfg.L = 100;
baseCfg.b = 1.6;
baseCfg.Ntotal = [];
baseCfg.NtotalScale = 1.4;
baseCfg.areaFraction = 1;
baseCfg.xA = .025;
baseCfg.fA = 4;
baseCfg.chainLengths = 20;
baseCfg.chainNumberFractions = 1;
baseCfg.chainCounts = [];        % optional exact counts, constrained by Kuhn budget

baseCfg.kAddition = 1;           % includes both equivalent supplied-chain ends
baseCfg.kChem = 1;
baseCfg.kEndEnd = 1;             % existing-junction multiplier
baseCfg.allowSelfLoops = false;
baseCfg.kSelfLoop = .01;
baseCfg.selfLoopDistanceFactor = 1;
baseCfg.kT = 4.14;
baseCfg.chainModel = 'cohen';
baseCfg.kRep = .828;
baseCfg.kRepAnh = 1.656;
baseCfg.rEqMode = 'fill';
baseCfg.nodeAreaFraction = 1;
baseCfg.eqIntegrator = 'damped2nd';
baseCfg.massEq = 1;
baseCfg.gammaEq = .5;

baseCfg.showLivePlot = false;
baseCfg.showNodeConnectivity = true; % used if live plotting/video is enabled
baseCfg.showSphericalNodes = true;   % shaded beads; false = flat circles
baseCfg.colorPreformedBonds = true;  % ascending N: black, red, burgundy
baseCfg.makeVideo = false;
baseCfg.showFinalHistogram = false;
baseCfg.saveFinalDistributions = false;
baseCfg.saveData = false;        % results_struct.mat saved by the sweep
baseCfg.saveLammps = true;
baseCfg.lammpsKLD = .414;
baseCfg.lammpsIncludeDanglingKuhn = true;
baseCfg.saveDumps = true;
baseCfg.dumpUsagePercent = 0:5:100;
baseCfg.dumpEveryEvents = 10;
baseCfg.dumpTimeInterval = 0;
baseCfg.dumpSaveNetwork = true;
baseCfg.dumpSaveDistributions = true;
baseCfg.dumpSaveFinal = true;

sweep = struct();
sweep.mode = 'replicates';      % 'composition' or 'replicated_composition' also supported
sweep.xA = [.01 .025 .05];
sweep.seeds = [11231 2015 3215464];
sweep.NtotalScaleByXA = [];     % one per xA; only in composition modes
sweep.parameterSets = {};      % empty = baseCfg mixture/functionality/rates
% To compare mixtures/functionality, uncomment and edit:
% sweep.parameterSets = {
%     struct('chainLengths',20,'chainNumberFractions',1,'chainCounts',[],'fA',4)
%     struct('chainLengths',[20 60],'chainNumberFractions',[.7 .3],'chainCounts',[],'fA',4)
%     struct('chainLengths',[20 60 120],'chainNumberFractions',[.5 .3 .2],'chainCounts',[],'fA',3)
% };
% Each parameter set is crossed with the selected compositions and seeds.
sweep.continueOnError = false;
sweep.saveResults = true;
plan = makeNetworkSweepPlan(baseCfg,sweep); disp(plan);
outputRoot = fullfile(projectRoot,'batch_runs_preformed');
[summary,batchDir] = runNetworkSweep(baseCfg,sweep,outputRoot);
disp(summary);
