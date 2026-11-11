% Case B single run. Edit this script; all implementation lives in shared SRC.
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(projectRoot,'SRC'));
cfg = defaultNetworkBConfig();

% Feed: xA = crosslinkers / (crosslinkers + Kuhn-segment budget).
cfg.L = 100;
cfg.b = 1.6;
cfg.Ntotal = [];                 % [] = computed from areaFraction, L, b, scale
cfg.NtotalScale = 1.75;
cfg.areaFraction = 1;
cfg.xA = .05;
cfg.fA = 4;
cfg.chainLengths = [10];           % e.g. [20 60] or [20 60 120]; max three
cfg.chainNumberFractions = [1];    % e.g. [.7 .3] or [.5 .3 .2]; sum = 1
cfg.chainCounts = [];            % optional exact integer counts; must fit budget

% Case B reactions. No dangling growth, dangling joining, or new-cluster nucleation.
cfg.kAddition = 1;               % k_p, area/time; includes both supplied-chain ends
cfg.kChem = 1;                   % fresh capture: fA*kChem*C_A per dangling chain
cfg.kEndEnd = 1;                 % existing-junction kernel multiplier on kChem
cfg.allowSelfLoops = false;      % dangling chain to own spare bare arm only
cfg.kSelfLoop = .01;             % independent loop multiplier on kChem
cfg.selfLoopDistanceFactor = 1;  % rate-only d_eff = factor*b*sqrt(N)
% At most one loop per junction and one strand per distinct junction pair.
cfg.randSeed = 11231;

% Shared mechanics. Fill-mode range uses nominal unscaled crosslinker count.
cfg.kT = 4.14;
cfg.chainModel = 'cohen';        % 'cohen' or 'gaussian'
cfg.kRep = .828;
cfg.kRepAnh = 1.656;
cfg.rEqMode = 'fill';
cfg.nodeAreaFraction = 1;
cfg.eqIntegrator = 'damped2nd';
cfg.massEq = 1;
cfg.gammaEq = .5;

cfg.showLivePlot = true;
cfg.showNodeConnectivity = true; % live/video nodes: light gray to black; false = all black
cfg.showSphericalNodes = true;   % shaded beads; false = flat circles
cfg.colorPreformedBonds = true;  % ascending N: black, red, burgundy
cfg.makeVideo = true;
cfg.showFinalHistogram = true;
cfg.saveFinalDistributions = true;
cfg.saveData = true;
cfg.saveLammps = true;
cfg.lammpsKLD = .414;
cfg.lammpsIncludeDanglingKuhn = true;
cfg.saveDumps = true;
cfg.dumpUsagePercent = 0:5:100; % crosslinker consumption, not chain capture
cfg.dumpEveryEvents = 10;      % also capture evolution after fresh A is exhausted
cfg.dumpTimeInterval = 0;
cfg.dumpSaveNetwork = true;
cfg.dumpSaveDistributions = true;
cfg.dumpSaveFinal = true;

feed = resolvePreformedFeed(cfg); disp(feed);
outputRoot = fullfile(projectRoot,'runs_preformed');
runLabel = sprintf('preformed_xA_%.5g_seed_%d',cfg.xA,cfg.randSeed);
[runDir,runName] = createUniqueRunFolder(outputRoot,runLabel);
cfg.outputDir = runDir;
results = runNetworkB_cleaned(cfg);
save(fullfile(runDir,'results_struct.mat'),'results');
fprintf('Saved run %s to %s\n',runName,runDir);
