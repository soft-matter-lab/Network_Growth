% Main driver for the cleaned Network A workflow.
% Edit parameters below, then run this file.

cfg = defaultNetworkAConfig();

% ---------------- User-facing inputs ----------------
cfg.Ntotal = 10000;     % total feed count (crosslinkers + B2)
cfg.xA = 0.05;         % crosslinker feed mole fraction
cfg.fA = 3;            % crosslinker functionality
cfg.fB = 2;            % linker functionality (kept as metadata)

cfg.kRep = 1;          % repulsion strength (increase if you want stronger anti-overlap)
cfg.eqIntegrator = 'damped2nd'; % 'damped2nd' or 'overdamped'
cfg.massEq = 1;
cfg.gammaEq = 0.5;

cfg.showLivePlot = true;
cfg.makeVideo = true;
cfg.videoName = 'networkA_formation_v2_cleaned';
cfg.dataFile = 'Network_PartA_v2_Cleaned_Data.mat';
% ----------------------------------------------------

outputRoot = fullfile(pwd, 'runs_pkg_cleaned');
runLabel = sprintf('single_xA_%s_N_%d', strrep(sprintf('%.4f', cfg.xA), '.', 'p'), cfg.Ntotal);
[runDir, runName] = createUniqueRunFolder(outputRoot, runLabel);
cfg.outputDir = runDir;

results = runNetworkA_cleaned(cfg);

save(fullfile(runDir, 'results_struct.mat'), 'results');
fprintf('Saved run to folder: %s\n', runDir);
fprintf('Run ID: %s\n', runName);
