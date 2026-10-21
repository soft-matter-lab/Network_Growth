% Main driver for the cleaned Network A workflow.
% Edit parameters below, then run this file.
% Radhika downloaded this from the original Network_Growth repo and cleaned it up to remove unnecessary dependencies and make it easier to run.

addpath(fullfile(fileparts(mfilename('fullpath')), 'SRC'));

cfg = defaultNetworkAConfig();

% ---------------- User-facing inputs ----------------
cfg.L = 50;            % box side length
cfg.areaFraction = 1.0; % Ntotal is computed: Ntotal*pi*b^2/4 = areaFraction*L^2
cfg.xA = 0.05;         % crosslinker feed mole fraction
cfg.fA = 4;            % crosslinker functionality
cfg.fB = 2;            % linker functionality (kept as metadata)

cfg.b = 1.6;
cfg.kT = 4.14;
cfg.chainModel = 'gaussian';       % 'cohen' (Langevin/Cohen) or 'gaussian'
cfg.kEndEnd = 10;               % end-end (cloud-cloud) rate multiplier on kChem; 1 = same as growth

cfg.randSeed = 99910;               % integer = reproducible; [] = different every run

cfg.nucleation = false;          % stochastic nucleation of new clusters (false = single seed)
cfg.kNuc = 0.005;                % nucleation rate constant: rate = kNuc*nA*nB/Abox
cfg.nucAvoidOccupied = true;    % place nuclei away from already-formed regions

cfg.kRep = 0.828;          % repulsion strength (increase if you want stronger anti-overlap)
cfg.kRepAnh = 1.656;        % anharmonic (quadratic-force) stiffening; 0 = harmonic only
cfg.rEqMode = 'fill';           % 'feed' (chain-size based) or 'fill' (crosslinker discs fill the box)
cfg.nodeAreaFraction = 0.9;     % 'fill': fraction of box covered by crosslinker discs

cfg.eqIntegrator = 'damped2nd'; % 'damped2nd' or 'overdamped'
cfg.massEq = 1;
cfg.gammaEq = 0.5;

cfg.showLivePlot = true;
cfg.makeVideo = true;
cfg.videoName = 'networkA_formation_v2_cleaned';
cfg.dataFile = 'Network_PartA_v2_Cleaned_Data.mat';
% ----------------------------------------------------

outputRoot = fullfile(pwd, 'runs_pkg_cleaned');
cfg.Ntotal = resolveNtotal(cfg);
fprintf('Ntotal = %d (L = %g, b = %g)\n', cfg.Ntotal, cfg.L, cfg.b);
runLabel = sprintf('single_xA_%s_N_%d', strrep(sprintf('%.4f', cfg.xA), '.', 'p'), cfg.Ntotal);
[runDir, runName] = createUniqueRunFolder(outputRoot, runLabel);
cfg.outputDir = runDir;

results = runNetworkA_cleaned(cfg);

save(fullfile(runDir, 'results_struct.mat'), 'results');
fprintf('Saved run to folder: %s\n', runDir);
fprintf('Run ID: %s\n', runName);
