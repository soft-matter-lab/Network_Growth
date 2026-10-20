function cfg = defaultNetworkAConfig()
% Default configuration for the cleaned Network A simulation package.

cfg = struct();

% Chemistry / composition
cfg.L = 100;
cfg.Ntotal = 1000;
cfg.xA = 0.10;
cfg.fA = 3;
cfg.fB = 2;
cfg.kChem = 1;

% Mechanics
cfg.b = 5;
cfg.kT = 4.14;
cfg.mobility = 1;
cfg.kRep = 1;
cfg.dtEq = 0.2;
cfg.maxMove = 0.1 * cfg.b;
cfg.tolV = 1e-2;
cfg.maxEq = 10000;
cfg.tolVfinal = 5e-3;
cfg.maxEqFinal = 8000;

% Equilibration integrator
% "damped2nd": m*x'' + gamma*x' = F(x) (default, less jittery visuals)
% "overdamped": x' = mobility*F(x) (legacy first-order mode)
cfg.eqIntegrator = 'damped2nd';
cfg.massEq = 1;
cfg.gammaEq = 0.5;

% Repulsion options
cfg.repelBondedPairs = true;      % also repel directly bonded node pairs
cfg.useRepulsionBinning = false;  % MATLAB vectorized all-pairs is faster for current problem sizes
cfg.repulsionFullPairThreshold = 120; % use full all-pairs below this N

% Visualization and output
cfg.showLivePlot = true;
cfg.plotEvery = 1;
cfg.eqPlotEvery = 3;
cfg.heatGrid = 150;
cfg.showFinalHistogram = true;

cfg.makeVideo = true;
cfg.videoName = 'networkA_formation_v2_cleaned';
cfg.videoFPS = 12;

cfg.saveData = true;
cfg.dataFile = 'Network_PartA_v2_Cleaned_Data.mat';
cfg.outputDir = '';  % optional run folder; when non-empty, data/movie save inside this folder
end
