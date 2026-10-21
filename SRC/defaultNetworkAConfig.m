function cfg = defaultNetworkAConfig()
% Default configuration for the cleaned Network A simulation package.

cfg = struct();

% Chemistry / composition
cfg.L = 100;
cfg.Ntotal = [];          % [] = computed from L, b, areaFraction (see resolveNtotal)
cfg.areaFraction = 0.8;   % Ntotal*pi*b^2/4 = areaFraction*L^2
cfg.xA = 0.10;
cfg.fA = 3;
cfg.fB = 2;
cfg.kChem = 1;
cfg.kEndEnd = 1;              % multiplier on kChem for end-end (chain-chain / cloud-cloud) bonding
% Nucleation: a free crosslinker is placed as a new, initially disconnected
% cluster (fA bare arms). Rate = kNuc * nA * nB / Abox  (= kNuc*C_A*C_B*Abox).
cfg.randSeed = [];               % integer for reproducible runs; [] = random each run
cfg.nucleation = false;       % false: single cluster grown from one seed (legacy)
cfg.kNuc = 0.01;              % nucleation rate constant (relative to kChem)
cfg.nucAvoidOccupied = true;  % bias new nuclei toward empty regions
cfg.nucCandidates = 20;       % random sites tried; the one farthest from existing junctions wins
cfg.allowSelfLoops = false;   % false: a crosslinker cannot bond two of its own arms together

% Mechanics
cfg.b = 5;
cfg.kT = 4.14;
% Chain tension model: 'cohen' (Cohen Pade approx. of inverse Langevin, finite
% extensibility at n*b) or 'gaussian' (linear spring, 3kT r/(n b^2))
cfg.chainModel = 'cohen';
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

% Repulsion radius: 'feed' = 0.5*b*sqrt(NB0/NA0) (half RMS chain size)
%                   'fill' = crosslinker discs of diameter rEq cover
%                            nodeAreaFraction of the box (2D plane filling)
cfg.rEqMode = 'feed';
cfg.nodeAreaFraction = 1.0;

% Repulsion options
% Force magnitude = kRep*s + kRepAnh*s^2 with s = rEq - d (0 => purely harmonic)
cfg.kRepAnh = 0;
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
