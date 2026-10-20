% Batch driver for multiple cleaned Network A runs.
% Each run gets its own unique folder; no files are overwritten.

addpath(fullfile(fileparts(mfilename('fullpath')), 'SRC'));

baseCfg = defaultNetworkAConfig();

% ---------------- Batch-level defaults ----------------
baseCfg.L = 100;
baseCfg.areaFraction = 0.8;         % Ntotal is computed from L, b, areaFraction
baseCfg.fA = 3;
baseCfg.fB = 2;
baseCfg.kRep = 1;
baseCfg.eqIntegrator = 'damped2nd';
baseCfg.massEq = 1;
baseCfg.gammaEq = 0.5;
baseCfg.showLivePlot = false;        % typically off for batch
baseCfg.makeVideo = true;           % turn on only if needed
baseCfg.showFinalHistogram = false;
baseCfg.saveData = true;
baseCfg.videoName = 'networkA_formation_v2_cleaned';
baseCfg.dataFile = 'Network_PartA_v2_Cleaned_Data.mat';
% ------------------------------------------------------

% Only composition and functionality vary between runs.
runOverrides = {
    struct('xA', 0.05, 'fA', 3, 'fB', 2), ...
    struct('xA', 0.10, 'fA', 3, 'fB', 2), ...  
    struct('xA', 0.15, 'fA', 3, 'fB', 2), ...
    struct('xA', 0.20, 'fA', 3, 'fB', 2)
};

batchRoot = fullfile(pwd, 'batch_runs_pkg_cleaned');
[batchDir, batchName] = createUniqueRunFolder(batchRoot, 'batch');
fprintf('Batch folder: %s\n', batchDir);

nRuns = numel(runOverrides);
summary = table('Size', [nRuns, 16], ...
    'VariableTypes', {'double','string','double','double','double','double','double', ...
        'double','double','string','double','double','double','double','double','double'}, ...
    'VariableNames', {'runIndex','runName','xA','L','Ntotal','fA','fB', ...
        'gammaEq','massEq','eqIntegrator','kRep','nNodes','nBonds', ...
        'nActiveEnds','ssaTime','runtimeSeconds'});

allResults = cell(nRuns,1);

for k = 1:nRuns
    cfg = mergeConfig(baseCfg, runOverrides{k});
    cfg.Ntotal = resolveNtotal(cfg);

    runLabel = sprintf('run_%03d_xA_%s_N_%d_fA_%d_fB_%d', ...
        k, strrep(sprintf('%.4f', cfg.xA), '.', 'p'), cfg.Ntotal, cfg.fA, cfg.fB);
    [runDir, runName] = createUniqueRunFolder(batchDir, runLabel);
    cfg.outputDir = runDir;

    fprintf('\n=== BATCH RUN %d/%d: %s ===\n', k, nRuns, runName);
    results = runNetworkA_cleaned(cfg);
    allResults{k} = results;

    save(fullfile(runDir, 'results_struct.mat'), 'results');

    summary.runIndex(k) = k;
    summary.runName(k) = string(runName);
    summary.xA(k) = cfg.xA;
    summary.L(k) = cfg.L;
    summary.Ntotal(k) = cfg.Ntotal;
    summary.fA(k) = cfg.fA;
    summary.fB(k) = cfg.fB;
    summary.gammaEq(k) = cfg.gammaEq;
    summary.massEq(k) = cfg.massEq;
    summary.eqIntegrator(k) = string(cfg.eqIntegrator);
    summary.kRep(k) = cfg.kRep;
    summary.nNodes(k) = size(results.nodes,1);
    summary.nBonds(k) = size(results.bonds,1);
    summary.nActiveEnds(k) = size(results.activeEnds,1);
    summary.ssaTime(k) = results.t;
    summary.runtimeSeconds(k) = results.runtimeSeconds;
end

writetable(summary, fullfile(batchDir, 'batch_summary.csv'));
save(fullfile(batchDir, 'batch_summary.mat'), 'summary', 'allResults', 'runOverrides', 'baseCfg');

fprintf('\nBatch complete: %s\n', batchName);
disp(summary);
