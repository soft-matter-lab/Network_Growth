% Batch driver for multiple cleaned Network A runs.
% Each run gets its own unique folder; no files are overwritten.

baseCfg = defaultNetworkAConfig();

% ---------------- Batch-level defaults ----------------
baseCfg.showLivePlot = false;        % typically off for batch
baseCfg.makeVideo = false;           % turn on only if needed
baseCfg.showFinalHistogram = false;
baseCfg.saveData = true;
baseCfg.videoName = 'networkA_formation_v2_cleaned';
baseCfg.dataFile = 'Network_PartA_v2_Cleaned_Data.mat';
% ------------------------------------------------------

% Define per-run overrides (edit this list as needed).
% Example below sweeps xA.
runOverrides = {
    struct('xA', 0.03), ...
    struct('xA', 0.05), ...
    struct('xA', 0.08), ...
    struct('xA', 0.10)
};

batchRoot = fullfile(pwd, 'batch_runs_pkg_cleaned');
[batchDir, batchName] = createUniqueRunFolder(batchRoot, 'batch');
fprintf('Batch folder: %s\n', batchDir);

nRuns = numel(runOverrides);
summary = table('Size', [nRuns, 10], ...
    'VariableTypes', {'double','string','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'runIndex','runName','xA','Ntotal','fA','nNodes','nBonds','nActiveEnds','ssaTime','runtimeSeconds'});

allResults = cell(nRuns,1);

for k = 1:nRuns
    cfg = mergeConfig(baseCfg, runOverrides{k});

    runLabel = sprintf('run_%03d_xA_%s_N_%d', ...
        k, strrep(sprintf('%.4f', cfg.xA), '.', 'p'), cfg.Ntotal);
    [runDir, runName] = createUniqueRunFolder(batchDir, runLabel);
    cfg.outputDir = runDir;

    fprintf('\n=== BATCH RUN %d/%d: %s ===\n', k, nRuns, runName);
    results = runNetworkA_cleaned(cfg);
    allResults{k} = results;

    save(fullfile(runDir, 'results_struct.mat'), 'results');

    summary.runIndex(k) = k;
    summary.runName(k) = string(runName);
    summary.xA(k) = cfg.xA;
    summary.Ntotal(k) = cfg.Ntotal;
    summary.fA(k) = cfg.fA;
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

