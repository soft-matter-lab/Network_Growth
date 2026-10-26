function [summary, batchDir] = runNetworkSweep(baseCfg, sweep, outputRoot)
% Sequential runs with separate folders and a summary saved after every run.
% Inspect makeNetworkSweepPlan(baseCfg,sweep) before launching a large batch.
if nargin < 2, sweep = struct(); end
if nargin < 3, outputRoot = fullfile(pwd,'batch_runs_pkg_cleaned_sweeps'); end
baseCfg = mergeConfig(defaultNetworkAConfig(),baseCfg);
[plan,sweep] = makeNetworkSweepPlan(baseCfg,sweep);
[batchDir,~] = createUniqueRunFolder(outputRoot,sweep.mode);
save(fullfile(batchDir,'batch_config.mat'),'baseCfg','sweep','plan');
writetable(plan,fullfile(batchDir,'run_plan.csv'));
summary = plan;
nRuns = height(plan);
summary.status = repmat("pending",nRuns,1);
summary.runDir = strings(nRuns,1);
summary.stopReason = strings(nRuns,1);
summary.errorMessage = strings(nRuns,1);
names = {'nNodes','nBonds','nActiveEnds','freeCrosslinkers','freeMonomers', ...
    'usagePercent','ssaEvents','ssaTime','runtimeSeconds','nDumps','rEq'};
for j = 1:numel(names), summary.(names{j}) = nan(nRuns,1); end
saveSummary(summary,batchDir);
fprintf('Batch folder: %s\n',batchDir);
for k = 1:nRuns
    cfg = baseCfg;
    cfg.xA = plan.xA(k); cfg.randSeed = plan.randSeed(k);
    cfg.NtotalScale = plan.NtotalScale(k);
    % Keep cfg.Ntotal nominal/unresolved: fill-mode rEq must not use scaled feed.
    label = sprintf('run_%03d_rep_%02d_xA_%.5g_seed_%d', ...
        k,plan.replicateIndex(k),cfg.xA,cfg.randSeed);
    [runDir,~] = createUniqueRunFolder(batchDir,label);
    cfg.outputDir = runDir;
    save(fullfile(runDir,'run_config.mat'),'cfg');
    summary.runDir(k) = string(runDir); summary.status(k) = "running";
    saveSummary(summary,batchDir);
    fprintf('\n=== SWEEP RUN %d/%d: xA = %g, seed = %d ===\n',k,nRuns,cfg.xA,cfg.randSeed);
    try
        results = runNetworkA_cleaned(cfg);
        if sweep.saveResults, save(fullfile(runDir,'results_struct.mat'),'results'); end
        summary.nNodes(k) = size(results.nodes,1);
        summary.nBonds(k) = size(results.bonds,1);
        summary.nActiveEnds(k) = size(results.activeEnds,1);
        summary.freeCrosslinkers(k) = results.nA;
        summary.freeMonomers(k) = results.nB;
        freeA0 = plan.initialCrosslinkers(k)-1;
        if freeA0 > 0, summary.usagePercent(k) = 100*(freeA0-results.nA)/freeA0; end
        summary.ssaEvents(k) = size(results.history,1);
        summary.ssaTime(k) = results.t;
        summary.runtimeSeconds(k) = results.runtimeSeconds;
        summary.nDumps(k) = results.nDumps;
        summary.rEq(k) = results.rEq;
        summary.stopReason(k) = string(results.stopReason);
        summary.status(k) = "completed";
        clear results; % do not accumulate entire networks across a large batch
    catch err
        summary.status(k) = "failed"; summary.errorMessage(k) = string(err.message);
        report = getReport(err,'extended','hyperlinks','off');
        save(fullfile(runDir,'run_error.mat'),'report');
        saveSummary(summary,batchDir);
        if ~sweep.continueOnError, rethrow(err); end
        warning('NetworkGrowth:SweepRunFailed','Run %d failed: %s',k,err.message);
    end
    saveSummary(summary,batchDir);
end
fprintf('\nBatch finished: %d completed, %d failed.\n', ...
    sum(summary.status=="completed"),sum(summary.status=="failed"));
end

function saveSummary(summary,batchDir)
save(fullfile(batchDir,'batch_summary.mat'),'summary');
writetable(summary,fullfile(batchDir,'batch_summary.csv'));
end
