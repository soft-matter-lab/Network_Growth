function state = initNetworkDumps(cfg, requestedConfig, NA0, NB0, rEq)
% Initialize append-only snapshots. Completed snapshots survive interrupted runs.
state = struct('enabled', logical(cfg.saveDumps), 'directory', '', 'count', 0);
if ~state.enabled, return; end
if ~isempty(cfg.dumpUsagePercent)
    validateattributes(cfg.dumpUsagePercent, {'numeric'}, ...
        {'real','finite','vector','>=',0,'<=',100});
end
validateattributes(cfg.dumpEveryEvents, {'numeric'}, ...
    {'real','finite','scalar','integer','nonnegative'});
validateattributes(cfg.dumpTimeInterval, {'numeric'}, ...
    {'real','finite','scalar','nonnegative'});
if ~(cfg.dumpSaveDistributions || cfg.dumpSaveNetwork)
    error('NetworkGrowth:EmptyDump', 'Enable at least one dump payload.');
end
if isempty(cfg.outputDir)
    error('NetworkGrowth:DumpOutputDir', 'Set outputDir when saveDumps is true.');
end
state.directory = fullfile(cfg.outputDir, 'dumps');
if isfolder(state.directory) && ~isempty(dir(fullfile(state.directory, '*.mat')))
    error('NetworkGrowth:ExistingDumps', 'Use a fresh outputDir; dumps already exist: %s', state.directory);
end
if ~isfolder(state.directory), mkdir(state.directory); end
state.cfg = cfg;
state.thresholds = unique(cfg.dumpUsagePercent(:)');
state.nextThreshold = 1;
state.lastTimeBin = 0;
state.initialFreeA = NA0 - 1;
metadata = struct('schemaVersion',1,'config',cfg, ...
    'requestedConfig',requestedConfig,'initialCrosslinkers',NA0, ...
    'initialFreeCrosslinkers',NA0-1,'initialMonomers',NB0, ...
    'rEq',rEq,'initialRngState',rng, ...
    'usageDefinition','100 * (initialFreeA - nA) / initialFreeA; seed excluded', ...
    'timeSampling','First post-event state at/after interval crossing; no interpolation', ...
    'nodesColumns',{{'x_unwrapped','y_unwrapped','degree'}}, ...
    'bondsColumns',{{'node1','node2','kuhnCount','imageX','imageY'}}, ...
    'activeEndsColumns',{{'x_wrapped','y_wrapped','parent','kuhnCount'}});
save(fullfile(state.directory,'metadata.mat'), 'metadata');
end
