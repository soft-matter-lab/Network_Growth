function state = writeNetworkDump(state, nodes, bonds, activeEnds, event, t, nA, nB, phase, stopReason)
% Capture one state for the union of usage, event-count and time triggers.
% Multiple milestones crossed by one reaction share a single snapshot.
if ~state.enabled, return; end
if nargin < 10, stopReason = ''; end
if state.initialFreeA > 0
    usage = 100 * (state.initialFreeA - nA) / state.initialFreeA;
else
    usage = NaN; % no free reservoir: percent consumption is undefined
end
targets = zeros(1,0);
while state.nextThreshold <= numel(state.thresholds)
    target = state.thresholds(state.nextThreshold);
    reached = (~isnan(usage) && usage + 1e-10 >= target) || ...
        (strcmp(phase,'initial') && target == 0);
    if ~reached, break; end
    targets(end+1) = target; %#ok<AGROW>
    state.nextThreshold = state.nextThreshold + 1;
end
eventDue = strcmp(phase,'event') && state.cfg.dumpEveryEvents > 0 && ...
    mod(event,state.cfg.dumpEveryEvents) == 0;
timeBins = zeros(1,0);
if strcmp(phase,'event') && state.cfg.dumpTimeInterval > 0
    bin = floor(t / state.cfg.dumpTimeInterval);
    if bin > state.lastTimeBin
        % Store the range, not an enormous vector after a long SSA time jump.
        timeBins = [state.lastTimeBin+1, bin];
        state.lastTimeBin = bin;
    end
end
isFinal = strcmp(phase,'final');
if isFinal && ~state.cfg.dumpSaveFinal, return; end
if ~(strcmp(phase,'initial') || isFinal || ~isempty(targets) || eventDue || ~isempty(timeBins))
    return;
end
state.count = state.count + 1;
snapshot = struct('schemaVersion',1,'index',state.count,'event',event, ...
    'ssaTime',t,'phase',phase,'stopReason',stopReason, ...
    'freeCrosslinkers',nA,'freeMonomers',nB,'usagePercent',usage, ...
    'usageTargets',targets,'eventIntervalTriggered',eventDue, ...
    'timeBinRange',timeBins,'nodeCount',size(nodes,1), ...
    'bondCount',size(bonds,1),'activeEndCount',size(activeEnds,1), ...
    'L',state.cfg.L,'b',state.cfg.b,'fA',state.cfg.fA);
if state.cfg.dumpSaveDistributions
    p = bonds(:,1); q = bonds(:,2);
    rv = nodes(q,1:2) + bonds(:,4:5)*state.cfg.L - nodes(p,1:2);
    snapshot.completed = struct('nodePairs',bonds(:,1:2), ...
        'kuhnCount',bonds(:,3),'contourLength',bonds(:,3)*state.cfg.b, ...
        'endToEnd',sqrt(sum(rv.^2,2)));
    growing = activeEnds(:,4) > 0;
    snapshot.dangling = struct('parent',activeEnds(growing,3), ...
        'kuhnCount',activeEnds(growing,4), ...
        'contourLength',activeEnds(growing,4)*state.cfg.b);
    snapshot.bareArmCount = sum(~growing);
end
if state.cfg.dumpSaveNetwork
    snapshot.nodes = nodes;
    snapshot.bonds = bonds;
    snapshot.activeEnds = activeEnds;
end
name = sprintf('snapshot_%06d.mat',state.count);
destination = fullfile(state.directory,name);
if isfile(destination)
    error('NetworkGrowth:ExistingSnapshot', 'Refusing to overwrite %s', destination);
end
% Rename only after save completes so readers ignore incomplete .part files.
partial = [destination '.part'];
save(partial,'snapshot','-v7');
[ok,msg] = movefile(partial,destination);
if ~ok, error('NetworkGrowth:DumpWrite','%s',msg); end
end
