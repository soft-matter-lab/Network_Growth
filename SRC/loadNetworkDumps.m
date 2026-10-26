function [snapshots, metadata, index] = loadNetworkDumps(runOrDumpDir, usagePercent)
% Load raw snapshots from either a run directory or its dumps subdirectory.
% Optional usagePercent retains only snapshots marked with those milestones.
% Omit it (or pass []) to load every snapshot, including the final state.
if nargin < 2, usagePercent = []; end
[files,metadata,dumpDir] = listNetworkDumpFiles(runOrDumpDir);
snapshots = cell(numel(files),1);
keep = false(numel(files),1); count = 0;
for k = 1:numel(files)
    data = load(fullfile(dumpDir,files(k).name),'snapshot');
    if isempty(usagePercent) || any(abs(data.snapshot.usageTargets(:)-usagePercent(:)') < 1e-9,'all')
        count = count+1; snapshots{count} = data.snapshot; keep(k) = true;
    end
end
snapshots = snapshots(1:count); files = files(keep);
if isempty(files)
    index = table(); return;
end
event = cellfun(@(s)s.event,snapshots);
ssaTime = cellfun(@(s)s.ssaTime,snapshots);
usagePercent = cellfun(@(s)s.usagePercent,snapshots);
phase = string(cellfun(@(s)s.phase,snapshots,'UniformOutput',false));
usageTargets = cellfun(@(s)s.usageTargets,snapshots,'UniformOutput',false);
file = string({files.name})';
index = table(file,event,ssaTime,usagePercent,phase,usageTargets);
end
