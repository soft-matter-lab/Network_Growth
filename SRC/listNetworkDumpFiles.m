function [files, metadata, dumpDir] = listNetworkDumpFiles(runOrDumpDir)
% Locate completed snapshots, ordered numerically, without loading their payloads.
dumpDir = char(runOrDumpDir);
if isfolder(fullfile(dumpDir,'dumps')), dumpDir = fullfile(dumpDir,'dumps'); end
metaFile = fullfile(dumpDir,'metadata.mat');
if ~isfile(metaFile)
    error('NetworkGrowth:NoDumps','No metadata.mat in %s. Select one run or its dumps folder.',dumpDir);
end
m = load(metaFile,'metadata'); metadata = m.metadata;
dumpDir = fileparts(metaFile); % remove a trailing separator before finding the run parent
files = dir(fullfile(dumpDir,'snapshot_*.mat'));
numbers = cellfun(@(name)sscanf(name,'snapshot_%d.mat',1),{files.name});
[~,order] = sort(numbers); files = files(order);
end
