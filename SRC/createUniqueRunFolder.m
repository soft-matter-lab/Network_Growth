function [runDir, runName] = createUniqueRunFolder(rootDir, label)
% Create a new unique output folder; never overwrites an existing run.

if nargin < 1 || isempty(rootDir)
    rootDir = pwd;
end
if nargin < 2 || isempty(label)
    label = 'run';
end

if ~exist(rootDir, 'dir')
    mkdir(rootDir);
end

label = regexprep(label, '[^A-Za-z0-9._-]+', '_');
stamp = datestr(now, 'yyyymmdd_HHMMSS_FFF');
baseName = sprintf('%s_%s', label, stamp);
runName = baseName;
runDir = fullfile(rootDir, runName);

suffix = 1;
while exist(runDir, 'dir')
    runName = sprintf('%s_%02d', baseName, suffix);
    runDir = fullfile(rootDir, runName);
    suffix = suffix + 1;
end

mkdir(runDir);
end

