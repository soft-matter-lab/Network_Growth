function info = buildLammpsIngredients(runOrDumpDir, options)
% Build LAMMPS ingredients from saved results, raw final MAT data, or snapshots.
% Automatic priority: results_struct.mat, final data MAT, final network dump.
% Output defaults to the selected simulation folder, even if saved paths are stale.
% Incomplete snapshots require the explicit allowLatestSnapshot option.
if nargin < 2, options = struct(); end
options = mergeConfig(struct('outputDir','','sourceFile','', ...
    'allowLatestSnapshot',false,'exportConfig',struct()),options);
validateattributes(options.allowLatestSnapshot,{'logical','numeric'},{'scalar','binary'});
validateattributes(options.exportConfig,{'struct'},{'scalar'});
allowed = {'lammpsKLD','lammpsNRho','lammpsRhoRange','lammpsZThickness','lammpsPeriodicPolicy'};
if ~all(ismember(fieldnames(options.exportConfig),allowed))
    error('NetworkGrowth:LammpsBuilderOptions', ...
        'exportConfig accepts only %s. Geometry and chain counts come from saved data.',strjoin(allowed,', '));
end
folder = char(runOrDumpDir);
if ~isfolder(folder)
    error('NetworkGrowth:LammpsBuilderFolder','Simulation folder not found: %s',folder);
end
[~,attr] = fileattrib(folder); folder = attr.Name;
[parent,name] = fileparts(folder);
if strcmpi(name,'dumps')
    runDir = parent; dumpDir = folder;
else
    runDir = folder; dumpDir = fullfile(folder,'dumps');
end
fallbackCfg = struct();
metaFile = fullfile(dumpDir,'metadata.mat');
if isfile(metaFile)
    m = load(metaFile,'metadata');
    if isfield(m,'metadata') && isfield(m.metadata,'config')
        fallbackCfg = m.metadata.config;
    end
end

% Inspect MAT variable names before loading potentially large unrelated files.
files = {}; kinds = {};
if ~isempty(options.sourceFile)
    file = char(options.sourceFile);
    if ~isAbsolute(file), file = fullfile(runDir,file); end
    if ~isfile(file)
        error('NetworkGrowth:LammpsBuilderSource','Source file not found: %s',file);
    end
    kind = sourceKind(file);
    if isempty(kind)
        error('NetworkGrowth:LammpsBuilderSource','No recognized saved network in %s.',file);
    end
    files = {file}; kinds = {kind};
else
    resultFile = fullfile(runDir,'results_struct.mat');
    if isfile(resultFile), files{end+1} = resultFile; kinds{end+1} = 'results'; end
    rawFiles = dir(fullfile(runDir,'*.mat'));
    rawPaths = {};
    for k = 1:numel(rawFiles)
        file = fullfile(runDir,rawFiles(k).name);
        if strcmp(sourceKind(file),'raw'), rawPaths{end+1} = file; end %#ok<AGROW>
    end
    defaultRaw = fullfile(runDir,'Network_PartA_v2_Cleaned_Data.mat');
    if ismember(defaultRaw,rawPaths)
        rawPaths = [{defaultRaw},setdiff(rawPaths,{defaultRaw},'stable')];
    elseif numel(rawPaths) > 1 && isempty(files)
        error('NetworkGrowth:LammpsBuilderAmbiguous', ...
            'Multiple final data MAT files found. Set options.sourceFile to choose one.');
    end
    files = [files,rawPaths]; kinds = [kinds,repmat({'raw'},size(rawPaths))];
    snapshots = dir(fullfile(dumpDir,'snapshot_*.mat'));
    numbers = cellfun(@(s)sscanf(s,'snapshot_%d.mat',1),{snapshots.name});
    [~,order] = sort(numbers,'descend'); snapshots = snapshots(order);
    for k = 1:numel(snapshots)
        files{end+1} = fullfile(dumpDir,snapshots(k).name); %#ok<AGROW>
        kinds{end+1} = 'snapshot'; %#ok<AGROW>
    end
end

state = []; source = struct(); reasons = {};
latest = []; latestSource = struct();
for k = 1:numel(files)
    [candidate,cfg,phase,reason] = readState(files{k},kinds{k},fallbackCfg);
    if ~isempty(reason)
        reasons{end+1} = sprintf('%s: %s',files{k},reason); %#ok<AGROW>
        continue;
    end
    candidate.config = cfg;
    candidateSource = struct('file',files{k},'kind',kinds{k},'phase',phase);
    if strcmp(kinds{k},'snapshot') && ~strcmp(phase,'final')
        if isempty(latest), latest = candidate; latestSource = candidateSource; end
        continue;
    end
    state = candidate; source = candidateSource; break;
end
if isempty(state) && ~isempty(latest)
    if ~options.allowLatestSnapshot
        error('NetworkGrowth:LammpsBuilderNoFinal', ...
            ['Only non-final network snapshots are available. To export the latest saved ' ...
             'checkpoint explicitly, set allowLatestSnapshot = true.']);
    end
    state = latest; source = latestSource;
    warning('NetworkGrowth:LammpsBuilderCheckpoint', ...
        'Using a non-final saved state (%s): %s',source.phase,source.file);
end
if isempty(state)
    details = strjoin(reasons,newline);
    error('NetworkGrowth:LammpsBuilderMissingData', ...
        ['No usable network state in %s. Need nodes, bonds with Kuhn counts and images, ' ...
         'activeEnds, and saved L and b. Distribution-only dumps or CSVs are insufficient.\n%s'],runDir,details);
end

cfg = mergeConfig(defaultNetworkAConfig(),state.config);
cfg = mergeConfig(cfg,options.exportConfig);
cfg.lammpsIncludeDanglingKuhn = true; % include all chain-bound segments, even with old export settings
cfg.outputDir = runDir;
if ~isempty(options.outputDir), cfg.outputDir = char(options.outputDir); end
reportFile = fullfile(cfg.outputDir,'lammps_build_report.mat');
if isfile(reportFile)
    error('NetworkGrowth:ExistingLammpsFiles','Refusing to overwrite %s. Use a fresh outputDir.',reportFile);
end
fprintf('Building LAMMPS ingredients from %s (%s).\n',source.file,source.phase);
info = writeLammpsNetwork(state.nodes,state.bonds,state.activeEnds,cfg);
info.sourceFile = source.file;
info.sourceKind = source.kind;
info.sourcePhase = source.phase;
info.runDir = runDir;
info.outputDir = cfg.outputDir;
info.config = cfg;
info.reportFile = reportFile;
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
save(reportFile,'info');
end

function kind = sourceKind(file)
variables = whos('-file',file); names = {variables.name};
kind = '';
if ismember('results',names)
    kind = 'results';
elseif ismember('snapshot',names)
    kind = 'snapshot';
elseif all(ismember({'nodes','bonds','activeEnds'},names))
    kind = 'raw';
end
end

function [state,cfg,phase,reason] = readState(file,kind,fallbackCfg)
state = []; cfg = fallbackCfg; phase = 'final'; reason = '';
switch kind
    case 'results'
        data = load(file,'results');
        if ~isfield(data,'results'), reason = 'missing results struct'; return; end
        state = data.results;
    case 'raw'
        variables = whos('-file',file); names = {variables.name};
        wanted = intersect(names,{'nodes','bonds','activeEnds','cfg','L','b'});
        state = load(file,wanted{:});
    case 'snapshot'
        data = load(file,'snapshot'); state = data.snapshot;
        phase = 'unknown';
        if isfield(state,'phase'), phase = char(state.phase); end
end
if ~isstruct(state) || ~isscalar(state) || ~all(isfield(state,{'nodes','bonds','activeEnds'}))
    reason = 'missing nodes, bonds, or activeEnds (full network payload required)'; return;
end
if isfield(state,'config'), cfg = mergeConfig(cfg,state.config); end
if isfield(state,'cfg'), cfg = mergeConfig(cfg,state.cfg); end
% Snapshot geometry takes precedence over metadata, which may be from an older run.
for name = {'L','b'}
    if isfield(state,name{1}), cfg.(name{1}) = state.(name{1}); end
end
if ~all(isfield(cfg,{'L','b'})) || isempty(cfg.L) || isempty(cfg.b)
    reason = 'missing saved L or b; physical geometry cannot be guessed'; return;
end
% Some older MAT files store a completely empty array without its column count.
if isempty(state.bonds), state.bonds = zeros(0,5); end
if isempty(state.activeEnds), state.activeEnds = zeros(0,4); end
end

function tf = isAbsolute(file)
tf = startsWith(file,{'/','\'}) || (numel(file)>=2 && file(2)==':');
end
