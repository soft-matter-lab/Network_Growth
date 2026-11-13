function [state,source]=loadFinalNetworkState(folder)
% Read a final dump/result/raw MAT only. Never substitute a usage checkpoint.
folder=char(folder);assert(isfolder(folder),'NetworkGrowth:FinalFolder','Folder not found.');
[~,attr]=fileattrib(folder);folder=attr.Name;
[parent,name]=fileparts(folder);if strcmpi(name,'dumps'),folder=parent;end
cfg=struct();metaFile=fullfile(folder,'dumps','metadata.mat');
if isfile(metaFile)
    d=load(metaFile,'metadata');
    if isfield(d.metadata,'config'),cfg=d.metadata.config;end
end
configFile=fullfile(folder,'growth_parameters.json');
if isfile(configFile),cfg=mergeConfig(cfg,jsondecode(fileread(configFile)));end
files={};kinds={};
if isfile(fullfile(folder,'final_state.mat'))
    files{end+1}=fullfile(folder,'final_state.mat');kinds{end+1}='snapshot';
end
dd=dir(fullfile(folder,'dumps','snapshot_*.mat'));
numbers=arrayfun(@(f)sscanf(f.name,'snapshot_%d.mat',1),dd);
[~,order]=sort(numbers,'descend');
for i=order(:)',files{end+1}=fullfile(dd(i).folder,dd(i).name);kinds{end+1}='snapshot';end
resultFile=fullfile(folder,'results_struct.mat');
if isfile(resultFile),files{end+1}=resultFile;kinds{end+1}='results';end
raw=dir(fullfile(folder,'*.mat'));rawFiles={};
for i=1:numel(raw)
    file=fullfile(raw(i).folder,raw(i).name);vars=whos('-file',file);names={vars.name};
    if all(ismember({'nodes','bonds','activeEnds'},names)),rawFiles{end+1}=file;end
end
preferred=fullfile(folder,'Network_PartA_v2_Cleaned_Data.mat');
if ismember(preferred,rawFiles),rawFiles=[{preferred},setdiff(rawFiles,{preferred},'stable')];end
if numel(rawFiles)>1 && isempty(files)
    error('NetworkGrowth:AmbiguousFinalState','Several raw final files exist in %s.',folder);
end
files=[files rawFiles];kinds=[kinds repmat({'raw'},1,numel(rawFiles))];
for i=1:numel(files)
    kind=kinds{i};file=files{i};
    switch kind
        case 'snapshot'
            d=load(file,'snapshot');
            if ~isfield(d,'snapshot')||~isfield(d.snapshot,'phase')||~strcmp(d.snapshot.phase,'final'),continue;end
            s=d.snapshot;
        case 'results'
            d=load(file,'results');s=d.results;
        case 'raw'
            s=load(file);
    end
    if ~all(isfield(s,{'nodes','bonds','activeEnds'})),continue;end
    if isfield(s,'phase')&&~strcmp(s.phase,'final'),continue;end
    effective=cfg;
    if isfield(s,'config'),effective=mergeConfig(effective,s.config);end
    if isfield(s,'cfg'),effective=mergeConfig(effective,s.cfg);end
    for key={'L','b','fA'}
        if isfield(s,key{1}),effective.(key{1})=s.(key{1});end
    end
    if ~isfield(effective,'L'),continue;end
    if ~isfield(effective,'growthMode'),effective.growthMode='growing';end
    if isempty(s.bonds),s.bonds=zeros(0,5);end
    if isempty(s.activeEnds),s.activeEnds=zeros(0,4);end
    s.L=effective.L;s.config=effective;s.phase='final';
    state=s;source=struct('file',file,'kind',kind,'phase','final','runDirectory',folder);
    return;
end
error('NetworkGrowth:NoFinalState','No complete saved final network found in %s; usage checkpoints are not final states.',folder);
end
