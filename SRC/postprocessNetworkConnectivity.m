function report=postprocessNetworkConnectivity(inputFolder,options)
% One saved run/dumps folder OR a recursively discovered batch of final runs.
% Shared by the growing/preformed drivers. No simulation state is modified.
if nargin<2,options=struct();end
options=mergeConfig(struct('outputRoot','','expectedGrowthMode','', ...
    'saveSnapshots',true,'snapshotResolution',600,'showNodeConnectivity',true, ...
    'rankTolerance',1e-9,'pairTolerance',1e-8,'imageRange',2, ...
    'checkRankStability',true,'continueOnError',true),options);
validateattributes(options.snapshotResolution,{'numeric'},{'scalar','integer','positive'});
folder=char(inputFolder);assert(isfolder(folder),'NetworkGrowth:ConnectivityFolder','Select a saved run or batch folder.');
[~,attr]=fileattrib(folder);folder=attr.Name;[parent,name]=fileparts(folder);
if strcmpi(name,'dumps'),folder=parent;end
folders=findRuns(folder);
if isempty(folders),error('NetworkGrowth:NoConnectivityRuns','No saved simulation folders found in %s.',folder);end
outputRoot=options.outputRoot;if isempty(outputRoot),outputRoot=folder;end
[outputDir,~]=createUniqueRunFolder(outputRoot,'connectivity_final');
metricsOptions=struct('rankTolerance',options.rankTolerance,'pairTolerance',options.pairTolerance, ...
    'imageRange',options.imageRange,'checkRankStability',options.checkRankStability);
summary=struct([]);
for i=1:numel(folders)
    [~,base,ext]=fileparts(folders{i});label=[base ext];
    sampleDir=fullfile(outputDir,sprintf('%03d_%s',i,label));
    row=struct('runDirectory',string(folders{i}),'outputDirectory',string(sampleDir), ...
        'status',"failed",'growthMode',"",'f',NaN,'k_self',NaN,'seed',NaN, ...
        'xA',NaN,'Nscale',NaN,'selfLoopsAllowed',NaN, ...
        'junctions',NaN,'physicalStrands',NaN,'selfLoops',NaN, ...
        'N_path_x',NaN,'N_path_y',NaN,'N_rigid',NaN,'wholeNetworkRigid',false, ...
        'largestRigidClusterFraction',NaN,'uncoveredRegions',NaN,'message',"");
    try
        [s,source]=loadFinalNetworkState(folders{i});c=s.config;
        if ~isempty(options.expectedGrowthMode)&&~strcmpi(c.growthMode,options.expectedGrowthMode)
            error('NetworkGrowth:ConnectivityMode','This is a %s run; use that mode driver.',c.growthMode);
        end
        fprintf('Final connectivity %d/%d: %s\n',i,numel(folders),folders{i});
        result=networkConnectivityMetrics(s.nodes,s.bonds,s.L,metricsOptions);
        result.source=source;result.growthConfig=c;
        mkdir(sampleDir);
        save(fullfile(sampleDir,'connectivity.mat'),'result','-v7');
        writeJson(fullfile(sampleDir,'connectivity.json'),result);
        exportMemberships(sampleDir,result,size(s.nodes,1));
        if options.saveSnapshots,plotNetworkConnectivity(s,result,sampleDir,options);end
        row.status="completed";row.growthMode=string(c.growthMode);
        row.f=configValue(c,'fA');row.k_self=configValue(c,'kSelfLoop');row.seed=configValue(c,'randSeed');
        row.xA=configValue(c,'xA');row.Nscale=configValue(c,'NtotalScale');
        row.selfLoopsAllowed=configValue(c,'allowSelfLoops');
        row.junctions=result.junctionCount;row.physicalStrands=result.physicalBondCount;row.selfLoops=result.selfLoopCount;
        row.N_path_x=result.N_path_x;row.N_path_y=result.N_path_y;row.N_rigid=result.N_rigid;
        row.wholeNetworkRigid=result.wholeNetworkRigid;
        row.largestRigidClusterFraction=result.rigidity.largestRigidClusterFraction;
        row.uncoveredRegions=result.rigidity.uncoveredRegionCount;
        fprintf('  N_path(x,y)=(%d,%d), N_rigid=%d, whole-network rigid=%d\n', ...
            row.N_path_x,row.N_path_y,row.N_rigid,row.wholeNetworkRigid);
    catch err
        row.message=string(err.message);
        if ~options.continueOnError,rethrow(err);end
        warning('NetworkGrowth:ConnectivityRunFailed','%s: %s',folders{i},err.message);
    end
    if isempty(summary),summary=row;else,summary(i,1)=row;end
    writetable(struct2table(summary),fullfile(outputDir,'connectivity_summary.csv'));
end
report=struct('outputDir',outputDir,'summary',struct2table(summary),'options',options, ...
    'completed',nnz([summary.status]=="completed"),'failed',nnz([summary.status]=="failed"));
save(fullfile(outputDir,'connectivity_summary.mat'),'report','-v7');
writeJson(fullfile(outputDir,'connectivity_summary.json'),summary);
fid=fopen(fullfile(outputDir,'README.txt'),'w');assert(fid>0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,['FINAL-STATE CONNECTIVITY: GROWING OR PREFORMED CHAINS\n' ...
    'Only saved final states are selected, never the first 100%%-usage checkpoint.\n' ...
    'N_path_x/y: maximum total periodic winding with unit strand capacity, then\n' ...
    'minimum total used Kuhn count. Routes may share nodes, never bonds. A circuit\n' ...
    'winding twice contributes two; N_path is not necessarily the number of circuits.\n' ...
    'N_rigid: maximal certified first-order rigid clusters with >=3 noncollinear\n' ...
    'junctions. Two-node bars are listed separately. Clusters can share junctions.\n' ...
    'wholeNetworkRigid means one cluster contains EVERY formed junction.\n' ...
    'The cell is fixed and displacement fields repeat after one periodic cell.\n' ...
    'This is bar-constraint geometry, not stiffness, prestress/contact rigidity,\n' ...
    'or variable-cell mechanics. No simulation, relaxation or modulus is computed.\n' ...
    'Self-loops and dangling chains are excluded from inter-junction constraints.\n' ...
    'Uncovered regions are connected components outside counted clusters, not modes.\n\n' ...
    'Each sample: connectivity.mat/json with certificates and memberships;\n' ...
    'paths_x/y.csv with circuit IDs, original bond rows and directions;\n' ...
    'rigid_cluster_nodes.csv, rigid_clusters.csv, junctions.csv. IDs are 1-based.\n' ...
    'Optional paths_x/y.png show ALL selected paths at equal width (red),\n' ...
    'with other strands gray and self-loops annotated orange; rigid_clusters.png\n' ...
    'shows all within-cluster strands blue and uncovered strands pale gray.\n' ...
    'Numerical failures/unfinished runs are explicit failed rows, never zero counts.\n' ...
    'Original files are unchanged. Optimization Toolbox is needed for linprog.\n']);
fprintf('Connectivity outputs: %s (%d completed, %d failed)\n',outputDir,report.completed,report.failed);
end

function folders=findRuns(folder)
% Stop at a run so raw_simulation copies/checkpoints are not counted twice.
known=isfile(fullfile(folder,'results_struct.mat'))||isfile(fullfile(folder,'final_state.mat'))|| ...
    isfolder(fullfile(folder,'dumps'))||isfile(fullfile(folder,'Network_PartA_v2_Cleaned_Data.mat'));
if ~known
    mat=dir(fullfile(folder,'*.mat'));
    for k=1:numel(mat)
        vars=whos('-file',fullfile(mat(k).folder,mat(k).name));names={vars.name};
        if all(ismember({'nodes','bonds','activeEnds'},names)),known=true;break;end
    end
end
if known,folders={folder};return;end
folders={};dd=dir(folder);
for i=1:numel(dd)
    name=dd(i).name;
    if ~dd(i).isdir||startsWith(name,'.')||startsWith(name,'connectivity_')|| ...
            ismember(name,{'engine_snapshot','SRC','node_modules'}),continue;end
    folders=[folders findRuns(fullfile(folder,name))]; %#ok<AGROW>
end
end
function value=configValue(cfg,name)
value=NaN;
if isfield(cfg,name)&&(isnumeric(cfg.(name))||islogical(cfg.(name)))&&isscalar(cfg.(name))
    value=double(cfg.(name));
end
end
function writeJson(file,data)
fid=fopen(file,'w');assert(fid>0);cleanup=onCleanup(@()fclose(fid));fwrite(fid,jsonencode(data),'char');
end
function exportMemberships(folder,r,V)
for a=1:2
    if a==1,p=r.pathsX;name='x';else,p=r.pathsY;name='y';end
    data=zeros(0,6);
    for k=1:numel(p.cycles)
        c=p.cycles{k};n=numel(c.edgeIndices);
        data=[data;repmat(k,n,1),(1:n)',c.sourceBondRows,c.edgeDirections,repmat(c.winding,n,1)]; %#ok<AGROW>
    end
    writetable(array2table(data,'VariableNames',{'circuit_id','order','source_bond_row','direction','winding_x','winding_y'}), ...
        fullfile(folder,['paths_' name '.csv']));
end
membership=zeros(0,2);clusters=zeros(0,4);
for k=1:numel(r.rigidity.clusters)
    ids=r.rigidity.clusters{k};membership=[membership;repmat(k,numel(ids),1),ids(:)]; %#ok<AGROW>
    clusters(end+1,:)=[k numel(ids) numel(r.rigidity.clusterSourceBondRows{k}) r.rigidity.clusterWindingRank(k)]; %#ok<AGROW>
end
writetable(array2table(membership,'VariableNames',{'cluster_id','junction_id'}),fullfile(folder,'rigid_cluster_nodes.csv'));
writetable(array2table(clusters,'VariableNames',{'cluster_id','junctions','internal_strands','winding_rank'}),fullfile(folder,'rigid_clusters.csv'));
writetable(table((1:V)',r.rigidity.junctionClusterMultiplicity,'VariableNames', ...
    {'junction_id','rigid_cluster_memberships'}),fullfile(folder,'junctions.csv'));
end
