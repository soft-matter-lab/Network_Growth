function style = networkDumpHistogramStyle(source, options)
% Shared bins and axes for completed-chain histograms over all supplied states.
% source: run/dumps directory, a snapshot struct, or a cell array of snapshots.
% Directory inputs are scanned one snapshot at a time to bound memory use.
if nargin < 2, options = struct(); end
options = mergeConfig(struct('kuhnXLim',[],'endToEndXLim',[], ...
    'histogramYLim',[],'kuhnYLim',[],'endToEndYLim',[],'endToEndBins',30, ...
    'displacementLimit',[],'displacementBins',40,'displacementCLim',[], ...
    'displacementColorPercentile',98, ...
    'weightXLim',[],'weightYLim',[]),options);
validateattributes(options.endToEndBins,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(options.displacementBins,{'numeric'},{'scalar','integer','positive','finite'});
validateattributes(options.displacementColorPercentile,{'numeric'},{'scalar','real','finite','>',0,'<=',100});
if ~isempty(options.displacementLimit)
    validateattributes(options.displacementLimit,{'numeric'},{'scalar','real','positive','finite'});
end
names = {'kuhnXLim','endToEndXLim','histogramYLim','kuhnYLim','endToEndYLim', ...
    'displacementCLim','weightXLim','weightYLim'};
for j = 1:numel(names)
    limits = options.(names{j});
    if ~isempty(limits)
        validateattributes(limits,{'numeric'},{'vector','numel',2,'real','finite','nonnegative','increasing'});
        options.(names{j}) = limits(:)';
    end
end
files = []; dumpDir = '';
if isstruct(source), source = {source}; end
if iscell(source)
    count = numel(source);
else
    [files,~,dumpDir] = listNetworkDumpFiles(source);
    count = numel(files);
end
maxN = 1; maxR = 1; maxComponent = 0;
for k = 1:count
    s = getSnapshot(source,files,dumpDir,k);
    if isfield(s,'completed')
        maxN = max([maxN;s.completed.kuhnCount(:)]);
        maxR = max([maxR;s.completed.endToEnd(:)]);
    end
    dv = networkDumpDisplacements(s);
    maxComponent = max([maxComponent;abs(dv(:))]);
end
style = struct();
style.normalization = 'probability';
style.kuhnXLim = options.kuhnXLim;
if isempty(style.kuhnXLim), style.kuhnXLim = [0 ceil(maxN)+1]; end
style.endToEndXLim = options.endToEndXLim;
if isempty(style.endToEndXLim), style.endToEndXLim = [0 1.05*maxR]; end
style.nEdges = (floor(style.kuhnXLim(1))-.5):1:(ceil(style.kuhnXLim(2))+.5);
style.weightXLim = options.weightXLim;
if isempty(style.weightXLim), style.weightXLim = style.kuhnXLim; end
style.weightEdges = (floor(style.weightXLim(1))-.5):1:(ceil(style.weightXLim(2))+.5);
style.weightYLim = options.weightYLim;
if isempty(style.weightYLim)
    peakWeight = 0;
    for k = 1:count
        s = getSnapshot(source,files,dumpDir,k);
        if isfield(s,'completed')
            fraction = chainWeightFractions(s.completed.kuhnCount,style.weightEdges);
            peakWeight = max([peakWeight,fraction]);
        end
    end
    style.weightYLim = [0 max(.05,ceil(110*peakWeight)/100)];
end
style.rEdges = linspace(style.endToEndXLim(1),style.endToEndXLim(2),options.endToEndBins+1);
style.histogramYLim = options.histogramYLim;
limit = options.displacementLimit;
if isempty(limit), limit = 1.05*maxComponent; end
if limit<=0, limit = 1; end
style.displacementLim = [-limit limit];
style.displacementEdges = linspace(-limit,limit,options.displacementBins+1);
style.displacementCLim = options.displacementCLim;
style.displacementColorPercentile = options.displacementColorPercentile;
if isempty(style.displacementCLim)
    % One shared robust cap: a sparse early state must not wash out later ones.
    % Retain occupied-bin probabilities only, not the full snapshot payloads.
    occupied = cell(count,1);
    for k = 1:count
        s = getSnapshot(source,files,dumpDir,k);
        dv = networkDumpDisplacements(s);
        if ~isempty(dv)
            pxy = histcounts2(dv(:,1),dv(:,2),style.displacementEdges,style.displacementEdges)/size(dv,1);
            occupied{k} = pxy(pxy>0);
        end
    end
    values = sort(vertcat(occupied{:}));
    cap = .001; % finite fallback for runs without any visible bonds
    if ~isempty(values)
        position = 1+(numel(values)-1)*options.displacementColorPercentile/100;
        lo = floor(position); hi = ceil(position); alpha = position-lo;
        cap = values(lo)*(1-alpha)+values(hi)*alpha;
    end
    style.displacementCLim = [0 cap];
end
if isempty(style.histogramYLim)
    peakN = 0; peakR = 0;
    for k = 1:count
        s = getSnapshot(source,files,dumpDir,k);
        if isfield(s,'completed')
            n = s.completed.kuhnCount; r = s.completed.endToEnd;
            if ~isempty(n), peakN = max([peakN,histcounts(n,style.nEdges)/numel(n)]); end
            if ~isempty(r), peakR = max([peakR,histcounts(r,style.rEdges)/numel(r)]); end
        end
    end
    style.kuhnYLim = [0 max(.05,ceil(110*peakN)/100)];
    style.endToEndYLim = [0 max(.05,ceil(110*peakR)/100)];
    style.histogramYLim = [0 max(style.kuhnYLim(2),style.endToEndYLim(2))]; % legacy summary
else
    style.kuhnYLim = style.histogramYLim;
    style.endToEndYLim = style.histogramYLim;
end
if ~isempty(options.kuhnYLim), style.kuhnYLim = options.kuhnYLim; end
if ~isempty(options.endToEndYLim), style.endToEndYLim = options.endToEndYLim; end
end

function s = getSnapshot(source,files,dumpDir,k)
if iscell(source)
    s = source{k};
else
    data = load(fullfile(dumpDir,files(k).name),'snapshot'); s = data.snapshot;
end
end
