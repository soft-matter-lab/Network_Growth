function style = networkDumpHistogramStyle(source, options)
% Shared bins and axes for completed-chain histograms over all supplied states.
% source: run/dumps directory, a snapshot struct, or a cell array of snapshots.
% Directory inputs are scanned one snapshot at a time to bound memory use.
if nargin < 2, options = struct(); end
options = mergeConfig(struct('kuhnXLim',[],'endToEndXLim',[], ...
    'histogramYLim',[],'kuhnYLim',[],'endToEndYLim',[],'endToEndBins',30),options);
validateattributes(options.endToEndBins,{'numeric'},{'scalar','integer','positive','finite'});
names = {'kuhnXLim','endToEndXLim','histogramYLim','kuhnYLim','endToEndYLim'};
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
maxN = 1; maxR = 1;
for k = 1:count
    s = getSnapshot(source,files,dumpDir,k);
    if isfield(s,'completed')
        maxN = max([maxN;s.completed.kuhnCount(:)]);
        maxR = max([maxR;s.completed.endToEnd(:)]);
    end
end
style = struct();
style.kuhnXLim = options.kuhnXLim;
if isempty(style.kuhnXLim), style.kuhnXLim = [0 ceil(maxN)+1]; end
style.endToEndXLim = options.endToEndXLim;
if isempty(style.endToEndXLim), style.endToEndXLim = [0 1.05*maxR]; end
style.nEdges = (floor(style.kuhnXLim(1))-.5):1:(ceil(style.kuhnXLim(2))+.5);
style.rEdges = linspace(style.endToEndXLim(1),style.endToEndXLim(2),options.endToEndBins+1);
style.histogramYLim = options.histogramYLim;
if isempty(style.histogramYLim)
    peakN = 0; peakR = 0;
    for k = 1:count
        s = getSnapshot(source,files,dumpDir,k);
        if isfield(s,'completed')
            peakN = max([peakN,histcounts(s.completed.kuhnCount,style.nEdges)]);
            peakR = max([peakR,histcounts(s.completed.endToEnd,style.rEdges)]);
        end
    end
    style.kuhnYLim = [0 max(1,ceil(1.1*peakN))];
    style.endToEndYLim = [0 max(1,ceil(1.1*peakR))];
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
