function style = networkColorStyle(source, options)
% Fixed visual scales from functionality and the full supplied-chain mixture.
% Never infer saturation or species rank from one frame's observed bonds.
if nargin<2, options = struct(); end
style = struct('showNodeConnectivity',true,'colorPreformedBonds',true,'showSphericalNodes',true, ...
    'fA',[],'isPreformed',false,'chainLengths',[],'bondColors',zeros(0,3));
for name = {'showNodeConnectivity','colorPreformedBonds','showSphericalNodes'}
    key = name{1};
    if isfield(source,key), style.(key) = source.(key); end
    if isfield(options,key), style.(key) = options.(key); end
    validateattributes(style.(key),{'logical','numeric'},{'scalar','binary'});
end
if isfield(source,'fA'), style.fA = source.fA; end
if isempty(style.fA)
    style.showNodeConnectivity = false; % incomplete legacy state: no invented scale
else
    validateattributes(style.fA,{'numeric'},{'scalar','finite','integer','>=',2});
end
style.isPreformed = isfield(source,'growthMode') && strcmp(source.growthMode,'preformed');
if isfield(source,'preformed') && isfield(source.preformed,'chainLengths')
    style.isPreformed = true;
    lengths = source.preformed.chainLengths;
elseif isfield(source,'chainLengths')
    lengths = source.chainLengths;
else
    lengths = [];
end
if style.isPreformed && style.colorPreformedBonds && ~isempty(lengths)
    style.chainLengths = sort(unique(lengths(:)'));
    if numel(style.chainLengths)>3
        error('NetworkGrowth:ChainSpecies','Bond colors support at most three prescribed lengths.');
    end
    palette = [0 0 0;1 0 0;128/255 0 32/255]; % black, red, burgundy (#800020)
    style.bondColors = palette(1:numel(style.chainLengths),:);
end
end
