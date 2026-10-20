function cfg = mergeConfig(defaults, overrides)
% Fill missing fields in overrides from defaults.

cfg = defaults;
if nargin < 2 || isempty(overrides)
    return;
end

names = fieldnames(overrides);
for i = 1:numel(names)
    cfg.(names{i}) = overrides.(names{i});
end
end

