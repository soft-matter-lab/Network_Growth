function cfg = resolveNetworkConfig(cfg)
% Select route defaults without resolving or mutating the requested feed.
if nargin < 1, cfg = struct(); end
mode = 'growing';
if isfield(cfg,'growthMode'), mode = validatestring(cfg.growthMode,{'growing','preformed'}); end
if strcmp(mode,'preformed'), defaults = defaultNetworkBConfig();
else, defaults = defaultNetworkAConfig(); end
cfg = mergeConfig(defaults,cfg);
cfg.growthMode = mode;
end
