function results = runNetworkA_cleaned(cfg)
% Compatibility entry point for Case A (growing chains).
if nargin < 1, cfg = struct(); end
cfg.growthMode = 'growing';
results = runNetwork_cleaned(cfg);
end
