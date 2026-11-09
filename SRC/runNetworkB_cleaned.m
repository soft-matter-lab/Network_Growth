function results = runNetworkB_cleaned(cfg)
% Case B entry point: end-link preformed chains of one to three lengths.
if nargin < 1, cfg = struct(); end
cfg.growthMode = 'preformed';
results = runNetwork_cleaned(cfg);
end
