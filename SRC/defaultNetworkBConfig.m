function cfg = defaultNetworkBConfig()
% Case B: end-link fixed-length preformed chains from a finite reservoir.
cfg = defaultNetworkAConfig();
cfg.growthMode = 'preformed';
cfg.chainLengths = 20;                  % 1--3 distinct positive integer lengths
cfg.chainNumberFractions = 1;           % fractions of CHAIN COUNTS, sum = 1
cfg.chainCounts = [];                   % [] = allocate from the Kuhn feed budget
cfg.kAddition = 1;                      % area/time, already includes both chain ends
cfg.nucleation = false;                 % Case B has exactly three reaction types
cfg.videoName = 'networkB_formation';
cfg.dataFile = 'Network_Preformed_Data.mat';
end
