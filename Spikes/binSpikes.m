function [spikeCounts, bins] = binSpikes(spikeMatrix, binSize, combineChannels)
% [spikeCounts bins] = BINSPIKES(spikeMatrix, binSize)
% SPIKEMATRIX is a C x T binary matrix where C is the number of channels
% and T is the number of time steps with a 1 indicating a spike occuring at
% that instant and location. This function will bin spikes to bins of size
% BINSIZE (measured in time steps). If flag COMBINECHANNELS is set to 1,
% different channels will be summed to create a 1D vector output.

if nargin < 3
    combineChannels = 0;
end

% Define the start of each bin and remove the last entry if less than a
% full bin
bins = 1:binSize:(size(spikeMatrix, 2) + binSize - 1);
[~, ind] = histc(1:size(spikeMatrix,2), bins);

% Loop over every channel and bin spikes using ACCUMARRAY
spikeCounts = zeros(size(spikeMatrix,1), length(bins) - 1);
for ichannel = 1:size(spikeMatrix,1)
    spikeCounts(ichannel, :) = ...
        accumarray(ind', spikeMatrix(ichannel, :), [], @sum);
end

% Crop off last entry of bins to make it equal lenght
bins = bins(1:(end-1));

% Combine channels if specified
if combineChannels == 1
    spikeCounts = sum(spikeCounts, 1);
end
