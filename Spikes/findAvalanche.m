function [avSizes, avLengths, avStart, avEnd] = ...
    findAvalanche(spikeMatrix, binSize, plotFlag)
% Function to calculate and optionally plot avalanche distribution
% statistics based on a C x T binary matrix SPIKEMATRIX, where C denotes
% different channels and T denotes time. BINSIZE is an optional integer
% argument to combine time steps into bins of BINSIZE steps. PLOTFLAG is an
% optional flag indicating that a histogram of the avalanches should be
% plotted.
% 
% All outputs are 1 x A vectors where A is the number of avalanches
% detected.
%   avSizes - The number of spikes in each avalanche
%   avLengths - The number of time steps each avalanche lasted for
%   avStart - The start time of each avalanche (in binned time coordinates)
%   avEnd - The end time of each avalanche (in binned time coordinates)

% We don't care about which channel the spikes came from
spikeMatrix = sum(spikeMatrix ~= 0);

% Bin spike vector to the required size
binMatSize = floor(size(spikeMatrix,2) / binSize);
binnedMatrix = zeros(1, binMatSize);
for ibin = 1:binMatSize
    binIndices = ((ibin - 1)*binSize + 1) : (ibin*binSize);
    binnedMatrix(ibin) = sum(spikeMatrix(binIndices));
end

% Find the start and end of each avalanche
[avStart, avEnd] = findRuns(binnedMatrix > 0);

% Store the size and length of each avalanche
avSizes = zeros(1,length(avStart));
avLengths = avSizes;
for iav = 1:length(avStart)
    iIndex = avStart(iav):avEnd(iav);
    avSizes(iav) = sum(binnedMatrix(iIndex));
    avLengths(iav) = length(iIndex)*binSize;
end

% Plot a histogram of avalanche sizes
if nargin >= 3 && plotFlag
    hist(avSizes, 30);
end
    
% Convert avalanche start and end indices back to the unbinned time series
avStart = fix( (avStart - 1) * binSize + 1);
avEnd = fix( (avEnd - 1) * binSize + 1);
    
end