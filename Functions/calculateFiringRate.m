function [firingRate, nActive] = calculateFiringRate(spikeMatrix, halfWindowSize, Fs, groupFiringFlag, singleUnitFlag, minRate)
% Calculates group firing rate based on SPIKEMATRIX using a sliding time
% window of length WINDOWSIZE. GROUPFIRINGFLAG is a binary flag indicating
% whether to calculate the group firing rate across all channels (1,
% default) or the individual firing rate for each channel (0).
% SINGLEUNITFLAG is a binary flag indicating whether to look for
% single-unit spikes (1, default) or multi-unit.
% MINRATE is an optional argument in Hertz that specifies a minimum average
% firing rate for a channel to be included in the group firing rate. If
% GROUPFIRINGRATE == 0, this argument does nothing.

if nargin < 4
    groupFiringFlag = 1;
end

if nargin < 5
    singleUnitFlag = 1;
end

if nargin < 6
    minRate = 0;
end

firingRate = zeros(size(spikeMatrix));
if ~singleUnitFlag
    spikeMatrix = -spikeMatrix;
end
spikeMatrix(spikeMatrix<0) = 0;

for itime = (halfWindowSize+1) : (size(spikeMatrix,2)-halfWindowSize-1)
    spikeWindow = spikeMatrix(:, (-halfWindowSize:halfWindowSize) + itime);
    firingRate(:, itime) = mean(spikeWindow, 2);
end

% Multiply by the sampling frequency to give a firing rate in Hz
firingRate = firingRate * Fs;
nActive = sum(sum(abs(spikeMatrix), 2) > 0);

% Turn the individual firing rates into a group firing rate if required
if groupFiringFlag
    % Divide by the number of channels that were actually active
    activeCh = find(mean(firingRate,2) > minRate);
    nActive = length(activeCh);
    firingRate = nanmean(firingRate(activeCh,:));
end