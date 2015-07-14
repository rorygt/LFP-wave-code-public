% Script to band-pass filter the LFPs between specified limits and take the
% Hilbert transform. Surrogate LFPs can also be created at this stage.

% Parameters should be set by the MAIN_findWaves.m file, or uncomment the
% following lines to set parameters to run this script independently.
%
% % Choose limits for band pass filtering in Hz
% fLow = 1;
% fHigh = 4;
% % Amount of time to discard (in seconds) from the start and end of the
% % recording
% discardTimeSecs = 5;
% % Location to save output file in. To save in the current directory, just
% % use an empty string
% outputLoc = './NewTest';

% Make surrogate data if required
if surrogateFlag
    LFPs = generateSurrogate(LFPs, 0, noiseSTD);
end

% Band pass filter LFPs between FLOW and FHIGH Hertz using an 8th order
% Butterworth filter
filteredLFPs = nan(size(LFPs));
% Don't filter channels with broken electrodes
for ichannel = setdiff(1:size(LFPs,1), badChannels)
    filteredLFPs(ichannel,:) = filterSignal(LFPs(ichannel,:),fLow,fHigh,Fs);
end

% Take Hilbert transform of filtered signal
hilbertLFPs = hilbert(filteredLFPs')';
% Cut off the ends of the signal that may experience boundary effects from
% the filtering or Hilbert transform
discardTimeSteps = fix(discardTimeSecs * Fs);
hilbertSection = discardTimeSteps:(size(hilbertLFPs,2)-discardTimeSteps);
hilbertLFPs = hilbertLFPs(:, hilbertSection);
recordingLength = size(hilbertLFPs,2) / Fs;
% Extract analytic phase and amplitude
phase = angle(hilbertLFPs);
amplitude = abs(hilbertLFPs);

% Save all variables to file
clearvars LFPs ichannel spikeMatrix spikeFs
%outputFileName2 = sprintf('%sfilteredLFPsHilbert_%d-%dHz_%s-%d.mat', ...
%    outputLoc, fLow, fHigh, experiment, file);
if exist('outputFile2Name', 'var')
    save(outputFile2Name)
end