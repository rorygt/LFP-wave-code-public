function [LFPs, Fs, spikeMatrix, spikeFs] = preprocessDataFile(fileName, singleUnitFlag)

% Function to open a raw data file and preprocess into a more usable
% format. INPUTFILENAME and OUTPUTFILENAME must include the paths to each
% file and must not include '.mat' on the end of the file.
%   - Loads files in located at INPUTFILENAME, where files must be in the
%     format strcat('FILEPREFIX','chN.mat') for channels N=1,2,3,...,100. 
%   - Load LFPs into a NxT array where N is the number of channels and T is
%     time and subtract the mean potential from each channel.
%   - Load spike times and create a NxT matrix SPIKEMATRIX with -1's
%     indicating the position of a multi-unit spike and 1, 2, ...
%     indicating the position of 1, 2, ... single-unit spikes.

% Load one file to find data length and initialize outputs
load(strcat(fileName, 'ch1.mat'), 'lfpstream')
cellRange = 1:length(lfpstream);
LFPs = zeros(100, length(cellRange)*length(lfpstream{1}));
spikeMatrix = zeros(size(LFPs));

% Store sampling frequencies
load(strcat(fileName, 'swep.mat'), 'thisSampleFreq');
Fs = thisSampleFreq(2);
spikeFs = thisSampleFreq(1);

% Loop over all channels
for ich = 1:100
    % Load spikes locations and LFPs
    ifileName = strcat(fileName, 'ch', num2str(ich), '.mat');
    
    if singleUnitFlag
        % Load LFPs and spiking data with single unit info
        load(ifileName, 'lfpstream', 'locs', 'sortCode')
        spikeTimesIndex = round(locs * Fs / spikeFs);
        % If sortCode is empty, all spikes are multi-unit
        if isempty(sortCode)
            sortCode = -1;
        else
            % Change multi-unit spikes to be -1 rather than 0
            sortCode(sortCode == 0) = -1;
        end
        spikeMatrix(ich, spikeTimesIndex) = sortCode;
    else
        % Load LFPs and spiking data without sorting
        load(ifileName, 'lfpstream', 'locs')
        spikeTimesIndex = round(locs * Fs / spikeFs);
        spikeMatrix(ich, spikeTimesIndex) = -1;
    end
    
    potential_i = cat(1,lfpstream{cellRange});
    % Rescale LFP to subtract its mean potential from the raw data
    LFPs(ich,:) = potential_i' - mean(potential_i);
    
end

end