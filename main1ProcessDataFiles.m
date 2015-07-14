% Main script to preprocess a raw data file off the server

% Parameters should be set by the MAIN_findWaves.m file, or uncomment the
% following lines to set parameters to run this script independently.
%
% % Maximum length of recording to store (in seconds). Files with longer
% % recording will be split into multiple files.
% maxFileLengthSecs = 5.5*60;
% % Flag to process single unit spikes if possible, otherwise will be
% % multi unit
% singleUnitFlag = 1;
% % Index of file to load (see list in switch statement below)
% fileIndex = 5;
% % Location to save output file in. To save in the current directory, just
% % use an empty string
% outputLoc = './NewTest';
%
% % Set up output file name
% outputFile1Prefix = strcat(outputLoc, '/LFPsHilbert', spikeType, 'Spikes_');

% Load LFPs and spikes from file
[LFPs, Fs, spikeMatrix, spikeFs] = ...
    preprocessDataFile(inputFileName, singleUnitFlag);

% Set all bad channels to NaN, including corner reference electrodes
if exist('badChannels', 'var')
    badChannels = unique([1 10 91 100 badChannels]);
else
    badChannels = [1 10 91 100];
end
LFPs(badChannels, :) = NaN;

% Split data into multiple files if recording time is longer than
% MAXFILELENGTHSECS
if exist('maxFileLengthSecs', 'var')
    maxFileLengthSteps = fix(maxFileLengthSecs*Fs);
end
if exist('maxFileLengthSecs', 'var') && size(LFPs,2) > maxFileLengthSteps
    display('Initial recording too long, splitting into multiple files...')
    % Store original variables before they are reduced to shorter
    longLFPs = LFPs;
    longSpikeMatrix = spikeMatrix;
    origExperiment = experiment;
    % Loop over each MAXFILELENGTHSECS segment that fits into the original
    for isegment = 1:ceil(size(LFPs,2)/maxFileLengthSteps)
        shortSegment = ((isegment-1)*maxFileLengthSteps + 1) : ...
            min(size(longLFPs, 2), isegment*maxFileLengthSteps);
        recordingLength = length(shortSegment)/Fs;
        LFPs = longLFPs(:,shortSegment);
        spikeMatrix = longSpikeMatrix(:,shortSegment);
        % Save to a file with a letter appended to the experiment string
        iexperiment = strcat(origExperiment, char((isegment-1)+'a'));
        outputFile1Name = sprintf('%s%s-%d.mat', outputFile1Prefix, ...
            iexperiment, file);
        save(outputFile1Name, 'experiment', 'file', 'recordingLength', ...
            'LFPs', 'Fs', 'spikeMatrix', 'spikeFs', 'badChannels', 'outputFile1Name')
        if exist('stimulusProtocol', 'var')
            save(outputFile1Name, 'stimulusProtocol', '-append')
        end
    end
    clearvars longLFPs longSpikeMatrix
else
    % Save all variables to a file
    recordingLength = size(LFPs,2)/Fs;
    if exist('outputFile1Prefix', 'var')
        outputFile1Name = sprintf('%s%s-%d.mat', outputFile1Prefix, ...
            experiment, file);
        save(outputFile1Name, 'experiment', 'file', 'recordingLength', ...
            'LFPs', 'Fs', 'spikeMatrix', 'spikeFs', 'badChannels', 'outputFile1Name')
        if exist('stimulusProtocol', 'var')
            save(outputFile1Name, 'stimulusProtocol', '-append')
        end
    end
end