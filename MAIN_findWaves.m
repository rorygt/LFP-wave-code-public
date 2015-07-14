% Main script to find patterns in the phase velocity field of LFP data

% clear all

%% Specify all parameters

% Flag indicating whether to recalculate and write over existing files if
% they have previously been created.
rewriteFiles = 0;

% MAIN1: INPUT FILE PARAMETERS
% Maximum length of recording to store (in seconds). Files with longer
% recording will be split into multiple files.
maxFileLengthSecs = 60+10;
% Flag to process single unit spikes if possible, otherwise will be
% multi unit
singleUnitFlag = 0;
% Type of recording stimulus protocol to look at. Options: 'Spontaneous',
% 'dotDirAndCorr', 'movingDotFine', 'CircleCorrectVersionOneFreq', 
% 'dotDirection'.
stimulusProtocol = 'short';
% % Index of file to load (see list in switch statement below)
fileIndex = 1;
% Location to save output file in. To save in the current directory, just
% use an empty string
outputLoc = '.';
% Flag to indicate whether program should ask for user input. If looping or
% running this program remotely, set to 0.
askForInput = 0;

% MAIN2: FILTERING PARAMETERS
% Flag to specify whether data should be phase randomized to create
% surrogate data
surrogateFlag = 0;
% Standard deviation of random phase noise to add to each channel
% independently. If set to 0, all channels will be randomized in the same
% way. Only applies if SURROGATEFLAG == 1.
noiseSTD = pi/2;
% Choose limits for band pass filtering in Hz
fLow = 1;
fHigh = 4;
% Amount of time to discard (in seconds) from the start and end of the
% recording
discardTimeSecs = 5;

% MAIN3: OPTICAL FLOW AND PATTERN DETECTION PARAMETERS

% Flag to indicate whether to calculate velocity field. Leave this at 1 for
% full analysis, set to 0 for fast, preliminary analysis
calculateVelocityField = 1;

% OPTICAL FLOW PARAMETERS
% Smoothness weighting parameter: alpha > 0, higher alpha gives a smoother
% velocity field
alpha = 10;
% Charbonnier penalty function scale
beta = 0.1;
% Display a message on the terminal once every NSTEPSDISPLAY time steps
% have been processed in creating the velocity field
nStepsDisplay = 10240;

% SYNCHRONIZATION AND TRAVELLING WAVE PARAMETERS
% Minimum order parameter for a plane wave to be detected
planeWaveThreshold = 0.85;
% Maximum velocity field magnitude for synchrony to be detected
% If set to zero, default is mean - 1*std over all time
synchronyThreshold = 0;

% PATTERN DETECTION PARAMETERS
% Minimum duration of a pattern for it to be stored (in seconds)
minDurationSecs = 0.01;
% Maximum duration between critical points (or synchrony/plane
% waves) for them to be counted as the same pattern (in seconds)
maxTimeGapSecs = 0.005;
% Maxiumum displacement between critical points between time steps for them
% to be counted as the same pattern (measured in grid spaces)
maxDisplacement = 1;
% Minimum spatial radius for a critical point to occupy for it to be
% counted, quantified by the winding number and dicergence/curl (in grid
% spaces)
minCritRadius = 2;


%% Choose input file
% Choose which file to open, each has unreliable channels that have been
% previously found
singleUnitExists = 0;
badChannels = [];
switch stimulusProtocol
    case 'Spontaneous'
        switch fileIndex
            case 1
                % 20 mins; Good set, but a number of channels that are out of sync with the
                % rest (71, 81, 82 at least)
                experiment = 'MA026';
                file = 5;
                badChannels = [80 86 90 96 71 81 82];
                
            case 2
                % 5 mins; Relatively chaotic but with good localized patterns
                experiment = 'MA027';
                file = 32;
                badChannels = [80 86 90 96];
                
            case 3
                % 20 mins; Standard set, good clear travelling waves, given to Paul
                experiment = 'my144';
                file = 101;
                badChannels = [28 44 77];
                singleUnitExists = 1;
                
            case 4
                % 20 mins
                experiment = 'MY147';
                file = 2;
                badChannels = [80 86 90 96];
                
            case 5
                % 5 mins; Other standard set
                experiment = 'MY147';
                file = 53;
                badChannels = [80 86 90 96];
                singleUnitExists = 1;
                
            case 6
                % Huge line of broken electrodes, not good
                experiment = 'my145';
                file = 10;
                badChannels = 33:64;
        end
    case 'movingDotFine'
        switch fileIndex
            case 1
                experiment = 'MA027';
                file = 41;
            case 2
                experiment = 'MY147';
                file = 42;
        end
        
    case 'dotDirAndCorr'
        switch fileIndex
            case 1
                experiment = 'my144';
                file = 111;
            case 2
                experiment = 'MY147';
                file = 31;
            case 3
                experiment = 'MY147';
                file = 5;
        end
        
    case 'CircleCorrectVersionOneFreq'
        switch fileIndex
            case 1
                experiment = 'MA027';
                file = 23;
            case 2
                experiment = 'MY147';
                file = 41;
        end
        
    case 'dotDirection'
        switch fileIndex
            case 1
                experiment = 'MY147';
                file = 3;
            case 2
                experiment = 'my144';
                file = 104;
            case 3
                experiment = 'my144';
                file = 106;
            case 4
                experiment = 'my144';
                file = 107;
            case 5
                experiment = 'MY145';
                file = 4;
            case 6
                experiment = 'MY145';
                file = 9;
        end
        
    case 'dirAndSpeed'
        switch fileIndex
            case 1
                experiment = 'my144';
                file = 110;
        end
        
    case 'gratDirAndCorr'
        switch fileIndex
            case 1
                experiment = 'MY145_NN';
                file = 113;
            case 2
                experiment = 'MY145';
                file = 114;
        end
    case 'short'
        switch fileIndex
            case 1
                experiment = 'shortMY147';
                file = 53;
        end
        
    otherwise
        error('Invalid file selection!')
        
end

% Only look for single units if they exist
singleUnitFlag = singleUnitFlag & singleUnitExists;

% Add corner channels as bad
badChannels = [badChannels, 1, 10, 91, 100];

% Set up path to raw data files and their name
if ~singleUnitFlag
%     inputFileName = sprintf(... % File path on STELLAR linux workstation
%         '///media/RORYUNI/MEA_data/%s/%s-%d/%s_Utah100-%d_',...
%         stimulusProtocol,experiment,file,experiment,file);
       inputFileName = sprintf('./%s-%d/%s_Utah100-%d_', ...
           experiment, file, experiment, file);
      spikeType = '';
else
    inputFileName = sprintf(...
        '///media/RORYUNI/MEA_data/%s/%s-%d_SingleUnit/%s_Utah100-%d_',...
        stimulusProtocol,experiment,file,experiment,file);
    spikeType = 'Sorted';
end

%% MAIN1: Extract LFPs and spikes from raw data files
% Define output file name
tic
outputFile1Prefix = strcat(outputLoc, '/LFPsHilbert', spikeType, 'Spikes_');
outputFile1Name = sprintf('%s%s-%d.mat', outputFile1Prefix, experiment, file);
% Check if processed file already exists
if ~rewriteFiles && exist(outputFile1Name, 'file')
    fprintf('Processed LFP file already exists, loading %s...\n', ...
        outputFile1Name)
    load(outputFile1Name)   
else
    % Check if file already exists split into different parts
    experimentPart = strcat(experiment, 'a');
    outputFile1NamePart = sprintf('%s%s-%d.mat', ...
        outputFile1Prefix, experimentPart, file);
    if ~rewriteFiles && exist(outputFile1NamePart, 'file')
        % Find how many different parts there are
        for isegment = 1:26
            iexperimentPart = strcat(experiment, char(isegment+'a'));
            ioutputFile1NamePart = sprintf('%s%s-%d.mat', ...
                    outputFile1Prefix, iexperimentPart, file);
                if ~exist(ioutputFile1NamePart, 'file')
                    break
                end
        end
        display('Processed LFP file already exists.')
        fprintf('Please choose a file between a and %s', char(isegment-1+'a'))
    else
        % Actually process raw files
        clearvars experimentPart outputFile1NamePart
        display('Pre-processing recording file...')
        main1ProcessDataFiles
        fprintf('LFPs and spikes saved to %s.\n', outputFile1Name)
    end
end
    
% Check if file has been split into multiple parts
loopAll = 0;
if exist('isegment', 'var')
    lastFile = char(isegment-1+'a');
    display('Recording length was too long, split into multiple files.')
    if askForInput
        fprintf('Please enter a file letter between a and %s to process,\n', ...
            lastFile)
        chosenFile = input(...
            'or just hit enter to process all files sequentially.', 's');
        if isempty(chosenFile)
            display('Analysis will sequentially loop over all files.')
            loopAll = 1;
        elseif ~isempty(regexp(chosenFile, sprintf('^[a-%s]$', lastFile)))
            experiment = strcat(experiment, chosenFile);
            outputFile1Name = sprintf('%s%s-%d.mat', outputFile1Prefix, experiment, file);
            fprintf('Opening file %s\n', outputFile1Name)
            load(outputFile1Name)
        end
    else
        loopAll = 1;
    end
end
toc

%% MAIN2: BAND-PASS FILTER LFPS AND EXTRACT PHASE AND AMPLITUDE
nLoops2 = 1;
if loopAll
    origExperiment = experiment;
    nLoops2 = isegment;
end

% Loop over analysis steps MAIN2 and MAIN3 multiple times if required
for idata2 = 1:nLoops2  %%%% CHANGE THIS BACK TO 1
if loopAll
    experiment = strcat(origExperiment, char(idata2-1+'a'));
    outputFile1Name = sprintf('%s%s-%d.mat', outputFile1Prefix, experiment, file);
    fprintf('Opening file %s\n', outputFile1Name)
    load(outputFile1Name, 'recordingLength', 'LFPs', 'spikeMatrix', ...
        'file', 'Fs')
end

tic
if surrogateFlag
    surrogateString = 'surr';
else
    surrogateString = '';
end
outputFile2Name = sprintf('%s/filteredLFPsHilbert_%d-%dHz_%s-%d%s.mat', ...
    outputLoc, fLow, fHigh, experiment, file, surrogateString);

% Check if file already exists
output2Exists = 0;
if ~rewriteFiles && exist(outputFile2Name, 'file')
    thisfLow = fLow; thisfHigh = fHigh; thisdiscardTimeSecs = discardTimeSecs;
    fprintf('Filtered LFPs and Hilbert transformation file already exists!\nLoading %s...\n', ...
        outputFile2Name)
    load(outputFile2Name, 'fLow', 'fHigh', 'discardTimeSecs', ...
        'filteredLFPs', 'hilbertLFPs', 'hilbertSection', ...
        'recordingLength', 'phase', 'amplitude')
    
    % TODO: display a warning if filtering parameters have been changed
        
else
    display('Band-pass filtering LFPs and taking Hilbert transformation...')
    main2FilterAndHilbert
    fprintf('Filtered LFPs and Hilbert transformation saved to %s.\n', outputFile2Name)
end
toc

%% MAIN3: CALCULATE PHASE VELOCITY FIELD AND FIND PATTERNS
clearvars spikeMatrix amplitude filteredLFPs hilbertLFPs LFPs
if calculateVelocityField
    time = clock;
    outputFile3Name = sprintf('%s/detectedPatterns_%d-%dHz_%s-%d%s_%s_%d%d.mat', ...
        outputLoc, fLow, fHigh, experiment, file, surrogateString, date, ...
        time(3), time(4));
    main3VelocityFieldAndPatterns
    printPatternStatistics(patterns)
end
end