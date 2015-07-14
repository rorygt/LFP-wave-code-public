% This script will filter the signals with desired cut-off frequencies,
% take the Hilbert transform of the filtered signals and then animate the
% analytic phase of all channels, with optional extra visualization at each
% step to check filtering and Hilbert transform.

alpha = 1;
beta = 0.001;
vectorScale = 1000;

shuffleFlag = 0;
surrogateFlag = 0;

% Flag for whether to load LFP data from file and filter or whether to use
% existing variables
loadNewData = 1;

% FILTERING VARIABLES
% File name to LFPs/spikes from
spikeFileName = './Processed_data/LFPsHilbertSpikes_MY147-53.mat';
% File name to load filtered LFPs from
filtFileName = './Processed_data/filteredLFPsHilbert_1-4Hz_MY147-53.mat';
% Input prefix
filePrefix = '///media/RORYUNI/MEA_data/';
% Name of stimulus protocol to use
stimulus = 'CircleCorrectVersionOneFreq';
% Alphanumerical code of experiments to load (most visualizations will only
% use EXPERIMENT1)
experiment = 'MY147';
% Numerical code of the file to load
file = 53;
% Low frequency cutoff in Hz. If fLow is 0, the filter function returns the
% original input signal.
fLow = 1;
% High frequency cutoff in Hz
fHigh = 4;

% DISPLAY VARIABLES
% Indicator of what data to visualize:
%   - 'phase' animates the phase grid only
%   - 'relativePhase' animates the phase grid relative to reference channel
%   - 'waveforms' animates phase grid and waveforms of 2 channels
%   - 'waveformAndSpectrogram' animates phase grid, spectrogram and 
%       waveforms of one channel (filtered and unfiltered)
%   - 'many' animates grids of the raw data, phase and amplitude
%   - 'phaseVelocity' animates phase as well as the optical flow velocities
%       computed between the current phase and next step as a quiver plot
%   - 'phaseAndSpikes' animates phase grid and grid of thresholded spikes
%   - 'phaseAndSpikesCOM' animates phase grid and the centre of mass of
%           spikes averaged over some time interval
%   - 'phaseVelocitySpikes' animates a phase grid, a phase velocity field
%           and a scatter plot of active spikes
%   - 'surrogates' animates phase grids of the original data and 3
%           surrogate LFPs
%   - 'spikes' animates phase grid and 3 different types of spikes (action
%           potentials, nLFPs and pLFPs)
visualizationType = 'spikes';
% Channel to use to show sample LFP traces
channel = 23;
% Second channel to display in waveforms visualization
secondChannel = 45;
thirdChannel = 67;
% Flag indicating whether to show raw data and spectrograms of the
% unfiltered and filtered signals
showFiltering = 1;
% Maximum frequency (in Hz) to display on spectrograms
maxFreq = 15;
% Flag indicating whether to also show a sample trance and amplitude of the
% Hilbert transform of the filtered data
showSampleHilbert = 1;
% Flag for whether to save an external copy of the phase animation
record = 1;
% Frame rate (frames per sec) for animation.
animationFramesPerSec = 24;
% Flag to resize images (smooth out the rough edges and interpolate bad
% channels)
resizeFlag = 1;
% Scale to resize images by
resizeScale = 2;
% Colour map to use for phase grids
%phaseCmap = pmkmp(256, 'edge');
% phaseCmap = pmkmp(256, 'CubicL');
% phaseCmap = pmkmp(256, 'Ostwald');
phaseCmap = pmkmp_new(256, 'ostwald_o');

% DATA LENGTH VARAIBLES
% The time in seconds to start the animation
startSec = 10;
%startSec = 30;
% The time in seconds to end the animation
endSec = 20;
%endSec = 37;
% The number of data point to skip to save compuatation time. Setting
% to 1 means no samples are skipped.
sampleReduction = 2;
% Flag for deleting large variables after they are not needed (LFPs,
% filteredLFPs, hilbertLFPs)
deleteLarge = 1;
% The index of the time to start the Hilbert transform (leave at 1 unless
% the startTime of the animation is very large)
hilbertStart = 1;

%% Load LFP data

% TODO: Clean this up to work better with new system
if loadNewData == 1
    if exist('spikeFileName', 'var')
        load(spikeFileName)
    else
        inputFileName = sprintf('%s%s/%s-%d/%s_Utah100-%d_', filePrefix, ...
            stimulus, experiment, file, experiment, file);
        singleUnitFlag = 0;
        fprintf('***** > Loading LFP data file, Timestamp: [%s]', datestr(now));
        main1ProcessDataFiles
    end
end

% Calculate surrogate data if required
if strcmp(visualizationType, 'surrogates')
    goodLFPs = LFPs(setdiff(1:100,allbad), :);
    LFPs = repmat(LFPs, [1, 1, 4]);
    LFPs(:,:,2) = multivariatePhaseSurrogate(LFPs(:,:,3), 0, 1);
    LFPs(:,:,3) = multivariatePhaseSurrogate(LFPs(:,:,3));
    LFPs(setdiff(1:100,allbad), :, 4) = multivariatePhaseSurrogate(goodLFPs')';
    clear goodLFPs
elseif surrogateFlag
    goodLFPs = LFPs(setdiff(1:100,allbad), :);
    goodLFPs = multivariatePhaseSurrogate(goodLFPs, pi/8);
    LFPs(setdiff(1:100,allbad), :) = goodLFPs;
    clearvars goodLFPs
end

% Calculate spiking matrices if required
if strcmp(visualizationType, 'spikes') && loadNewData == 1
    filtSize = 2;
    spikeSmoothSpan = 21;
    if min(spikeMatrix(:)) < 0
        spikeMatrix = (spikeMatrix ~= 0);
    end
    if size(spikeMatrix, 2) == 1
        spikeMatrix = reshape(spikeMatrix, 100, size(spikeMatrix, 1)/100);
    end
    pLFPMatrix = findPositiveLFPSpikes(LFPs, filtSize);
    nLFPMatrix = findPositiveLFPSpikes(-LFPs, filtSize);
    
    % Smooth spikes
    spikeRate = zeros(size(spikeMatrix));
    for ii = 1:size(spikeMatrix, 1)
        spikeRate(ii,:) = smooth(spikeMatrix(ii,:), spikeSmoothSpan);
        pLFPMatrix(ii,:) = smooth(pLFPMatrix(ii,:), spikeSmoothSpan);
        nLFPMatrix(ii,:) = smooth(nLFPMatrix(ii,:), spikeSmoothSpan);
    end
    
    % Convert to 10x10xT matrices and SI units
    spikeMatrix = vector2grid(spikeRate)*Fs;
    nLFPMatrix = vector2grid(nLFPMatrix)*Fs;
    pLFPMatrix = vector2grid(pLFPMatrix)*Fs;
    
    % Find maximum firing rates of each spike matrix
    maxRates = [max(spikeMatrix(:)), max(nLFPMatrix(:)), max(pLFPMatrix(:))];
end


% SECONDARY VARIABLES - these are determined by other variables above
% The index of the first sample to look at.
startIndex = fix(Fs*startSec);
% The index of the last sample to look at
endIndex = fix(Fs*endSec);
% The index of the time to end the Hilbert transform, currently set to 60
% seconds after the end time of the animation
hilbertEnd = endIndex + 60*Fs;

%% Band pass filtering

% Only apply a new filter if data has not already been filtered with the
% same cutoffs.

if loadNewData == 1
    % Band pass filter LFPs between FLOW and FHIGH Hertz using an 8th order
    % Butterworth filter
    filteredLFPs = nan(size(LFPs));
    % Don't filter channels with broken electrodes
    for ichannel = setdiff(1:size(LFPs,1), badChannels)
        filteredLFPs(ichannel,:) = filterSignal(LFPs(ichannel,:),fLow,fHigh,Fs);
    end
end

%% Figure preparation

% Information about data and filtering to be displayed on the animation
dataDetails = sprintf('%s-%d, Filtered %g to %g Hz.',...
    experiment,file,fLow,fHigh);

% Make new figure
close all
fig = figure('Units','normalized','Position',[0.25 0.25 .5 .5],...
    'Name',dataDetails);
set(gcf,'color','w','defaultTextFontSize', 14, 'defaultAxesFontSize', 16)

% Define section of data to display
movieSection = startIndex:sampleReduction:endIndex;
time = movieSection/Fs;


%% Show single channel LFP signal and spectrogram before and after filtering 

if showFiltering == 1
    
    if ~exist('LFPs','var') || ~exist('filteredLFPs','var')
        disp('LFP trace cannot be shown! Please run again with loadNewData = 1.')
        
    else
    
    % Plot raw LFP signal
    subplot(2,2,1)
    plot(time, LFPs(channel, movieSection))
    xlim([time(1) time(end)])
    xlabel('Time (s)')
    ylabel('Signal amplitude')
    title(sprintf('Unfiltered signal (channel %d)',channel));
    
    % Draw spectrogram of raw signal
    subplot(2,2,2)
    STFTspectrogram(LFPs(channel, movieSection), sampleReduction,...
        maxFreq,startSec);
    title('STFT Spectrogram (unfiltered data)')
    
    % Plot filtered LFP signal
    subplot(2,2,3)
    plot(time,filteredLFPs(channel, movieSection))
    xlim([time(1) time(end)])
    xlabel('Time (s)')
    ylabel('Signal amplitude')
    title(sprintf('Filtered signal (channel %d)',channel));
    
    % Draw spectrogram of filtered signal
    subplot(2,2,4)
    STFTspectrogram(filteredLFPs(channel, movieSection), sampleReduction,...
        maxFreq, startSec);
    title('STFT Spectrogram')
    
    % Add a title to the whole figure
    titleString = strcat([dataDetails, ' Press button to continue...']);
    ann = annotation('textbox', [0.05 0.9 1 0.1], 'String', titleString,...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    % Display figure until continue button is pressed
    uicontrol('String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    clf
    
    end
end


%% Apply Hilbert transform to signal and plot if required

% If necessary, take a new Hilbert transform of the LFPs
if loadNewData == 1 || (~exist('shortHilbert','var') && ~exist('phase','var'))
    fprintf('Hilbert transformation\n')
    % Take long Hilbert transform of the data
    hilbertLFPs = zeros([size(LFPs,1), length(hilbertStart:hilbertEnd), size(LFPs, 3)]);
    for ifile = 1:size(LFPs,3)
        hilbertLFPs(:,:,ifile) = ...
            hilbert(squeeze(LFPs(:,hilbertStart:hilbertEnd,ifile))')';
    end

    % Store the relevant section of LFPs, filteredLFPs and hilbertLFPs in
    % shorter variables
    shortLFPs = LFPs(:, movieSection, :);
    shortFilteredLFPs = filteredLFPs(:, movieSection, :);
    shortHilbert = hilbertLFPs(:, movieSection, :);
    
    
    % Shuffle spatial electrodes if required
%     if strcmp(visualizeData, 'surrogates')
%         shortHilbert(:,:,2) = shortHilbert(randperm(100), :, 2);
%     else
    if shuffleFlag
        shortHilbert = shortHilbert(randperm(100), :);
    end
    
    % Delete large variables to save RAM (if specified by the deleteLarge
    % flag)
    if deleteLarge
        clearvars hilbertLFPs LFPs filteredLFPs;
    end
end

% Plot analytic phase and amplitude of selected channel if required
if showSampleHilbert == 1
    
    % Find the analytic phase and amplitude
    phase = angle(shortHilbert);
    amplitude = abs(shortHilbert);
    
    % Plot filtered signal
    subplot(3,1,1)
    plot(time, shortFilteredLFPs(channel, :))
    xlim([time(1) time(end)])
    xlabel('Time (s)')
    ylabel('Signal amplitude')
    title(sprintf('Filtered signal (channel %d)',channel));
    
    % Plot analytic phase of filtered signal
    subplot(3,1,2)
    plot(time, phase(channel, :))
    xlim([time(1) time(end)])
    xlabel('Time (s)')
    ylabel('Phase (radians)')
    title('Analytic phase');
    
    % Plot analytic amplitude of filtered signal
    subplot(3,1,3)
    plot(time, amplitude(channel, :))
    xlim([time(1) time(end)])
    xlabel('Time (s)')
    ylabel('Amplitude')
    title('Analytic amplitude');
    
    % Add a title to the whole figure
    titleString = strcat([dataDetails, ' Press button to continue...']);
    ann = annotation('textbox', [0.05 0.9 1 0.1], 'String', titleString,...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    % Display figure until continue button is pressed
    uicontrol('String','Continue','Callback','uiresume(gcbf)');
    uiwait(gcf);
    clf
    
    % Delete temporary variables
    clearvars amplitude phase
    
end

% Display a blank figure for resizing if no plots have been created yet
if showFiltering ~= 1 && showSampleHilbert ~= 1
    title('Resize figure and then press button to begin movie')
    uicontrol('String','Begin','Callback','uiresume(gcbf)');
    uiwait(gcf);
    clf
end

%% Set up animation

% Set up video file for recording phase if required
if record == 1
  
    % Define elements of title
    dataDetailsShort = sprintf('%s%s_%d_%gto%gHz', stimulus, experiment,file,fLow,fHigh);
    sectionDetails = sprintf('Start%dSkip%d',startIndex,sampleReduction);
    
    % Set title and open video file
    movieTitle = strcat(dataDetailsShort,'_',visualizationType, ...
                '_', sectionDetails, '_', datestr(now,'mmdd'), '_', ...
                 datestr(now,'HHMM'), '.avi');
    vidObj = VideoWriter(movieTitle);
    vidObj.FrameRate = animationFramesPerSec;
    open(vidObj);
end

% Add button to stop animation
stopFlag = 0;
h = uicontrol('String','Stop','Callback','stopFlag = 1');

% Define annotation as title over many subplots (slows down animation)
% ann = annotation('textbox', [0 0.9 1 0.1]);
switch visualizationType

    case 'waveformAndSpectrogram'
    % Define waveform to visualize
    unfiltWaveform = shortLFPs(channel, :);
    filtWaveform = shortFilteredLFPs(channel, :);
    
    % Define axis limits for spectrogram and waveforms
    specLims = [0 0.01];
    unfiltLims = [min(unfiltWaveform) max(unfiltWaveform)];
    filtLims = [min(filtWaveform) max(filtWaveform)];
    allLims = [specLims; unfiltLims; filtLims];

    case 'waveforms'
    % Define waveforms to visualize
    channels = [channel, channel, secondChannel, thirdChannel];
    waveforms = [shortLFPs(channel, :); shortFilteredLFPs(channel,:);...
        shortFilteredLFPs(secondChannel,:);shortFilteredLFPs(thirdChannel,:)];
    
    % Define axis limits for waveforms
    waveformHandles = [];
    lims = [min(waveforms(:)) max(waveforms(:))];
    allLims = [lims; lims; lims; lims];
    
    case 'phaseVelocity'
    % Calculate phase
    phase = angle(vector2grid(shortHilbert));
    % Calculate phase velocity field
    [vx, vy, allConvSteps] = opticalFlow(phase, ...
        101-badChannels, alpha, beta);
    
    case 'phaseVelocitySpikes'
        % Calculate maximum Hilbert amplitude for plotting
        maxAmp = max(abs(shortHilbert(:)));
    
end

%% Update animation at every time step
tic
for t = 1:length(movieSection)
    
    % Take the instant step 
    hilbertStep = squeeze(shortHilbert(:,t,:));    
    hilbertStep(badChannels,:) = NaN;
    refPhase = angle(hilbertStep(channel));
    
    % Remove dead electrodes
    if resizeFlag == 1
        for isur = 1:size(hilbertStep,2)
            hilbertStep(:,isur) = interpolateDeadElectrodes(hilbertStep(:,isur));
        end
    end
    
    hilbertStep = vector2grid(hilbertStep);
    
    % Resize data if required
    if resizeFlag == 1
        hilbertStep = imresize(hilbertStep, resizeScale, 'bilinear');
    end
    
    switch visualizationType
        
        case 'phase'
            % Only animate a large analytic phase representation
            
            % Extract phase
            phaseStep = angle(hilbertStep);

            % Plot phase
            displayGrid(phaseStep, phaseCmap, [-pi pi], 1, 0);
            
            % Remove axis tick labels
            set(gca, 'XTick', [0, 2.5, 5, 7.5, 10],...
                     'YTick', [0, 2.5, 5, 7.5, 10])
            set(gca, 'XTickLabel', {'0','1','2','3','4'}, ...
                     'YTickLabel', {'4','3','2','1','0'})
            xlabel('Electrode position (mm)')
            %ylabel('Electrode spacing (mm)')
            
            % Add legends
            timeStr = sprintf('t = %0.2f s', time(t));
            nameStr = sprintf('%s-%d', experiment, file);
            limx = xlim;
            text(1.1*limx(2), 1.07*limx(2), timeStr, 'clipping', 'off', 'fontsize', 16);
            %text(1.05*limx(2), 5.5, nameStr, 'clipping', 'off', 'fontsize', 16);
            title('Instantaneous phase of filtered local field potentials')

        case 'relativePhase'
            % Animate phase relative to reference channel
            
            phaseStep = angle(hilbertStep);
            phaseStepRel = anglesubtract(phaseStep, refPhase);
            
%             % Plot absolute phase
%             subplot(1,2,1)
%             displayGrid(phaseStep, phaseCmap, [-pi pi], 1, 0);
%             % Add title
%             titleString = sprintf('Analytic phase at %0.3f s, step %d.', ...
%                 time(t), t);
%             title(titleString)
            
            % Plot relative phase
            %subplot(1,2,2)
            displayGrid(phaseStepRel, phaseCmap, [-pi pi], 1, 0);
            % Add title
            titleString = sprintf('Phase relative to channel %d at time %0.3f.', ...
                channel, time(t));
            title(titleString)
            
        case 'waveforms'
            % Show phase grid and waveforms from 2 channels
            phaseStep = angle(hilbertStep);
            waveformHandles = showPhaseAndWaveforms(phaseStep, ...
                waveforms, t, time, channels, Fs, ...
                allLims, waveformHandles, phaseCmap);
        
        case 'waveformAndSpectrogram'
            % Show phase grid as well as updating spectrogram and waveform
            % traces 
            
            phaseStep = angle(hilbertStep);
            showPhaseWaveformAndSpectrogram(phaseStep, unfiltWaveform, ...
                filtWaveform, t, time, Fs, allLims)
            
        case 'many'
            % Show a 2x2 grid of 10x10 grids of various quantities
            
            showMany(hilbertStep, shortFilteredLFPs(t,:))
            
        case 'pca'
            % Show phase grid, the first 2 principal eigenimages and the
            % reconstructed image from the the first 2 eigenimages
            
            currIndex = showPCA(longerPhase, hilbertStep, t, ...
                movieSection, time, sampleReduction, currIndex, phaseCmap, Fs);
            
        case 'phaseVelocity'
            % Show phase grid as well as a quiver plot of optical flow
            % computed between the current time step and the next step
            
            % Stop if last time step is reached, and velocity field will
            % have one less step
            if t == length(movieSection)
                break
            end
            
            % Display phase
            subplot(1,2,1)
            displayGrid(phase(:,:,t), phaseCmap, [-pi pi], 1, 0);
            % Display velocity
            subplot(1,2,2)
            centres = 0.5 : (size(phase, 1) - 0.5);
            quiver(centres, centres, vectorScale*vx(:,:,t), ...
                vectorScale*vy(:,:,t), 'Color', 'black')
            xlim([0 10]); ylim([0 10]);
            % Add title
            titleString = sprintf('Velocity field at %0.3f s, step %d.',...
                time(t), t);
            title(titleString)
            
        case 'phaseAndSpikes'
            % Show a grid of phase and another grid of thresholded spikes
            
            % Extract phase
            phaseStep = angle(hilbertStep);

            % Plot phase
            displayGrid(phaseStep, phaseCmap, [-pi pi], 1, 0);
            
            % Add title
            titleString = sprintf('Analytic phase at %0.3f s, step %d.', ...
                time(t), t);
            title(titleString)
            
            % Also show spikes
            spikeStep = spikeMatrix(:,movieSection(t));
            singleSpikes = 100 - find(spikeStep ~= 0);
            spikeX = floor(singleSpikes/10) + 0.5;
            spikeY = mod(singleSpikes, 10) + 0.5;
            hold on
            scatter(spikeX, spikeY, 100, [0.7 0.7 0.7], 'fill')
            hold off

        case 'phaseAndSpikesCOM'
            % Show a grid of phase overlayed with the centre of mass of
            % spikes
            
            halfWindowSize = 4;
            
            % Extract phase
            phaseStep = angle(hilbertStep);

            % Plot phase
            displayGrid(phaseStep, phaseCmap, [-pi pi], 1, 0);
            
            % Add title
            titleString = sprintf('Analytic phase at %0.3f s, step %d.', ...
                time(t), t);
            title(titleString)
            
            % Add spikes' centre of mass
            if t > halfWindowSize && t < (length(movieSection) - halfWindowSize)
                iwindow = (t-halfWindowSize) : (t+halfWindowSize);
                [rowCOM, colCOM] = ...
                    spikeCentreOfMass(vector2grid(spikeMatrix(:,iwindow)));
                hold on
                scatter(rowCOM, colCOM, 300, [0.7 0.7 0.7], 'fill')
                hold off
            end
            
        case 'phaseVelocitySpikes'
            % Show phase grid, velocity field and active spikes in seperate
            % subfigures
            
            % Plot phase
            subplot(2,8,5:8)
            phaseStep = angle(hilbertStep);
            displayGrid(phaseStep, phaseCmap, [-pi pi], 1, 1);
            % Change axes to show lengths rather than channel numbers
            set(gca, 'XTick', 0:2.5:10,...
                     'YTick', 0:2.5:10)
            set(gca, 'XTickLabel', {'0','1','2','3','4'}, ...
                     'YTickLabel', {'4','3','2','1','0'})
            ylabel('Electrode position (mm)')
            title(sprintf('Phase map, t=%0.3f s', time(t)))
            
            axis square
            tvelocity = movieSection(t)-section(1)+1;
            
            % Plot amplitude
%             subplot(2,2,2)
%             ampStep = abs(hilbertStep);
%             displayGrid(ampStep, bone, [0 maxAmp], 1, 1);
%             axis square

            % Plot order parameters and active patterns
            hactive = subplot(2,8,13:16)
            % Work out section of data to plot
            traceSection = tvelocity + (fix(-0.25*Fs):fix(0.75*Fs));
            plotTimes = 1:length(traceSection);
            timeLims = (traceSection([1 end]) + section(1) - 1) / Fs;
            % Plot patterns active
            imagesc(plotTimes, 1:-1:0, patternsActiveSingle(traceSection)')
            colormap(patternCMap)
            hpatcb = colorbar;
            % Plot order parameters
            hold on
            plot(plotTimes, phi(traceSection), 'k', ...
                 plotTimes, v0(traceSection), 'r')
            xlim(plotTimes([1 end]))
            ylim([0 1])
            % Add line indicating current time
            currT = fix(0.25*Fs);
            line([currT, currT], [0, 1], 'Color', [1 0 0]);
            hold off
            % Label plot
            set(hpatcb, 'YTick', 7/16 + (0:7/8:(7-7/8)), ...
                'YTickLabel', shortNames)
            set(gca, 'XTick', [])
            title('Structure analysis')
            xlabel('Time (1s)')
            ylabel('Order parameter')
            
            % Plot velocity field
            subplot(2,8,9:11)
            ivx = velocityX(:,:,tvelocity);
            ivy = velocityY(:,:,tvelocity);
            vScale = 2000;
            quiver(0.5:9.5, 0.5:9.5, vScale*ivx, vScale*ivy, 0)
            set(gca,'YDir','reverse', 'XTick', 0:2.5:10, 'YTick', 0:2.5:10, ...
                'XTickLabel', [], 'YTickLabel', []);
            axis([0 10 0 10])
            axis square
            title('Phase velocity field')
            
            % Plot spikes
            subplot(2,8,1:3)
            spikeStep = spikeMatrix(:,movieSection(t));
            singleSpikes = 100 - find(spikeStep ~= 0);
            spikeX = floor(singleSpikes/10) + 0.5;
            spikeY = mod(singleSpikes, 10) + 0.5;
            scatter(spikeX, spikeY, 50, [0 0 0], 'k', 'fill')
            set(gca, 'XTick', 0:2.5:10, 'YTick', 0:2.5:10, ...
                'XTickLabel', {'0','1','2','3','4'}, ...
                'YTickLabel', {'0','1','2','3','4'})
            ylabel('Electrode position (mm)')
            axis([0 10 0 10])
            axis square
            box on
            title('Multi-unit activity')
            
            % Add overall title
            %suptitle('Full analysis of neural array recordings')
            
            % Fix up active patterns plot (which gets messed up by
            % suptitle)
            
        case 'surrogates'
            % Show 4 phase grids of different surrogate data methods
            for isur = 1:4
                subplot(2,2,isur)
                displayGrid(angle(hilbertStep(:,:,isur)), phaseCmap, [-pi pi], 0, 1);
                axis square
                switch isur
                    case 1
                        title(sprintf('Original phase map, t=%0.3f s', time(t)))
                    case 2
                        title('Spatially shuffled phase map')
                    case 3
                        title('Surrogate with spatial correlations')
                    case 4
                        title('Surrogate with temporal correlations')
                end
            end
            
        case 'spikes'
            % Show moving averaged spikes (along with phase) for action
            % potentials, nLFPs and pLFPs
            subplot(2,2,1)
            displayGrid(angle(hilbertStep), phaseCmap, [-pi pi], 1, 1);
            axis square
            title('Delta-band LFP phase')
            
            realIndex = movieSection(t);
            subplot(2,2,2)
            displayGrid(spikeMatrix(:,:,realIndex), hot, [0, maxRates(1)], 1, 1);
            title('Action potential firing rate (Hz)')
            
            subplot(2,2,3)
            displayGrid(pLFPMatrix(:,:,realIndex), hot, [0, maxRates(2)], 1, 1);
            title('pLFP firing rate (Hz)')
            
            subplot(2,2,4)
            displayGrid(nLFPMatrix(:,:,realIndex), hot, [0, maxRates(3)], 1, 1);
            title('nLFP firing rate (Hz)')
                        
                        
        otherwise
            error('Invalid visualizeData string') 
    end
    
    % Save video frame if required
    if record == 1 && t>1
        writeVideo(vidObj, im2frame(zbuffer_cdata(fig)));
    end
    
    % Stop the animation if the stop button is pressed
    if stopFlag == 1
        break
    end
    
    drawnow
end
toc

% Close video file if required
if record == 1
    close(vidObj);
end
