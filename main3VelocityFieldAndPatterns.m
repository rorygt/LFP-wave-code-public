% Script to compute velocity field from a phase map, find all critical
% points in the velocity field and find patterns from critical points and
% the velocity field.

%% Set parameters
% Parameters should be set by the MAIN_findWaves.m file, or uncomment the
% following lines to set parameters to run this script independently.
% 
% % OPTICAL FLOW PARAMETERS
% % Smoothness weighting parameter: alpha > 0, higher alpha gives a smoother
% % velocity field
% alpha = 20;
% % Charbonnier penalty function scale
% beta = 0.001;
% 
% % SYNCHRONIZATION AND TRAVELLING WAVE PARAMETERS
% % Minimum order parameter for a plane wave to be detected
% planeWaveThreshold = 0.85;
% % Maximum velocity field magnitude for synchrony to be detected
% % If set to zero, default is mean - 1*std over all time
% synchronyThreshold = 0;
% 
% % PATTERN DETECTION PARAMETERS
% % Minimum duration of a pattern for it to be stored (in seconds)
% minDurationSecs = 0.01;
% % Maximum duration between critical points (or synchrony/plane
% % waves) for them to be counted as the same pattern (in seconds)
% maxTimeGapSecs = 0.005;
% % Maxiumum displacement between critical points between time steps for them
% % to be counted as the same pattern (measured in grid spaces)
% maxDisplacement = 1;
% % Minimum spatial radius for a critical point to occupy for it to be
% % counted, quantified by the winding number and dicergence/curl (in grid
% % spaces)
% minCritRadius = 2;
% 
% nStepsDisplay = 50000;


%%
% Use data that has already been loaded if it exists, otherwise choose the
% exeriment and file to load
if ~exist('experiment', 'var') || ~exist('file', 'var') || ...
        ~exist('phase', 'var') || ~exist('Fs', 'var') || ...
        ~exist('badChannels', 'var')
    experiment = 'MY147';
    file = 53;
    fLow = 1;
    fHigh = 4;
    sprintf('Loading file %s-%d to find patterns', experiment, file)
    load(sprintf('filteredLFPsHilbert_%d-%dHz_%s-%d.mat', ...
        fLow, fHigh, experiment, file))
end

% Set parameters structure
params.minDuration = fix(minDurationSecs*Fs);
params.maxTimeGap = fix(maxTimeGapSecs*Fs);
params.maxDisplacement = maxDisplacement;

% Change 100xT phase matrix to 10x10xT and define the index of the
% electrodes to be ignored
if size(phase, 1) == 100
    phase = vector2grid(phase);
end
badChannelsGrid = 101 - badChannels;

% Compute phase velocity field from phase map using optical flow
display('Computing phase velocity field...');
tic
[velocityX, velocityY, allConvSteps] = opticalFlow(phase, badChannelsGrid, alpha,...
    beta, nStepsDisplay);
toc
clearvars phase

% Find synchronized activity and travelling waves
display('Finding plane waves and synchrony...')
tic
[phi, v0, v_direction] = orderParameter(velocityX, velocityY);
% If threshold for synchronization is not defined, use the one standard
% deviation under the mean TODO: check this
if synchronyThreshold == 0
    synchronyThreshold = mean(v0) - std(v0);
end
% Find all plane waves
[patterns, pwActive] = addToPatternsStructure(...
    'planeWave', phi>=planeWaveThreshold, [], params);
% Find all periods of synchronous activity
[patterns, syActive] = addToPatternsStructure(...
    'synchrony', v0<=synchronyThreshold, [], params, patterns);
toc


% Detect all critical points
display('Finding critical points...')
tic
critpointStruct = findAllCriticalPoints(velocityX, velocityY);
toc

% Save all variables and delete large ones to save RAM
save(outputFile3Name)

% Remove invalid critical points
%
% Patterns are invalid if they are too close to the edge (specified by
% MINRADIUS), if they occur at the same time as there is global synchrony
% (specified by SYNCHRONYTHRESHOLD) or if they are not spatially extended
% as defined by the winding number or the divergence for nodes and curl for
% foci.
display('Removing invalid critical points...')
tic
reducedCritpointStruct = reduceCritpointStruct(critpointStruct, ...
    minCritRadius, velocityX, -velocityY, syActive);
toc

clearvars velocityX velocityY 

% Find critical point patterns that persist over time
display('Finding critical point patterns...')
tic
critPTypes = {'stableNode', 'unstableNode', 'stableFocus', ...
    'unstableFocus', 'saddle'};
patternTypes = {'sink', 'source', 'spiralIn', 'spiralOut', 'saddle'};
for itype = 1:length(critPTypes)
    istruct = reducedCritpointStruct.(critPTypes{itype});
    [patterns, iactive] = addToPatternsStructure(patternTypes{itype}, ...
        istruct.time, istruct.coords, params, patterns);
    patternCentreStruct.(patternTypes{itype}) = selectStructEntries(...
        istruct, iactive);
end
toc

% Save all results
save(outputFile3Name, 'reducedCritpointStruct', 'patterns', ...
    'critPTypes', 'patternTypes', 'patternCentreStruct', '-append')
fprintf('Saved to file %s\n', outputFile3Name)