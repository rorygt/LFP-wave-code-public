function [patternStats, pTypes] = printPatternStatistics(patterns, recordingLength)
% Prints pattern statistics for input PATTERNS. Patterns can either be a
% structure (produced by CREATECRITPOINTSTRUCTURE), an array of
% structures or a string containing a filename to open and extract a
% variable called PATTERNS from. RECORDINGLENGTH is an optional double
% that indicates the total length of recording in seconds, so that counts can be
% expressed as a rate rather than an absolute number.

% Process if input is a cell array of structures
cellFlag = 0;
if ischar(patterns)
    load(patterns, 'patterns', 'recordingLength')
elseif ~isequal(size(patterns), [1 1])
    cellFlag = 1;
    allPatterns = patterns(:);
    patterns = allPatterns{1};
end

nStructs = length(patterns);


pTypes = { 'synchrony', 'planeWave', 'sink', 'source', ...
    'spiralIn', 'spiralOut', 'saddle' };

%pTypes = unique(patterns.type)';
patternStats = zeros(3, length(pTypes));

% Loop over each pattern type
for iType = 1:length(pTypes)
    % Print pattern name
    disp(pTypes{iType})
    
    % Loop over each patterns structure and collate data
    iCounts = [];
    iDurations = [];
    iStepDisplacement = [];
    for jStruct = 1:nStructs
        if cellFlag
            patterns = allPatterns{jStruct};
        end
        if ~any(ismember(patterns.type, pTypes{iType}))
            disp('   No patterns of this type exist!')
            continue
        end
        jPattern = cropStructure(patterns, pTypes{iType});
        iCounts = [iCounts; length(jPattern.type)];
        iDurations = [iDurations; jPattern.duration];
        iStepDisplacement = [iStepDisplacement; jPattern.stepDisplacement];
    end
    
    % Print statistics
    if exist('recordingLength', 'var')
        fprintf('   Total number per second = %f\n', sum(iCounts)/recordingLength)
    else
        fprintf('   Total number = %i\n', sum(iCounts))
    end
    fprintf('   Mean duration = %f, SD = %f\n', mean(iDurations),...
        std(iDurations))
    fprintf('   Mean displacement = %f, SD = %f\n', ...
        mean(iStepDisplacement), std(iStepDisplacement))
    patternStats(:, iType) = [sum(iCounts); mean(iDurations);...
        mean(iStepDisplacement)];
    
end