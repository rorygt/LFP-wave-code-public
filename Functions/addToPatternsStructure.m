function [pStruct, reducedActiveTimes] = addToPatternsStructure(...
    typeName, activeTimes, positions, params, exisitingPStruct)
% Function ADDTOPATTERNSSTRUCTURE processes critical point times and
% locations into a readable structure of patterns.
%
% Inputs:
%   TYPENAME: String of the name of the type of pattern being stored
%   ACTIVETIMES: Either a sorted 1D vector of indices or a binary vector
%       indicating when the critical point is active.
%   POSITIONS: An Nx2 matrix indicating the row and column coordinates of
%       each critical point detected (N is the number of critical points,
%       N=length(ACTIVETIMES>0) ). If adding a pattern without a location
%       (e.g. plane wave or synchrony), this should be an empty array.
%   PARAMS: Structure of parameters:
%       -params.minDuration, the minimum duration (in time steps) for 
%           a pattern to be recorded.
%       -params.maxTimeGap, the maximum time gap (in time steps) between
%           critical points being present for them to be counted as the 
%           same pattern.
%       -params.maxDisplacement, the maximum spatial displacement (in grid
%           units) between consecutive critical points for them to be
%           counted as the same pattern.
%   EXISTINGPSTRUCT: Optional argument, if this is input then output
%       PSTRUCT will be appended to the existing structure.
%
%   Outputs:
%       PSTRUCT: A structure comprised of the following fields, where M is
%       the number of patterns detected:
%           - pStruct.type, Mx1 cell array of the string TYPENAME.
%           - pStruct.startTime, Mx1 array of the starting time of each
%               pattern (in time steps).
%           - pStruct.endTime, Mx1 array of the end time of each pattern.
%           - pStruct.duration, Mx1 array of the duration of each pattern.
%           - pStruct.startPosition, Mx2 array of the position of the
%               critical point that starts the pattern ([row, column]).
%           - pStruct.stepDisplacment, Mx1 array of the average
%               displacement (in grid units) between successive critical
%               points.
%           - pStruct.MSD, Mx1 array of the mean squared displacement from
%               the starting position.
%       REDUCEDACTIVETIMES: Shows which times in ACTIVETIMES were actually
%       used. If ACTIVETIMES is a binary array, REDUCEDACTIVETIMES is a
%       binary array of the same size with only 1's that were used in
%       patterns remaining. If ACTIVETIMES is a vector of indices,
%       REDUCEDACTIVETIMES is a binary array of the same size with 1's
%       indicating that that particular index was used.
%   TODO: Currently REDUCEDACTIVETIMES is always a full length
%   binary array even if 



%% Case where positions are not supplied (not critical point pattern)
if isempty(positions)  % Pattern type is not a critical point pattern
    wasIndexInput = 0;
    if max(activeTimes) > 1 % ACTIVETIMES is not a binary vector
        % Convert ACTIVETIMES to a binary vector
        wasIndexInput = 1;
        binaryTemp = zeros(max(activeTimes), 1);
        binaryTemp(activeTimes) = 1;
        activeTimes = binaryTemp;
    end
    
    % Skip 
    
    % Find all runs at least MINDURATION steps long, ignoring gaps of
    % MAXTIMEGAP
    [startTimes, endTimes] = findRuns(activeTimes, params.minDuration, ...
        0, params.maxTimeGap);
    
    % Create REDUCEDACTIVETIMES by storing only times when an occurance of
    % a pattern is recorded
    reducedActiveTimes = zeros(size(activeTimes));
    if ~isempty(startTimes)
        for iocc = 1:length(startTimes)
            reducedActiveTimes(startTimes(iocc) : endTimes(iocc)) = 1;
        end
    end
    
    % Initialize structure entries
    startPositions = nan(length(startTimes), 2);
    stepDisplacements = nan(size(startTimes));
    MSDs = stepDisplacements;
    

else
    %% Case where positions are supplied (critical point pattern)
    % Make sure that ACTIVETIMES is in index format
    wasIndexInput = 1;
    if max(activeTimes) <= 1
        wasIndexInput = 0;
        activeTimes = find(activeTimes == 1);
    end
    
    % Initialize binary array indicating if a critical point has already
    % been rejected or added to a different pattern
    pointIsSearched = zeros(size(activeTimes));
    % Initialize binary array indicating if critical points have actually
    % been used to make a pattern
    pointIsUsed = pointIsSearched;
    % Initialize output vectors
    ioutput = 1;
    % Note: If every point is used in a pattern, there would be this number
    % of patterns detected       v       v          v       v
    startTimes = zeros(ceil(length(activeTimes)/params.minDuration), 1);
    endTimes = startTimes;
    startPositions = repmat(startTimes, 1, 2);
    stepDisplacements = startTimes;
    MSDs = startTimes;
    
    % Loop over each critical point
    for icrit = 1:length(activeTimes)
        % Recursively find a chain of critical points starting from the
        % current point using the FINDCRITPATTERNS subfunction below
        pointsInPattern = findCritPatterns(icrit, activeTimes, ...
            positions, params, pointIsSearched);
        
        % Mark points as having been searched
        pointIsSearched(pointsInPattern) = 1;
        
        % Store pattern only if it lasts long enough
        if length(pointsInPattern) > params.minDuration;
            % Store results in output vectors
            startTimes(ioutput) = activeTimes(pointsInPattern(1));
            endTimes(ioutput) = activeTimes(pointsInPattern(end));
            iposition = positions(pointsInPattern(1), :);
            startPositions(ioutput,:) = iposition;
            stepDisplacements(ioutput) = mean( findDist(...
                positions(pointsInPattern(1:(end-1))), ...
                positions(pointsInPattern(2:end)) ));
            MSDs(ioutput) = mean( findDist( positions(pointsInPattern,:), ...
                repmat(iposition, length(pointsInPattern), 1)) );
            
            % Increment output index
            ioutput = ioutput + 1;
            % Mark points as actually used
            pointIsUsed(pointsInPattern) = 1;
            
        end
        
    end
    
    % Crop off zero entries in output vectors
    outSection = 1:(ioutput-1);
    startTimes = startTimes(outSection);
    endTimes = endTimes(outSection);
    startPositions = startPositions(outSection, :);
    stepDisplacements = stepDisplacements(outSection);
    MSDs = MSDs(outSection);
    
    % Save POINTISUSED binary array for output
    reducedActiveTimes = pointIsUsed;
    
end

%% Combine statistics into a structure array
% Create new structure or use existing structure
if nargin < 5
    oldP = struct('type', [], 'startTime', [], 'endTime', [], ...
        'duration', [], 'startPosition', [], 'stepDisplacement', [], ...
        'MSD', []);
else
    % Shorten variable name to make typing up next section more easy!
    oldP = exisitingPStruct;
end

% Create vector of a string of the type of the current pattern
types = repmat({typeName}, length(startTimes), 1);

% Store results as fields in PSTRUCT
pStruct.type = [oldP.type; types];
pStruct.startTime = [oldP.startTime; startTimes];
pStruct.endTime = [oldP.endTime; endTimes];
pStruct.duration = [oldP.duration; endTimes-startTimes+1];
pStruct.startPosition = [oldP.startPosition; startPositions];
pStruct.stepDisplacement = [oldP.stepDisplacement; stepDisplacements];
pStruct.MSD = [oldP.MSD; MSDs];

% Make sure that REDUCEDACTIVETIMES is a binary array
reducedActiveTimes = reducedActiveTimes > 0;

end


function pointsInPattern = findCritPatterns(thisIndex, activeTimes, ...
    positions, params, pointIsSearched)
% Subfunction to recursively find a chain of critical points starting at
% the point indicated by THISINDEX that are separated in time by less than
% PARAMS.MAXTIMEGAP and in space by less than PARAMS.MAXDISPLACEMENT

thisTime = activeTimes(thisIndex);
thisPos = positions(thisIndex, :);

% Find the critical point in the next PARAMS.MAXTIMEGAP time steps that
% occurs the soonest after THISTIME and has displacement less than
% PARAMS.MAXDISPLACEMENT.
foundIndex = 0;
lastTime = thisTime;
smallestDisp = inf;

% Iteratively check the next points
for ipoint = thisIndex : length(activeTimes)
    
    % Find the time difference
    itimeDiff = activeTimes(ipoint) - thisTime;
    % Stop looking if the time difference is too great or if a point has
    % already been found in a previous time step
    if itimeDiff > params.maxTimeGap || ...
            (foundIndex>0 && (activeTimes(ipoint) - lastTime) > 0)
        break
    end
    
    % Check if candidate has not been used and is not in the same step
    if ~pointIsSearched(ipoint) && itimeDiff > 0
        % Find displacement
        idistance = findDist(thisPos, positions(ipoint, :));
        % Check if pattern is closer than the maximum distance and closer
        % than any other found patterns
        if idistance < params.maxDisplacement && idistance < smallestDisp
            foundIndex = ipoint;
            smallestDisp = idistance;
        end
    end
    
    lastTime = activeTimes(ipoint);
end

% If no pattern is found, return just the input index
if foundIndex == 0
    pointsInPattern = thisIndex;
else
    % Recursively continue the chain from point FOUNDINDEX
    furtherPoints = findCritPatterns(foundIndex, activeTimes, ...
        positions, params, pointIsSearched);
    pointsInPattern = [thisIndex; furtherPoints];
end

end

function d = findDist(pos1, pos2)
% Subfunction to find the Euclidean distance between positions POS1 and
% POS2, where POS1 and POS2 are 1x2 vectors in the form of [row, column].
% If POS1 and POS2 are Nx2 vectors, distance is calculated between each
% row.
d = sqrt(sum((pos1-pos2).^2, 2));

end
