function [velocityX, velocityY, allConvSteps] = ...
    opticalFlow(video, badChannels, alpha, beta, nStepsDisplay)
% OPTICALFLOW estimates the velocity of pixels between frames in the 
% (x x y x t) matrix VIDEO. 
%
% INPUTS:
%   ALPHA is the smoothness weighting. 
%   BADCHANNELS is a vector containing the indices of all channels to be
%       ignored in VIDEO.
%   NSTEPSDISPLAY is an optional integer specifying the number of steps
%       before a message is displayed. If NSTEPSDISPLAY is not defined or
%       is set to 0, no messages will be displayed, otherwise the progress
%       of the calculation will be updated every NSTEPSDISPLAY steps.

%TEMP
angleFlag = 0;

% Initialize result structures
nrows = size(video, 1);
ncols = size(video, 2);
nframes = size(video, 3);
ivx = zeros(nrows, ncols);
ivy = ivx;
velocityX = zeros(nrows, ncols, nframes-1);
velocityY = velocityX;
allConvSteps = nan(1, nframes-1);

% Initialize temporary variables
prevFrame = [];
frame = video(:,:,1);
nextFrame = video(:,:,2);
next2Frame = video(:,:,3);

% To construct the coefficient matrix A for the system of linear equations
% in OPTICALFLOWSTEP, it is necessary to be able to find the location of
% the four surrounding pixels around any given pixel in index form. Set
% this matrix up now to save computational time in the loop
surroundLocs = zeros(4, nrows, ncols);
weightFactors = surroundLocs+1;
% Select arbitrary valid point to use as default index
defInd = min(setdiff(1:(nrows*ncols), badChannels));
for irow = 1:nrows
    for icol = 1:ncols
        % This is basically the SUB2IND function, but allows subscripts
        % outside the edge of the array
        isurroundLocs = [irow-1; irow+1; irow; irow] + ...
            ncols * ([icol; icol; icol-1; icol+1] - 1);
        % When surrounding locations are off the edge of the array or are a
        % location in BADCHANNELS, duplicate the entry in the other
        % direction. If this occurs, change the WEIGHTFACTOR from 1 to 2 to
        % indicate it has been duplicated.
        if irow == 1 || ismember(isurroundLocs(1), badChannels)
            isurroundLocs(1) = isurroundLocs(2);
            weightFactors(1:2,irow,icol) = 2;
        end
        if irow == nrows || ismember(isurroundLocs(2), badChannels)
            % If there is an invalid row on both sides, use averaging from
            % columns only
            if  isurroundLocs(1) == isurroundLocs(2)
                isurroundLocs(1) = defInd;
                isurroundLocs(2) = defInd;
                weightFactors(1:2,irow,icol) = 0;
                weightFactors(3:4,irow,icol) = 2;
            else
                isurroundLocs(2) = isurroundLocs(1);
                weightFactors(1:2,irow,icol) = 2;
            end
        end
        if icol == 1 || ismember(isurroundLocs(3), badChannels)
            isurroundLocs(3) = isurroundLocs(4);
            weightFactors(3:4,irow,icol) = 2*weightFactors(3:4,irow,icol);
        end
        if icol == ncols || ismember(isurroundLocs(4), badChannels)         
            % If there is an invalid column on both sides, use averaging from
            % rows only
            if isurroundLocs(3) == isurroundLocs(4)
                isurroundLocs(3) = defInd;
                isurroundLocs(4) = defInd;
                weightFactors(3:4,irow,icol) = 0;
                weightFactors(1:2,irow,icol) = 2*weightFactors(1:2,irow,icol);
            else
                isurroundLocs(4) = isurroundLocs(3);
                weightFactors(3:4,irow,icol) = 2*weightFactors(3:4,irow,icol);
            end
        end
        surroundLocs(:,irow,icol) = isurroundLocs;
        
    end
end

% Loop over all time steps
for it = 4 : ( size(video, 3) + 2 )
    
    % Calculate optical flow
    [ivx, ivy, convSteps] = opticalFlowStep(frame, nextFrame, ...
        badChannels, surroundLocs, weightFactors, alpha, ...
        beta, 0, ivx, ivy, prevFrame, next2Frame, angleFlag);
    
    % Store results
    allConvSteps(it-3) = convSteps;
    velocityX(:,:,it-3) = ivx;
    velocityY(:,:,it-3) = ivy;
    
    % Display the current step every MSTEPSDISPLAY steps
    if nargin >=5 && nStepsDisplay > 0
        if mod(it, nStepsDisplay) == 0
            display(sprintf('Calculating velocity, step %d', it))
        end
    end
    
    % Next set of frames
    prevFrame = frame;
    frame = nextFrame;
    nextFrame = next2Frame;
    if it <= size(video, 3)
        next2Frame = video(:, :, it);
    else
        next2Frame = [];
    end
    
end