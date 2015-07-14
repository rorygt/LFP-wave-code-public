function reducedCritpointStruct = reduceCritpointStruct(critpointStruct, ...
    minRadius, velocityX, velocityY, syncActive)
% REDUCECRITICALCARRAY crops structures in CRITICALCARRAY, produced by
% findallcriticalpoints.m, according to the parameters given. VELOCITYX and
% VELOCITYY define the velocity field that patterns are being detected in.
%
%  - MINRADIUS: critical points must extend for a radius of
%       MINRADIUS to be counted. This is evaluated for all points by
%       looking at the Poincare index taken around squares with radii from
%       1 to MONRADIUS. Additionally, stable nodes must lie in a circle of
%       radius MINRADIUS with negative divergence, unstable nodes must lie
%       in a region of positive divergence, and foci must lie in a region
%       of positive or negative curl.
%  - SYNCACTIVE: a binary array indicating when global synchronization is
%       active in the system. Any critical points active at the same time
%       as synchronization are deleted.

types = fieldnames(critpointStruct);

% Counters for the numbers of points discarded by each condition
edgeCount = zeros(size(types));
syncCount = edgeCount;
windingCount = edgeCount;
divCurlOrderCount = edgeCount;
originalCritCount = edgeCount;
reducedCritCount = edgeCount;

% Loop over all types of patterns in CRITPOINTSTRUCT
for ii = 1:length(types)
    
    % Extract the array of time steps and coordinates for the current type
    iType = types{ii};
    iTimestep = critpointStruct.(iType).time;
    iCoords = critpointStruct.(iType).coords;
    
    originalCritCount(ii) = length(iTimestep);
    
    %% Remove points when global synchrony is active
    if nargin >=5
        duringSync = syncActive(iTimestep);
    else
        duringSync = zeros(size(iTimestep));
    end
    
    %% Remove patterns too close to edges
    lowerLimit = 1 + minRadius;
    upperLimit = 10 - minRadius;
    goodRows = iCoords(:,1) > lowerLimit & iCoords(:,1) < upperLimit;
    goodCols = iCoords(:,2) > lowerLimit & iCoords(:,2) < upperLimit;
    goodPosition = goodRows & goodCols;
    
    %% Remove patterns with inappropriate winding number or div/curl
    nPoints = length(iTimestep);
    discardWinding = zeros(nPoints, 1);
    discardDivCurl = zeros(nPoints, 1);
    
    for jp = 1:nPoints
        % Skip if point is already identified as invalid
        if ~goodPosition(jp) || duringSync(jp)
            continue
        end
        
        % Find winding number of squares around pattern centre
        jvx = velocityX(:,:,iTimestep(jp));
        jvy = velocityY(:,:,iTimestep(jp));
        wNum = zeros(1, minRadius);
        for irad = 1:minRadius
            sqrIndices = counterclkSqr(iCoords(jp,1), iCoords(jp,2), irad);
            wNum(irad) = windingNumber(jvx(sqrIndices), jvy(sqrIndices));
        end
        
        switch iType
            case 'saddle'
                % Pattern is a saddle, winding number must be negative
                discardWinding(jp) = any(wNum > -1);
                
            case {'stableNode', 'unstableNode'}
                % Pattern is a stable or unstable node, winding number must be
                % positive
                discardWinding(jp) = any(wNum < 1);
                
                % Calculate divergence of the velocity field
                jDiv = divergence(jvx, jvy);
                
                % Flip sign of divergence if pattern is a stable node
                if strcmp(iType, 'stableNode')
                    jDiv = -jDiv;
                end
                
%                 % TEMP: Remove divergence and curl measures
%                 % Mark pattern to be discarded if it is not in a positive
%                 % region
%                 discardDivCurl(jp) = ~positiveArea(jDiv, ...
%                     iCoords(jp,1), iCoords(jp,2), minRadius);
%                 

            case {'stableFocus', 'unstableFocus'}
                % Pattern is a stable or unstable focus, winding number
                % must be positive
                discardWinding(jp) = any(wNum < 1);
                
                % Calculate curl of the velocity field
                jCurl = curl(jvx, jvy);
                
%                 % TEMP: Remove divergence and curl measures
%                 % Mark pattern to be discarded if it is not in a positive
%                 % or negative region of curl/vorticity
%                 isPositive = positiveArea(jCurl, iCoords(jp,1),...
%                     iCoords(jp,2), minRadius);
%                 isNegative = positiveArea(-jCurl, iCoords(jp,1),...
%                     iCoords(jp,2), minRadius);
%                 discardDivCurl(jp) = ~isPositive & ~isNegative;

        end
        
    end
    
    %% Remove invalid points
    reducedCritpointStruct.(iType) = selectStructEntries(critpointStruct.(iType), ...
        goodPosition & ~duringSync & ~discardDivCurl & ~discardWinding);
    
    % Count number of points removed by each criteria
    edgeCount(ii) = sum(~goodPosition);
    syncCount(ii) = sum(duringSync);
    windingCount(ii) = sum(discardWinding);
    divCurlOrderCount(ii) = sum(discardDivCurl);
    reducedCritCount(ii) = sum(goodPosition & ~duringSync & ...
        ~discardDivCurl & ~discardWinding);
    
end
    
    % OPTIONAL: Show the number of critical points discarded
    disp('Critical point structure reduced, summary table below.')
    disp('Columns: saddles, stable nodes, unstable nodes,')
    disp('stable foci, unstable foci.')
    disp('Rows: initial number of points,') 
    disp('too close to edge, during synchrony,')
    disp('winding number invalid, div/curl invalid,')
    disp('final number.')
    disp([originalCritCount, edgeCount, syncCount, windingCount, ...
        divCurlOrderCount, reducedCritCount]')
end

    function isPositive = positiveArea(grid, cRow, cCol, radius)
        % Searches 2D array GRID for an area of radius RADIUS centred at (CROW,
        % CCOL) that contains only positive values. ISPOSITIVE = 1 if the area
        % is all positive, = 0 otherwise.
        
        isPositive = 1;
        
        % Loop over all cells nearby
        for irow = ceil(cRow - radius) : floor(cRow + radius)
            for icol = ceil(cCol - radius) : floor(cCol + radius)
                % If cell is close enough centre and is not positive, change
                % flag and return result
                if (irow - cRow) ^ 2 + (icol - cCol) ^ 2 <= radius && ...
                        grid(irow, icol) <= 0
                    isPositive = 0;
                    return
                end
            end
        end
    end
    
    function edgeIndices = counterclkSqr(row, col, radius)
    % Find the indices in a 10x10 grid of a square of MINRADIUS around a given
    % location with points given in a counterclockwise direction, starting from
    % the top left corner.
    sqrSide = 2*radius;
    edgeIndices = sub2ind([10, 10], ... % Convert to index
        [1:sqrSide, sqrSide * ones(1, sqrSide-2), ... % Rows
        sqrSide:-1:1, ones(1, sqrSide-2)] + ceil(row-radius-1), ...
        [ones(1, sqrSide), 2:(sqrSide-1), ... % Columns
        sqrSide * ones(1, sqrSide), (sqrSide-1):-1:2] + ceil(col-radius-1));
    end