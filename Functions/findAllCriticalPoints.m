function critpointStruct = findAllCriticalPoints(velocityX, velocityY)
% FINDALLCRITICALPOINTS finds the critical points at each point in time for
% the velocity fields determined by the 10x10xtime matrices VELOCITYX and
% VELOCITYY. CRITPOINTSTRUCT is a structure containing the type,
% location and Jacobian of every critical point detected in the velocity
% field vectors.
%
% FORMAT OF CRITPOINTSTRUCT
%        critpointStruct.(type).(statistic)
%   (type): saddle, stableFocus, stableNode, unstableFocus or unstableNode
%   (statistic): time, coords or jacobian
%       -time gives a Nx1 array of the integer time step of each point
%       -coords gives a Nx2 matrix of [row coordinate, column coordinate]
%       -jacobian gives a Nx2x2 matrix of the local Jacobian
% N is the number of critical points found

nSteps = size(velocityX, 3);

% Initialize output structure
types = {'saddle', 'stableFocus', 'stableNode', 'unstableFocus', ...
    'unstableNode'};
for jj = 1:length(types)
    jtype = types{jj};
    % If there are more than two critical points per time step in the
    % velocity field, the function will return an error. Change NSTEPS*2 in
    % the following statements to a higher value.
    critpointStruct.(jtype).time = nan(nSteps*2, 1);
    critpointStruct.(jtype).coords = nan(nSteps*2, 2);
    critpointStruct.(jtype).jacobian = nan(nSteps*2, 2, 2);
    critpointStruct.(jtype).currIndex = 1;
end

for istep = 1:nSteps
    % Find critical points at the current time step
    vx = velocityX(:,:,istep);
    vy = velocityY(:,:,istep);
    [rowcoords, colcoords, jacobians] = classifyCrit(vx, vy);
    
    % Store results only if there are any critical points
    if ~isempty(rowcoords)
        
        % Store each critical point in the structure
        for icrit = 1:size(rowcoords)
            
            % Sort pattern type by Jacobian trace and determinant
            ijac = jacobians(:,:,icrit);
            
            if det(ijac) < 0
                % Saddle point
                itype = 'saddle';
                
            elseif trace(ijac)^2 > 4*det(ijac)
                if trace(ijac) < 0
                    % Stable node
                    itype = 'stableNode';
                else
                    % Unstable node
                    itype = 'unstableNode';
                end
            else
                if trace(ijac) < 0
                    % Stable focus
                    itype = 'stableFocus';
                else
                    % Unstable focus
                    itype = 'unstableFocus';
                end
            end
            
            % Store pattern details
            iindex = critpointStruct.(itype).currIndex;
            critpointStruct.(itype).time(iindex) = istep;
            critpointStruct.(itype).coords(iindex,:) = ...
                [rowcoords(icrit), colcoords(icrit)];
            critpointStruct.(itype).jacobian(iindex,:,:) = ijac;
            
            % Increment storage array indices
            critpointStruct.(itype).currIndex = iindex + 1;
            
            % TODO: Add in code to double storage length if necessary
        end
    end
    
end

% Crop off empty rows of each matrix
for jj = 1:length(types)
    jtype = types{jj};
    lastIndex = critpointStruct.(jtype).currIndex;
    critpointStruct.(jtype).time(lastIndex:end) = [];
    critpointStruct.(jtype).coords(lastIndex:end,:) = [];
    critpointStruct.(jtype).jacobian(lastIndex:end,:,:) = [];
    critpointStruct.(jtype) = rmfield(critpointStruct.(jtype),'currIndex');
end