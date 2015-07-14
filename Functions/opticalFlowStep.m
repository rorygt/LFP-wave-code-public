function [u, v, convergenceLoop] = opticalFlowStep(im1, im2, ...
    nanIndices, surroundLocs, weightFactors, alpha, beta,...
    showFlag, u0, v0, im0, im3, angleFlag)
% Find the optical flow U, V between two images IM1 and IM2 from a video
% sequence. 

% Maximum fractional change between iterations to be counted as a fixed
% point
maxChange = 0.01;

% Relaxation parameter for sucessive overrelaxation method
% relaxation = 1;

% Default inputs
if nargin<3
    nanIndices = find(isnan(im1));
end
if nargin<6
    alpha=1; % Smoothing factor
end
if nargin<7
     beta = 0.001; % Charbonnier penalty weighting
end
if nargin < 8
    showFlag = 0; % Do not show output
end
if nargin < 11
    u0 = zeros(size(im1(:,:,1))); % Use empty matrices as starting guesses
    v0 = zeros(size(im2(:,:,1)));
end

u = u0;
v = v0;

%% Find derivatives

% Spatial derivatives
[Ex1, Ey1] = anglegradientnan(im1, nanIndices, angleFlag);
[Ex2, Ey2] = anglegradientnan(im2, nanIndices, angleFlag);
Ex = (Ex1+Ex2)/2;
Ey = (Ey1+Ey2)/2;

% Temporal derivative
if nargin < 10 || isempty(im0) || isempty(im3)
    % Take centred difference
    Et = anglesubtract(im2, im1, angleFlag);
else
    % Use 5 point stencil
    Et = anglesubtract(1/12 * anglesubtract(im0, im3),...
        2/3 * anglesubtract(im1, im2), angleFlag);
end

% Define filter to estimate spatial mean
% meanFilt=[1/12 1/6 1/12;1/6 0 1/6;1/12 1/6 1/12];


%fixedPointFlag = 0;
dataP = inf(size(im1));
smoothP = dataP;

% Temporarily turn off warnings so that the singular matrix warning
warning off all

% Loop over different non-linear penalty functions until a fixed point is
% reached
for convergenceLoop = 1:1000
    
    % Compute the first order error in data and smoothness
    dataE = Ex.*u + Ey.*v + Et;
    [upx, upy] = anglegradientnan(u, nanIndices, 0);
    [vpx, vpy] = anglegradientnan(v, nanIndices, 0);
    smoothE = upx.^2 + upy.^2 + vpx.^2 + vpy.^2;
    
    % Compute nonlinear penalty functions
    lastDataP = dataP;
    dataP = (1 + dataE.^2 / beta^2) .^ (-1/2);
    
    lastSmoothP = smoothP;
    smoothP = (1 + smoothE / beta^2) .^ (-1/2);
    %smoothPAvg = conv2(smoothP, meanFilt, 'same');
    
    % Check if dataP and smoothP have reached a fixed point
    dataPChange = abs(dataP-lastDataP) ./ abs(dataP);
    smoothPChange = abs(smoothP-lastSmoothP) ./ abs(smoothP);
    
    if max(max(dataPChange)) < maxChange &&...
            max(max(smoothPChange)) < maxChange
%         display(convergenceLoop)
        break
    end
    
%     % OLD METHOD: Iteratively solve equations for u and v using successive
%     % over-relaxation
%     for ii = 1:niterations
%         
%         % Compute local averages of the flow components
%         uAvg = conv2(u, meanFilt, 'same');
%         vAvg = conv2(v, meanFilt, 'same');
%         
%         % Update estimation for flow vectors
%         u = (1-relaxation) * u + relaxation * ...
%             ( (smoothP + smoothPAvg)/2 .* uAvg - ...
%             1/alpha * Ex .* dataP .* (Ey .* v + Et) ) ./ ...
%             (1/alpha * dataP .* Ex.^2 + (smoothP + smoothPAvg)/2);
%         v = (1-relaxation) * v + relaxation * ...
%             ( (smoothP + smoothPAvg)/2 .* vAvg - ...
%             1/alpha * Ey .* dataP .* (Ex .* u + Et) ) ./ ...
%             (1/alpha * dataP .* Ey.^2 + (smoothP + smoothPAvg)/2);
%         
%         
%     end

    % NEW METHOD: manipulate the system of linear equations into standard
    % A*x=b form. x is a 200x1 vector conisting of u and then v.
    sigma = dataP ./ smoothP;
    delta = alpha + sigma.*Ex.^2 + sigma.*Ey.^2;
    N = size(im1, 1);
    A = zeros(2*N^2);
    b = zeros(2*N^2, 1);
    for irow = 1:N
        for jcol = 1:N
            % Define an index to turn row i and column j to a single number
            index = irow + N*(jcol-1);
            
            % Set row to zero if pixel is invalid
            if ismember(index, nanIndices)
                continue
            end
            
            % Initialize vectors uij and vij
            uij = zeros(1, 2*N^2);
            uij(index) = -delta(index);
            vij = zeros(1, 2*N^2);
            vij(index+N^2) = -delta(index);

            % Define index of surrounding spaces and weighting factor
            surr = surroundLocs(:, irow, jcol);
            factor = weightFactors(:, irow, jcol);
            
            % Fill in the matrix A
            % Values for uij
            uij(surr) = ( sigma(surr) .* Ey(surr).^2 + alpha) /4 .* factor;
            uij(surr+N^2) = - (sigma(surr).*Ex(surr).*Ey(surr)) ...
                /4 .* factor;
            A(index, :) = uij;
            
            % Values for vij
            vij(surr) = - (sigma(surr).*Ex(surr).*Ey(surr)) /4 .* factor;
            vij(surr+N^2) = ( sigma(surr) .* Ex(surr).^2 + alpha) /4 .* factor;
            A(index+N^2, :) = vij;
            
            % Fill in the vector b
            b(index) = sigma(index) * Ex(index) * Et(index);
            b(index+N^2) = sigma(index) .* Ey(index) .* Et(index);
        end
    end
    % Solve this system of linear equations
    A = sparse(A);
    xexact = A\b;
    
    % Reshape back to 10x10 grids
    u = reshape(xexact(1:N^2), N, N);
    v = reshape(xexact((1:N^2)+N^2), N, N);
%     quiver(u,v)
%     drawnow
    
    
end

% Turn warnings back on
warning on all

% Plot results
if showFlag == 1
    subplot(1,2,1)
    imagesc(im1)
    hold on
    quiver(u,v)
    hold off
    subplot(1,2,2)
    imagesc(im2)
    colorbar
end

end

    
