function [rowcoords, colcoords, jacobians] = classifyCrit(vx, vy)
% CLASSIFYCRIT finds and classifies the critical points in the vector field
%   defined by VX and VY. Outputs 2 1xN vectors (where N is the number of
%   critical points detected) expressing the row and  column coordinates of
%   each point, and a 2x2xN matrix expressing the estimated Jacobian at
%   each point.

% Find critical points
[rowcoords, colcoords] = criticalpointsbilinear(vx, vy);

jacobians = zeros(2,2,length(rowcoords));

for ic = 1:length(rowcoords)
    % Find partial derivatives of the 4 corners that the critical point
    % resides in
    ix = rowcoords(ic);
    iy = colcoords(ic);
    corners = sub2ind(size(vx), ...
        [floor(ix) floor(ix); ceil(ix) ceil(ix)], ...
        [floor(iy) ceil(iy); floor(iy) ceil(iy)]);
    cornersGradx = singleanglegradientnan(vx, corners, 0);
    cornersGrady = singleanglegradientnan(vy, corners, 0);
    
    % Estimate Jacobian at the critical point through bilinear
    % interpolation
    ixdec = ix - floor(ix);
    iydec = iy - floor(iy);
    dxx = [1-ixdec ixdec] * cornersGradx(:,:,1) * [1-iydec iydec]';
    dxy = [1-ixdec ixdec] * cornersGradx(:,:,2) * [1-iydec iydec]';
    dyx = [1-ixdec ixdec] * cornersGrady(:,:,1) * [1-iydec iydec]';
    dyy = [1-ixdec ixdec] * cornersGrady(:,:,2) * [1-iydec iydec]';
    jacobians(:,:,ic) = [dxx dxy; dyx dyy];
    
end
    
end