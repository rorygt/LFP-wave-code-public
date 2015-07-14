function index = windingNumber(vx, vy)
% WINDINGNUMBER computes the winding number of a sequence of vectors given
% by 1 dimensional VX and VY. VX and VY are assumed to be samples taken
% from a counterclockwise path in a vector field.

% Reject if VX and VY are different lengths or not 1D
if any(size(vx) ~= size(vy)) || numel(vx) ~= max(size(vx))
    error('Invalid size input vector VX or VY')
end

% Duplicate the first point so that the path starts and ends at the same
% spot, and only look at the sign
posvx = [vx(:); vx(1)] > 0;
posvy = [vy(:); vy(1)] > 0;

% Find where VY is positive for two steps in a row
posvy = posvy(1:(end-1)) & posvy(2:end);

% The index is the number of times V crosses North in a counterclockwise
% direction minus the number of crossings in a clockwise direction
index = sum(posvy & ( diff(posvx) == -1)) - ...
    sum(posvy & (  diff(posvx) == +1));