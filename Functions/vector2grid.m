function grid = vector2grid(vector)
% Converts a 1x100, 100x1 or 100xt vector into a 10x10 or 10x10xt grid,
% rotated correctly. channelDim optionally specifies the dimension
% specifying channel number, which must have length 100.
%
% NOTE: THIS TRANSFORMATION CHANGES THE INDEXING OF THE DATA
% A row vector is stored as [1 2 3 4 ... 99 100]. However, in a 10x10 grid
% format, Paul's group stores the configuration of electrodes as:
%       [100  ......  10 ]
%       | .            . |
%       | .            . |
%       | .            2 |
%       [91   ......   1 ]
% This is 180 degrees rotated from the default MATLAB reshaping. The
% relationship between the index of a value is a row or column vector IV
% and the index of the same value in this rotated matrix IM is:
%    IV = 101 - IM

[nx, ny] = size(vector);

% Check which dimension has length 100 (preferentially choosing the first
% dimension)
if nx == 100
    t = ny;
elseif ny == 100
    t = nx;
    vector = vector';
elseif nx == 10 && ny == 10
    return
else
    error('vectorLength', 'Invalid input vector size!')
end


% Reshape to 10x10 grid and rotate
if t==1
    grid = rot90(reshape(vector,10,10), 2);
    
else
    % Reshape to 10x10 grid at every time step
    grid = zeros(10,10,t);
    for it = 1:t
        grid(:,:,it) = rot90(reshape(vector(:,it),10,10), 2);
    end

end