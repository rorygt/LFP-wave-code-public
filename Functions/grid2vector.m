function vector = grid2vector(grid)
% Converts a 10x10xt grid into a 100xt vector, rotated correctly.
% For details on the transformation, see the related function vector2grid.

[nrow, ncol, t] = size(grid);

% Check that grid is the right size
if nrow ~= 10 || ncol ~= 10
    error('Grid must be 10x10xt.')
end

% Rotate grid 180 degrees and then reshape to a vector
if t == 1
    vector = reshape(rot90(grid, 2), 100, 1);
else
    vector = zeros(100, t);
    for it = 1:t
        vector(:,it) = reshape(rot90(grid(:,:,it), 2), 100, 1);
    end
end

end