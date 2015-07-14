function grad = singleanglegradientnan(f, index, angleFlag)
% SINGLEANGLEGRADIENTNAN computes the gradient of matrix F at index INDEX,
%   correctly adjusting for walls of the matrix and NAN values. ANGLEFLAG
%   determines whether to treat the data in F as angular or not; by default
%   it is set to 1.

if nargin < 3
    angleFlag = 1;
end

ndim = ndims(f);
grad = nan([length(index(:)), ndim]);

% Loop over every index required
for iindex = 1:length(index(:));

igrad = zeros(1, ndim);
perm = [2:ndim 1];

target = cell(1, ndim);
[target{:}] = ind2sub(size(f), index(iindex));

% Apply the appropriate formula to find gradient in both dimensions
for idim = 1:ndims(f); % Loop over all dimensions
    it = target{1};
    frow = f(:, target{2:end});
    
    if it+1 > length(frow) || isnan(frow(it+1))
        if it-1 < 1 || isnan(frow(it-1))
            % Node has an edge or nan on both sides, no valid gradient
            igrad(idim) = nan;
        else
            % Apply forward difference formula
            igrad(idim) = anglesubtract(frow(it), frow(it-1), angleFlag);
        end
    elseif it-1 < 1 || isnan(frow(it-1))
        % Apply forward difference formula
        igrad(idim) = anglesubtract(frow(it+1), frow(it), angleFlag);
    elseif it+2 > length(frow) || it-2 < 1 || isnan(frow(it+2)) || ...
            isnan(frow(it-2))
        % Apply centred difference formula
        igrad(idim) = anglesubtract(frow(it+1), frow(it-1), angleFlag);
    else
        % Apply 5 point stencil
        igrad(idim) = anglesubtract( ...
           1/12 * anglesubtract(frow(it-2), frow(it+2), angleFlag), ...
           2/3 * anglesubtract(frow(it-1), frow(it+1), angleFlag));
    end
    
    % Permute for next pass so that the idim dimension is the first index
    f = permute(f, perm);
    target = target(perm);
    
end % End loop over dimensions

grad(iindex, :) = igrad;

end% End loop over indices

% Swap first and second dimension so that 1st is x and 2nd is y
temp = grad(:,1);
grad(:,1) = grad(:,2);
grad(:,2) = temp;

% Reshape grad to be have the same dimensions as INDEX with an extra
% trailing dimension
grad = reshape(grad, [size(index) ndim]);
grad = squeeze(grad);

% Transpose if input was a column vector
if size(index, 1) == 1 && ndims(index) == 2
    grad = grad';
end


end