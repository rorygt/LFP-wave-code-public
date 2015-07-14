function B = nanmedfilt2(A, varargin)
% NANMEDFILT2: 2-D median filtering with sparse NaN's.
% NaN's in A are replaced by the median value of all non-NaN values in A,
% then MATLAB's medfilt2 function is applied.

Anans = isnan(A);
A(Anans) = median(A(~Anans));

B = medfilt2(A, varargin{:});

B(Anans) = nan;