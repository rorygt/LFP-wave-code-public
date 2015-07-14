function varargout = anglegradientnan(f, nanIndices, angleFlag, singleIndex, varargin)
%  ANGLEGRADIENTNAN Approximate gradient of angular data f which may 
%     contain NaN values (which can have indices given in nanIndices to
%     save time).
%
%   Adapted from the standard MATLAB gradient function.
%   BUT custom spacing cannot be given and matrix F must be 1D or 2D.
%
%   Can also work for non-angular data, just set ANGLEFLAG = 0.
%
%   SINGLEINDEX is an optional integer argument that flags that the
%   gradient should only be calculated for a single entry in F with index
%   SINGLEINDEX.
%
%   [FX,FY] = GRADIENT(F) returns the numerical gradient of the
%   matrix F. FX corresponds to dF/dx, the differences in x (horizontal) 
%   direction. FY corresponds to dF/dy, the differences in y (vertical) 
%   direction. The spacing between points in each direction is assumed to 
%   be one. When F is a vector, DF = GRADIENT(F)is the 1-D gradient.
%
%   [FX,FY] = GRADIENT(F,H), where H is a scalar, uses H as the
%   spacing between points in each direction.
%
%
%   Note: The first output FX is always the gradient along the 2nd
%   dimension of F, going across columns.  The second output FY is always
%   the gradient along the 1st dimension of F, going across rows.  For the
%   third output FZ and the outputs that follow, the Nth output is the
%   gradient along the Nth dimension of F.
%

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 5.17.4.8 $  $Date: 2011/05/17 02:22:19 $

if nargin < 2
    nanIndices = find(isnan(f));
end

if nargin < 3
    angleFlag = 1;
end

[err,f,ndim,loc,rflag] = parse_inputs(f,varargin);
if err, error(message('MATLAB:gradient:InvalidInputs')); end

% Convert indices of channels to skip to matrix subscripts
if ndim == 1
    skipSubs = nanIndices;
elseif ndim == 2
    [skipx skipy] = ind2sub(size(f), nanIndices);
    skipSubs = [skipx(:)'; skipy(:)'];
else
    error('HighDim', 'Input matrix must have 2 or less dimensions.');
end


% Loop over each dimension. Permute so that the gradient is always taken along
% the columns.

if ndim == 1
  perm = [1 2];
else
  perm = [2:ndim 1]; % Cyclic permutation
end

for k = 1:ndim
   [n,p] = size(f);
   
   g  = zeros(size(f),class(f)); % case of singleton dimension
   
   % Take forward differences on left and right edges
   if n > 1
      g(1,:) = anglesubtract(f(2,:), f(1,:), angleFlag);
      g(n,:) = anglesubtract(f(n,:), f(n-1,:), angleFlag);
   end

   % Take centered differences on interior points
   if n > 2
      g(2:n-1,:) = anglesubtract(f(3:n,:), f(1:n-2,:), angleFlag);
   end
   
   % EXPERIMENTAL: Take 5 point stencil if possible
   if n > 3
       fivepoint = anglesubtract(...
           1/12 * anglesubtract(f(1:n-4,:), f(5:n,:), angleFlag), ...
           2/3 * anglesubtract(f(2:n-3,:), f(4:n-1,:), angleFlag));
       current = g(3:n-2, :);
       isValid = ~isnan(fivepoint);
       current(isValid) = fivepoint(isValid);
       g(3:n-2, :) = current;
   end
   
   % Take forward differences around NaN values.
   if k == 1
       xs = skipSubs(1,:);
       ys = skipSubs(2,:);
   else
       xs = skipSubs(2,:);
       ys = skipSubs(1,:);
   end
   
   % Behind NaN values
   xBehind = xs - 2;
   [gBehind isGoodBehind] = processnans(f, xBehind, ys, angleFlag);
   igBehind = sub2ind(size(g), xBehind(isGoodBehind)+1, ys(isGoodBehind));
   g(igBehind) = gBehind;
   
   % In front of NaN values
   xFront = xs + 1;
   [gFront isGoodFront] = processnans(f, xFront, ys, angleFlag);
   igFront = sub2ind(size(g), xFront(isGoodFront), ys(isGoodFront));
   g(igFront) = gFront;

   % Store results
   varargout{k} = ipermute(g,[k:max(ndim,2) 1:k-1]);

   % Set up for next pass through the loop
   f = permute(f,perm);
   
end 

% Swap 1 and 2 since x is the second dimension and y is the first.
if ndim>1
  tmp = varargout{1};
  varargout{1} = varargout{2};
  varargout{2} = tmp;
end

if rflag, varargout{1} = varargout{1}.'; end
end

%% Take forward differences around NaN values.
function [nangrad xgood] = processnans(f, xlower, y, angleFlag)
% PROCESSNANS
%   NANGRAD = PROCESSNANS(F, XLOWER, Y) returns the forward difference on
%   the matrix F between XLOWER and XLOWER+1 at columns Y.
    sf = size(f);
    xgood = xlower>0 & xlower<10;
    if ndims(f) == 1
        lowerCoords = xlower(xgood);
        higherCoords = xlower(xgood)+1;
    else
        lowerCoords = sub2ind(sf, xlower(xgood), y(xgood));
        higherCoords = sub2ind(sf, xlower(xgood)+1, y(xgood));
    end
    nangrad = anglesubtract(f(higherCoords), f(lowerCoords), angleFlag);
end


%% Parse inputs
function [err,f,ndim,loc,rflag] = parse_inputs(f,v)
%PARSE_INPUTS
%   [ERR,F,LOC,RFLAG] = PARSE_INPUTS(F,V) returns the spacing
%   LOC along the x,y,z,... directions and a row vector
%   flag RFLAG. ERR will be true if there is an error.

err = false;
loc = {};
nin = length(v)+1;

% Flag vector case and row vector case.
ndim = ndims(f);
vflag = 0; rflag = 0;
if iscolumn(f)
   ndim = 1; vflag = 1; 
elseif isrow(f) % Treat row vector as a column vector
   ndim = 1; vflag = 1; rflag = 1;
   f = f.';
end;
   
indx = size(f);

% Default step sizes: hx = hy = hz = 1
if nin == 1, % gradient(f)
   loc = cell(1, ndims(f));
   for k = 1:ndims(f)
      loc(k) = {1:indx(k)};
   end;

elseif (nin == 2) % gradient(f,h)
   % Expand scalar step size
   if (length(v{1})==1)
      loc = cell(1, ndims(f)); 
      for k = 1:ndims(f)
         h = v{1};
         loc(k) = {h*(1:indx(k))};
      end;
   % Check for vector case
   elseif vflag
      loc(1) = v(1);
   else
      err = true;
   end

elseif ndims(f) == numel(v), % gradient(f,hx,hy,hz,...)
   % Swap 1 and 2 since x is the second dimension and y is the first.
   loc = v;
   if ndim>1
     tmp = loc{1};
     loc{1} = loc{2};
     loc{2} = tmp;
   end

   % replace any scalar step-size with corresponding position vector
   for k = 1:ndims(f)
      if length(loc{k})==1
         loc{k} = loc{k}*(1:indx(k));
      end;
   end;

else
   err = true;
end
end
