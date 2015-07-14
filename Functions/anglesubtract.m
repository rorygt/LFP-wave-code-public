function xydiff = anglesubtract(x, y, angleFlag)
% Subtracts the matrix y from the equally sized matrix x and returns a
% value between -pi and +pi. This is valid if x and y both contain angular
% data in radians. If angleFlag == 0, return a non-angular subtraction.

if nargin < 3
    angleFlag = 1;
end

if angleFlag == 1
    % METHOD 1: Modulo method, gives same results as method 2
    xydiff = mod(x - y + pi, 2*pi) - pi;
    
    
    % METHOD 2: Trig method, works but is slower
    % tic
    % xydiff2 = atan2(sin(x-y), cos(x-y));
    % t2 = toc;
    
    
    % METHOD 3: Old method, currently doesn't work and is also slower
    % % Convert matrices to vectors
    % xv = x(:);
    % yv = y(:);
    %
    % % Find the minimum absolute difference between x and y
    % differences = [xv-yv, xv+2*pi-yv, xv-yv-2*pi]';
    % [~, mindex] = min(abs(differences));
    % xydiff3 = differences(mindex);
    %
    % % Reshape to a matrix
    % xydiff3 = reshape(xydiff3, size(x));
    
else
    xydiff = x - y;
end

end