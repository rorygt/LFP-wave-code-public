function [hgrid, hcb] = ...
    displayGrid(signal, colorMapSpec, signalLimits, barFlag, subfigFlag)

% Takes a n^2 channel signal (as a n^2 vector or nxn grid) and plots it
% as a nxn grid with specified colormap. colorMapSpec provides the input
% to the colormap function, signalLimits is a 1x2 array giving the limits
% of the signal [min max], barFlag is a boolean determining whether to plot
% a colorbar or not, subfigFlag is a boolean determining whether or not
% the current figure is a subfigure (and should have colours frozen) or a
% single axis and firstFlag is a boolean establishing if it is the first
% time the grid is being displayed (in which case axes, grid lines etc need
% to be established).

[nrow, ncol] = size(signal);

% Use default scaled limits if absolutes are not specified
if nargin < 3 || isempty(signalLimits)
    defaultScale = 1;
else
    defaultScale = 0;
end

if nrow == 1 && mod(ncol, sqrt(ncol)) < 1e-10 % Signal is 1xn^2 row vector
    % Convert vector to grid
    n = round(sqrt(ncol));
    signalGrid = vector2grid(signal);
elseif ncol == 1 && mod(nrow, sqrt(nrow)) < 1e-10 % Signal is column vector
    % Convert vector to grid
    n = round(sqrt(ncol));
    signalGrid = vector2grid(signal);
elseif nrow == ncol % Signal is n x n grid
    n = nrow;
    signalGrid = signal;
else
    error('signalSizeDisplay', 'That size signal cannot be displayed!')
end


% Display signal
offset = 0.5 / (n/10);
if defaultScale == 1
    hgrid = imagesc([offset, 10-offset], [offset, 10-offset], signalGrid);
else
    hgrid = imagesc([offset, 10-offset], [offset, 10-offset], signalGrid, signalLimits);
end
colormap(colorMapSpec)

% Freeze colormap
if nargin >= 5 && subfigFlag == 1
    freezeColors
end

% Add grid lines
hold on
for kk = 1:10
    plot([kk kk],[0 10],'k-')
    plot([0 10],[kk kk ],'k-')
end
axis square
hold off

% Display and freeze colorbar if required
if nargin >= 4 && barFlag == 1
    hcb = colorbar;
    % Change labelling if signal is phase data from -pi to +pi
    if isequal(signalLimits, [-pi pi])
        set(hcb, 'YTick', [-pi pi], 'YTickLabel', {'-p', '+p'}, ...
            'FontName', 'symbol')
    end
    if nargin >= 5 && subfigFlag == 1
        cbfreeze(hcb)
    end
end

% Label figure
set(gca,'Xtick', 0.5:9.5);
set(gca,'Ytick', 0.5:9.5);
set(gca,'XTickLabel',num2str((90:-10:0).'),'tickLength', [0 0]);
set(gca,'YTickLabel',num2str((10:-1:1).'),'tickLength', [0 0]);


end