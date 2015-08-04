% FINDAVALANCHEDIST finds spiking avalanches is a CxT binary matrix (where
% C is the number of channels and T is the number of time steps), then
% plots the distribution of avalanche sizes and tries to fit a power-law to
% the distribution.

%% Set variables and load spiking matrix

% This script will act on the variable called SPIKEMATRIX in the current
% workspace by default. To load from a file instead, uncomment the line
% below and choose a file path.
%load('./my144-101_spikes.mat', 'spikeMatrix', 'Fs')

% SPIKEMATRIX will be cut into bins of BINSIZE time steps before avalanches
% are detected
binSize = 1;

% LOGHISTSCALE defines the size of each histogram bin for plotting in
% logarithmic scale
logHistScale = 0.2;

%% Find avalanches and plot as a log-log histogram
% Find avalanche sizes and lengths
[avSizes, avLengths] = findAvalanche(spikeMatrix, binSize);

% Find histogram edges and counts in logarithmic scale
edges = 10.^(0:logHistScale:log10(max(avSizes)));
h = histc(avSizes, edges);

% Plot and label
scatter(log10(edges), log10(h))
xlabel('Log10 of avalanche size')
ylabel('Log10 of number of avalanches')


%% Fit power law and plot
% Fit power-law parameters and plot a survivor plot of the avalanche data
% along with the line of best fit
[xmin, xmax, alpha] = pwrfit_minmax(avSizes, 1, max(avSizes)*0.7, 1);
fprintf('Best fit is alpha = %0.3f, xmin = %i, xmax = %i\n', alpha, xmin, xmax)

