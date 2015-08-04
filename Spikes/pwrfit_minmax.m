function [xmin,xmax,alpha,KS] = pwrfit_minmax(x,xmin0,xmax0,plotFlag)

% Function which finds the parameters for the power-of-law-best-fit. It
% also optionally plots the fit against the data as a complementary
% cumulative distribution function. It is slow and inelegant but it gets
% the job done. It is a good idea to use a couple of starting guesses for
% xmin0, give or take a decade either way. This is because the optimization
% routine can sometimes converge to a local minimum rather than the global
% minimum.
% Note that the potential values for alpha are hard coded below: this
% version first looks for the best alpha in the range 1:0.1:3, then narrows
% this down to the nearest 0.005. A result of alpha=0.095 or =3.05
% indicates the best fit lies outside this range.
% INPUT
% x: Data to which we are attempting to apply a fit
% xmin0: User guess at what value the power-law starts
% xmax0: User guess at what value the power-law ends
% plotFlag: Optional binary flag indicating whether to plot the fit
%           (default 0)
% OUTPUT
% xmin: Lower scale parameter of power-law-of-best-fit
% xmax: Upper scale parameter of power-law-of-best-fit
% alpha: Shape parameter of power-law-of-best-fit
% EXAMPLE
% z = pwrrand(1e6,2.5,1) + exprnd(1,1,1e6) asymptotes to a power-law with
% alpha = 2.50 at high values but behaves like an exponential distribution
% at low values. Using xmin0 = 10 as a first guess we obtain xmin = 11.08,
% alpha = 2.63.


global dist

% Allows the data to be used by the sub-function ksfun
dist = x;

% Find the best value of xmin using xmin0 as a starting point - decrease
% the TolX value for finer precision (the minimum value of xmin that the
% optimization scheme can increment before it decides it has found the
% minimum)
xlims0 = [xmin0 xmax0];
options = optimset('TolX',1e-6);
[xlims, KS] = fminsearch(@ksfun,xlims0,options);

% Throw away data not between the scaling parameters xmin and xmax
x = x(x>=xlims(1) & x<=xlims(2));

% Adjust the scaling parameters so that they matches the smallest value in
% the remaining data
xmin = min(x);
xmax = max(x);

% Find rough estimate of shape parameter alpha
% NOTE: This is also hard coded into the KSFUN subfunction
% TODO: Allow ALPHA_RANGE to be set in the input variables and use it as a
% global variable in KSFUN.
alpha_range = 1:0.1:3;
alpha = alpha_mle(x, alpha_range);

% Find better estimate of shape parameter alpha
alpha_range = (alpha-0.05):0.005:(alpha+0.05);
alpha = alpha_mle(x, alpha_range);

% Calculate KS statistic again
[data_cdf,rng] = ecdf(x);
fit_cdf = cumsum( sum((xmin:xmax).^-alpha) * rng.^(-alpha)) ;
KS = max(abs(data_cdf - fit_cdf));

% Plot the data above xmin together with the best fit in a survival
% function (i.e. a complementary cumulative distribution function)
if nargin >= 4 && plotFlag ~= 0
    [data_ccdf,rng] = ecdf(x,'function','survivor');
    fit_ccdf = (rng/xmin).^(-alpha+1);
    loglog(rng,data_ccdf,'b.-',rng,fit_ccdf,'r');
    xlim([min(rng),1.05*max(rng)]);
    xlabel('x');
    ylabel('P(x<X)');
    legend('Sample distribution','Power-law-of-best-fit');
end

end

function KS = ksfun(xlims)

% Calculates the scaling parameter (alpha) and thereby calculates the
% Kolmogorov-Smirnov statistic for a given value of xmin and xmax. The
% optimization scheme calls this function and varies xmin and xmax until it
% finds a value which minimizes the Kolmogorov-Smirnov statistic.

global dist

% Throw away data outside the scaling thresholds indicated by xlims
x = dist(dist>=xlims(1) & dist<=xlims(2));
% Return a large value for KS if limits remove all data in x
if isempty(x)
    KS = 1e20;
    return
end

% Calculate the empirical cdf values (data_cdf) of our data and their range
% (rng)
[data_cdf,rng] = ecdf(x);

% Test a range of alpha values
alpha_range = 1:0.1:3;

% Find best estimate for alpha using maximum likelihood
[alpha, C] = alpha_mle(x, alpha_range);

% Calculate the resultant cumulative distribution function of this
% power-law-of-best fit
fit_cdf = cumsum(C * rng.^(-alpha)) ;

% Calculate the Kolmogorov-Smirnov distance
KS = max(abs(data_cdf-fit_cdf));

end

function [alpha, C] = alpha_mle(x, alpha_range)
% Calculates the maximum likelihood estimate of the shape parameter alpha
% of the power-law-of-best-fit to discrete data X. Also returns the 
% normalization constant C, where the power law is defined as
%       p(x) = C * x .^ -alpha,
% for min(x) <= x <= max(x). Only chooses the best value of alpha from the
% vector ALPHA_RANGE.

% Define the range of x
xrange = min(x):max(x);

% Calculate reciprocal of the normalization constant C for each alpha value
alpha_range = alpha_range(:)';
% This may run out of memory if xrange is big, so iteratively calculate C
% values if this happens
try
    rC = sum(repmat(xrange', 1, length(alpha_range)) .^ ...
        -repmat(alpha_range, length(xrange), 1));
catch
    rC = zeros(size(alpha_range));
    for ii = 1:length(alpha_range)
        rC(ii) = sum(xrange .^ -alpha_range(ii));
    end
end

% Find best value for alpha by maximum likelihood estimation
L = -length(x)*log(rC) - alpha_range*sum(log(x));
[~, I] = max(L);
alpha = alpha_range(I);
C = 1/rC(I);

end

