function [alpha, D] = leastSq(x, y, varargin)
%% [alpha, D] = leastSq(log(x), log(y), ax)
% Least squares estimator from Ling 2019, eq. 2.4. For regions where
% log(MSD) is expected to be linear with log(τ). If a third argument is
% provided, plot the fit (on axis ax if ax is non-empty, or gca otherwise)
% Note: This gives correct D for 1D case. For 2D ÷ 2, 3D ÷ 3
validateattributes(x,{'numeric'},{'column'})
validateattributes(y,{'numeric'},{'column'})

xbar = mean(x);
ybar = mean(y);

alpha = sum((y - ybar).*(x - xbar)) ./ sum( (x - xbar).^2 );
D = 0.5 * exp( ybar - alpha * xbar);

if nargin > 2
    if ~isempty(varargin{1})
        ax = varargin{1};
    else
        ax = gca;
    end
    plot(ax, exp(x), 2*D*exp(x).^alpha, 'k:','LineWidth',2)

end