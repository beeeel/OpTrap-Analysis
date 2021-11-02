function [alpha, D] = leastSq(x, y)
% Least squares estimator from Ling 2019, eq. 2.4. For regions where
% log(MSD) is expected to be linear with log(τ).
validateattributes(x,{'numeric'},{'column'})
validateattributes(y,{'numeric'},{'column'})

xbar = mean(x);
ybar = mean(y);

alpha = sum((y - ybar).*(x - xbar)) ./ sum( (x - xbar).^2 );
D = 0.5 * exp( ybar - alpha * xbar);

end