function [alpha, D] = leastSq(x, y)

validateattributes(x,{'numeric'},{'column'})
validateattributes(y,{'numeric'},{'column'})

xbar = mean(x);
ybar = mean(y);

alpha = sum((y - ybar).*(x - xbar)) ./ sum( (x - xbar).^2 );
D = 0.5 * exp( ybar - alpha * xbar);

end