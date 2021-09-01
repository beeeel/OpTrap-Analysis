% Numerical derivative, need to find reference for 5 point stencil
function dxdt = stencil5log(tau, x)
%% Compute the numerical derivative of each column of x for logspaced data
% dxdt = stencil5log(tau, x)
% If x is a function sampled at points defined in tau, and tau is
% logarithmically spaced ( tau(n+1) = a * tau(n) ), this alternate 5-point
% stencil will accurately calculate the gradient.

% Lazy hack for array shapes
if ~(size(tau,1) == size(x, 1) && size(tau, 2) == 1)
    error('Input array size mismatch')
    % haha useless error message to troubleshoot
end

tr = tau(2:end)./tau(1:end-1);
a = mean(tr(end/2:end));

n = ((a.^2 + a.^-2 - 2) ./ (a + a.^-1 - 2)).^3;

k = size(x,1);
dxdt = zeros(size(x));
% First and last 2 points contain 0. Roughly first quarter of MSD has poor
% compliance to x(n+1) = a x(n).
for i = 3:k-2
    dxdt(i,:) = ( x(i - 2,:) - n * x(i - 1,:) + n * x(i + 1,:) - x(i + 2,:) ) ./ ( tau(i) .* (n .* (a + a.^-1 - 2) - (a.^2 + a.^-2 - 2)));
end