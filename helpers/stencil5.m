% Numerical derivative, need to find reference for 5 point stencil
function dxdt = stencil5(x, dt)
%% Compute the numerical derivative of each column of x
% dxdt = stencil5(x, dt)

% Lazy hack for array shapes
if size(dt,1) == size(x, 2) && size(dt, 2) == 1
    dt = dt.';
end

k = size(x,1);
dxdt = zeros(size(x));
% First and last 2 points contain 0.
for i = 3:k-2
    dxdt(i,:) = ( x(i - 2,:) - 8 * x(i - 1,:) + 8 * x(i + 1,:) - x(i + 2,:) ) ./ ( 12 * dt );
end