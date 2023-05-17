function area = piecewiseIntegration(t, x)
%% area = piecewiseIntegration(t, x)
% Perform piecewise integration of non-uniformly spaced function -
% basically draw a bar for each point you sample.

N = size(t,1);
if N == 1 && size(t,2) ~= 1
    N = size(t,2);
end
dt = diff(t);

area = dt(1) * x(1);

for idx = 2:N-1
    area = area + 0.5 * (dt(idx-1) + dt(idx)) * x(idx);
end

area = area + dt(end) * x(end);