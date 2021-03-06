function PlotEllipseOverlay(Maj, Min, Orient, Centre,varargin)
%% PlotEllipseOverlay(MajorAxisLength, MinorAxisLength, Orientation, Centre)
% Take ellipse parameters (lengths in same unit as plot axis, orientation
% in radians) and centre co-ordinates (x, y) and plot them on the current
% axis.
if ~isempty(varargin)
    LSpec = varargin;
else
    LSpec = {'k--','LineWidth', 2};
end

theta = 0:0.01:2*pi;
X = 0.5 * Maj .* cos(theta) .* cos(Orient) ...
    - 0.5 * Min .* sin(theta) .* sin(-Orient)...
    + Centre(1);
Y = 0.5 * Maj .* cos(theta) .* sin(-Orient) ...
    + 0.5 * Min .* sin(theta) .* cos(Orient)...
    + Centre(2);
plot(X, Y, LSpec{:})
plot(Centre(1),Centre(2),'kx')