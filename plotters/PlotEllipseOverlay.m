function PlotEllipseOverlay(Maj, Min, Orient, Centre)
%% PlotEllipseOverlay(MajorAxisLength, MinorAxisLength, Orientation, Centre)
% Take ellipse parameters (lengths in same unit as plot axis, orientation
% in radians) and centre co-ordinates (x, y) and plot them on the current
% axis.
theta = 0:0.01:2*pi;

plot(0.5 * Maj .* cos(theta) .* cos(Orient) ...
    - 0.5 * Min .* sin(theta) .* sin(-Orient)...
    + Centre(1),... % x values end here
    0.5 * Maj .* cos(theta) .* sin(-Orient) ...
    + 0.5 * Min .* sin(theta) .* cos(Orient)...
    + Centre(2),'k--','LineWidth',2)
plot(Centre(1),Centre(2),'kx')