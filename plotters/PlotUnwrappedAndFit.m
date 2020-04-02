function PlotUnwrappedAndFit(unwrapped, Ia, info, FitEqn, FSs,frame)
%% PlotUnwrappedAndFit(unwrapped, Ia, info, FitEqn, FontSizes,frame)
% Plot the unwrapped data with fitted data overlayed
imagesc(1:360,(1:size(unwrapped,1))*0.07,unwrapped(:,:,frame))
hold on
plot(0.07 * Ia(:,:,frame),'r.')
plot(0.07 * FitEqn(info(frame).uMajorAxisLength, info(frame).uMinorAxisLength, ...
    info(frame).uOrientation, 1:360),'k:','LineWidth',3)
hold off
title({'Unwrapped cell with' 'fitting data and Centroid fitted data'},'FontSize',FSs.TFontSize)
xlabel('CW angle from +x (degrees)','FontSize',FSs.XFontSize)
ylabel('Radius (\mum)','FontSize',FSs.YFontSize)
legend('Column max','Fitted ellipse')
% add grid, change xticks
ax = gca;
ax.XTick = 0:90:360;
grid on
end