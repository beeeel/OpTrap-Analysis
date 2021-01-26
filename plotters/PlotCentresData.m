function PlotCentresData(xCentres, yCentres, timeVec

fh.Name = fName;
clf

% Histogram of xCentres and yCentres in units um
subplot(3,1,1)
hold on
histogram(xCentresHPcrop.*1e9,'Normalization','probability')
histogram(yCentresHPcrop.*1e9,'Normalization','probability')
xlabel('Centre position (nm)')
ylabel('Bin probability')
title({['Histogram of ' num2str(fpass) 'Hz high-passed centres'], ['trap stiffness kx = ' num2str(xStiffHP./1e-6) ' pN/\mu m, ky = ' num2str(yStiffHP./1e-6) ' pN/\mu m']})
legend('X','Y')

% Scatterplot of each centre in units um
subplot(3,2,3)
plot(xCentresHPcrop.*1e6,yCentresHPcrop.*1e6,'.')
title({'Scatterplot of high-passed centres', [num2str(diff(cropTHP)+1) ' frames']})
xlabel('X Centre position (\mu m)')
ylabel('Y Centre position (\mu m)')
axis equal

% Histogram of residuals
subplot(3,2,4)
hold on
histogram(1e6.*(xCentresHPcrop - xCentresM(cropTHP(1):cropTHP(2))))
histogram(1e6.*(yCentresHPcrop - yCentresM(cropTHP(1):cropTHP(2))))
title({'Histogram of residuals', 'after high-pass filtering'})
xlabel('Correction amount (\mu m)')
ylabel('Count')
legend('X','Y')

% Time traces of X and Y in units um
timeVecHP = 1e-3.*timeVec(cropT(1):cropT(2));
timeVecHP = timeVecHP(cropTHP(1):cropTHP(2));
subplot(3,2,5)
plot(timeVecHP, xCentresHPcrop.*1e6,'.')
xlabel('Time (s)')
ylabel('X centre (\mu m)')
title({'High-passed X position' ' time trace'})
subplot(3,2,6)
plot(timeVecHP, yCentresHPcrop.*1e6,'.')
xlabel('Time (s)')
ylabel('Y centre (\mu m)')
title({'High-passed Y position' 'time trace'})
drawnow
end