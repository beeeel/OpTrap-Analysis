%% Make a cool graphic showing how unwrap works
% Save file name. Don't put a file extension
fname = '~/png/test';

% Parameters for the ellipse
Maj = 20;
Min = 10;
Orient = 0;
Centre = [0,0];
LSpec = {'k','LineWidth',2};

ElEqn = @(a, b, phi, x) sqrt((a * cos(0.0175*x + phi)).^2 + (b * sin(0.0175*x + phi)).^2 );
Theta = 0:360;

fh = figure(9);
clf
subplot(1,2,1)
hold on
PlotEllipseOverlay(Maj, Min, Orient, Centre,LSpec{:})
axis image
xlim(0.55 * Maj * [-1 1])
ylim(0.55 * Maj * [-1 1])
title('Cartesian co-ordinates')
subplot(1,2,2)
hold on
% plot(ElEqn(Maj, Min, Orient, Theta(1)),LSpec{:})
ylim([0,Maj*1.1])
xlim([0,360])
title('Polar co-ordinates')

N_AddToGif(1,fh,fname)
for theta = 1:30:390
    subplot(1,2,1)
    N_DrawAngleLine([0,0],theta,1.05*Maj)
    subplot(1,2,2)
    N_DrawAngleLine([theta,0],90,1.05*Maj)
    plot(ElEqn(Maj, Min, Orient, Theta(1:theta)),LSpec{:})

    N_AddToGif(0,fh,fname)
end

function N_DrawAngleLine(start,angle,length)
% start is [x0, y0]. angle is angle above x-axis in deg.
xdata = start(1) * [1 1] + [0 length * cosd(angle)];
ydata = start(2) * [1 1] + [0 length * sind(angle)];
plot(xdata,ydata,'b--','LineWidth',1.5)
end

function N_AddToGif(frame,fh,fname)
drawnow
fr = getframe(fh);
im = frame2im(fr);
[imind,cm] = rgb2ind(im,256);
% Write to the GIF File
if frame == 1
    imwrite(imind,cm,strcat(fname,'.gif'),'gif', 'Loopcount',inf);
else
    imwrite(imind,cm,strcat(fname,'.gif'),'gif','WriteMode','append');
end
end