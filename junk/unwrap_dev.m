%% Measure cell by radial unwrapping
% Start with the centre of the cell, as found by find_cell, then draw lines
% radially to create line plots of intensity (interpolating as necessary).
% Find edge of cell from intensity plots, and then fit to parametric
% equation for ellipse.
%% Load
% Load imstack and info file - this contains the results from find_cell
CellType = 'LS174T';
Set = 'normoxia';
Num = '11';
[Imstack, info, meta] = LoadImstackInfoMeta(CellType, Set, Num);
%% Unpack stuff from info
% These were made by find_cell_v2
centres = [info.mCentres];
radii = min([info.mCentres] - min(size(Imstack{1}{1,1})));
%% Unwrap and fit to radial equation for ellipse
% Parameters
sc_up = 1.0;
n_theta = 360;
n_reps = 10;
tol = 0.15;

% Preallocate
r = (1:round(sc_up * max(radii)))';
fits = zeros(3,length(Imstack{1}));
theta = linspace(single(1),single(360),n_theta);
% From parametric equation of ellipse (x = a cos(t), y = b sin(t)), take to
% polar, r = (a.cos(f * theta + phi))^2 +  (b.sin(f * theta + phi))^2,
% where f is the frequency - twice per 2pi, theta is the angle (x below),
% and phi is the phase offset (orientation of ellipse)
fiteqn = @(a, b, c, x) (a.^2) *cos(0.0175*x + c).^2 + (b.^2)*sin(0.0175*x + c).^2 ;
%fiteqn = @(a, b, c, x) (max(a,b).^2) *cos(0.0175*x + c).^2 + (min(a,b).^2)*sin(0.0175*x + c).^2 ; % This makes worse fits

tic
% For memory efficiency, everything is in a single call to interp3, which
% is then cast to uint16. The first argument is the whole imagestack, the
% remaining are x, y, z query points, all arguments are cast to single for
% memory efficiency
unwrapped = uint16(interp3(...
    single(cat(3,Imstack{1}{:,1})),...
    single(reshape(centres(1,:),1,1,length(centres)) + r.*cos(theta*pi/180)),...
    single(reshape(centres(2,:),1,1,length(centres)) + r.*sin(theta*pi/180)),...
    single(repmat(reshape(1:length(Imstack{1}),1,1,length(Imstack{1})),length(r),n_theta))));
fprintf('Finished unwrapping in %gs\n',toc)
% Find the maxes and apply filtering based on deviation from median (will
% need large tolerance for large deformation)
[~, Ia] = max(unwrapped);
idxa = Ia > (1-tol)*median(Ia,2) & Ia < (1+tol)*median(Ia,2);
% For each frame, take angle corresponding to fit points and perform the
% fit
for frame = 1:length(Imstack{1})
    % Take corresponding angles, repeat them and unwrap
    th_fit = repmat(theta(idxa(:,:,frame)),1,n_reps) + ...
        reshape(360*(0:n_reps-1).*ones(length(theta(idxa(:,:,frame))),1),1,[]);
    fitobj = fit(double(th_fit'),repmat(Ia(:,idxa(:,:,frame),frame)',n_reps,1),fiteqn,...
        'Lower', [0 0 0], 'Upper', [inf inf pi], 'Start', [10, 10, 3]);
    fits(:,frame) = [fitobj.a, fitobj.b, fitobj.c];
    prog = ceil(100 * frame / length(Imstack{1}));
    fprintf('%s\r',['[' repmat('=',1,prog) repmat(' ',1,100-prog) ']'])
end
% Fix the equivalent fitting equations problem - if b>a
fits(3,:) = fits(3,:) + (fits(2,:) > fits(1,:)) * pi/2;
fits(3,:) = mod(fits(3,:)+pi/2, pi) - pi/2;
fits(1:2,:) = [max(fits(1:2,:)); min(fits(1:2,:))];
fprintf('%s\n',repmat(' ',1,104))
fprintf('Fitted data in %gs\n',toc)
%% Analysis of fits
figure(34)
subplot(4,1,1)
hold off
plot(0.07*fits(1,:).^2)
title('Fitted major and minor axis')
hold on
plot(0.07*fits(2,:).^2)
ylabel('axis (\mum)')
% ylim([110 120])

subplot(4,1,2)
hold off
plot(abs(fits(1,:) - fits(2,:))./(fits(1,:)+fits(2,:)))
title('Taylor deformation parameter')
hold on

subplot(4,1,3)
hold off
plot(0.07*centres')
title('Centres')
ylabel('position (\mum)')

subplot(4,1,4)
plot(fits(3,:),'.')
title('Orientation')
ylabel('\theta (rad)')
%% 1 frame plot
% Shows raw image and unwrapped image with the column maxes and fitted line
frame = 300;
xq = centres(1,frame) + r .* cos(theta*pi/180);
yq = centres(2,frame) + r .* sin(theta*pi/180);
figure(41)
subplot(2,1,1)
hold off
imagesc(Imstack{1}{frame,1})
hold on
plot(centres(1,frame), centres(2,frame),'k.','MarkerSize',16)
plot(centres(1,frame)+radii(frame)*cos(thetaR), ...
    centres(2,frame)+radii(frame)*sin(thetaR),':k','LineWidth',2)
plot(xq(:,1:30:end),yq(:,1:30:end),'k--')
axis image off
title('cell')
subplot(2,1,2)
title('Unwrapped cell')
hold off
a = imagesc(1:360,r*0.07,unwrapped(:,:,frame)); a.Parent.YAxis.Direction = 'normal';
hold on
plot(repmat((1:30:360),2,1),repmat(0.07*r([1,end]),1,12),'k--')
[~,I] = max(unwrapped(:,:,frame));
idx =  I > (1-tol)*median(I,2) & I < (1+tol)*median(I,2);
plot(theta(~idx), 0.07*I(~idx),'xr','MarkerSize',4)
plot(theta(idx), 0.07*I(idx),'.k')
plot(1:360, 0.07*fiteqn(fits(1,frame), fits(2,frame), fits(3,frame), 1:360),...
    'r-')
legend off
xlabel('Angle (degrees)'), ylabel('radius (\mum)'); 

%% Calculate the minimum deformation threshold to be detectable above noise
% First find at standard deviations for the first 30 frames
n_frames = 30;
n_sims = 1000;
n_theta = 360;
n_reps = 10;
Drange = [0 0.1];

theta = linspace(single(1),single(360),n_theta);
Dvals = linspace(Drange(1),Drange(2),n_sims);
sim_data = zeros(n_theta,n_sims);
sim_fits = zeros(3,n_sims);
R = sqrt(mean(Ia(idxa(:,:,1:n_frames)),'all'));
orientn = - pi/2 + pi*rand(1,n_sims);
rng('shuffle')
devs = zeros(n_frames,1);

for frame = 1:n_frames
    devs(frame) = std(Ia(:,idxa(:,:,frame),frame));
end
% Assume this deviation is the noise on a circle
noise = mean(devs);
%%
% Create some fictitious data using the radial equation and noise
% comparable to the noise found above, for a range of D values
tic
fprintf('\n')
for frame = 1:n_sims
    D = Dvals(frame);
    sim_data(:,frame) = fiteqn(R,R*(1-D)./(1+D),orientn(frame),...
        linspace(1,360,n_theta)) + normrnd(0,noise,1,n_theta);
    th_fit = repmat(theta,1,n_reps) + ...
        reshape((0:n_reps-1).*ones(n_theta,1),1,[]);
    fitobj = fit(double(th_fit'),repmat(sim_data(:,frame),n_reps,1),fiteqn,...
        'Lower', [0 0 -pi/2], 'Upper', [inf inf pi/2], 'Start', [10, 10, 3]);
    sim_fits(:,frame) = [fitobj.a, fitobj.b, fitobj.c];
    prog = ceil(100 * frame / n_sims);
    fprintf('%s\r',['[' repmat('=',1,prog) repmat(' ',1,100-prog) ']'])
end
% Fix the equivalent fitting equations problem - if b>a
sim_fits(3,:) = sim_fits(3,:) + (sim_fits(2,:) > sim_fits(1,:)) * pi/2;
sim_fits(3,:) = mod(sim_fits(3,:)+pi/2, pi) - pi/2;
sim_fits(1:2,:) = [max(sim_fits(1:2,:)); min(sim_fits(1:2,:))];
fprintf('%s\n',repmat(' ',1,104))
fprintf('Simulated data in %gs\n',toc)
%% Show results of simulations
sim_Dval = abs((sim_fits(1,:) - sim_fits(2,:)) ./ (sim_fits(1,:) + sim_fits(2,:)));
figure(113)
set(gcf,'WindowStyle','docked')
subplot(3,1,1)
hold off
plot(Dvals,sim_Dval,'.r')
hold on
plot(Dvals,Dvals,'k')
xlabel('D value'), ylabel('D value');
legend('fitted','ground truth','Location','best')
title('Fitted deformation compared to simulation truth')
subplot(3,1,2)
plot(Dvals,abs(sim_Dval./Dvals - 1),'.')
title('Relative error')
xlabel('D value'),ylabel('Relative error')
ylim([0 0.5])
subplot(3,1,3)
plot( sim_fits(3,:) - orientn,'.')
% subplot(3,2,5)
% plot(Dvals, sim_fits(3,:))
% title('Fitted orientation')
% subplot(3,2,6)
% plot(Dvals, orientn)
xlabel('D value'), ylabel('Orientation')
set(gca,'YTick',[-pi, -pi/2, 0, pi/2, pi],'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
title('Error on orientation')
%% 
figure(64)
frs = randi(n_sims, 1, 25);
for fr = 1:25
    subplot(5,5,fr)
%     frs = 780;
%     fr = 1;
    plot(sim_data(:,frs(fr))), hold on
    plot(fiteqn(sim_fits(1,frs(fr)), sim_fits(2,frs(fr)), sim_fits(3,frs(fr)),1:360))
    hold off
    title(num2str(frs(fr)))
end
