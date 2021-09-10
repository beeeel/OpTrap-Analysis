%% Have a look at the event camera data
% You need to set the path to the DAT file
filePath = '~/Data/phd/Bead_Stage_SineWave_td.dat';

% This is where I have the Prophesee MATLAB scripts
addpath ~/junk/prophesee-matlab-scripts/

% This data loading function comes from the Prophesee MATLAB scripts and
% needs to be on your MATLAB path
td_data = load_cd_events(filePath)
% td_data = load_atis_data(filePath)

% Have a look the first 5 events:
% Time is an integer, number of us. X and Y are pixel numbers. Polarity is
% either +1 or -1 (brighter or darker)
disp('Time  X     Y      Polarity')
disp([td_data.ts(1:5) td_data.x(1:5) td_data.y(1:5) td_data.p(1:5)] )

% Get a line of best fit to find average event rate
tFit = fit(td_data.ts,[1:length(td_data.ts)]','poly1');

%% Have a quick look at the distributions
figure(1)
clf

% Event rate
subplot(3,1,1)
plot(td_data.ts,1:length(td_data.ts),'.') 
hold on
plot(tFit)
legend('Event timestamps',['Average event rate = ' num2str(tFit.p1, 3) ' event/us'],'Location','best')
xlabel('Event number')
ylabel('Event timestamp (\mu s)')
title('When events are happening')

% X distribution of events
subplot(3,1,2)
histogram(td_data.x)
title('X Distribution of events')
xlabel('Pixel number')
ylabel('Count') 

% Y distribution of events
subplot(3,1,3)      
histogram(td_data.y)
title('Y Distribution of events')
xlabel('Pixel number')
ylabel('Count')

%% Reconstruct video in a "windowed" manner
% Probably not the most efficient way to reconsruct as it is not vectorised 

% Sample the video at fixed times (in us)
timesArr = 0:1e6:max(td_data.ts);

% Prepare an array which will contain the change to the image in each time
% window
windowedImage = zeros(480, 640, length(timesArr));

% First time window starts at first event
firstEventIdx = 1;
% For each time window
for timeIdx = 1:length(timesArr)
    % Except for the first window, add to the total from the previous
    % window
    if timeIdx ~= 1
        windowedImage(:,:,timeIdx) = windowedImage(:,:,timeIdx-1);
    end
    % Find the end of the window
    lastEventIdx = find(td_data.ts > timesArr(timeIdx),1) - 1;
    % For each event in the window
    for pIdx = firstEventIdx:lastEventIdx
        % Get the array indices for that pixel
        y = td_data.y(pIdx) + 1;
        x = td_data.x(pIdx) + 1;
        % Add that polarity to the cumulative image
        windowedImage(y, x, timeIdx) = windowedImage(y, x, timeIdx) + td_data.p(pIdx);
    end
    % Start the next window after this one
    firstEventIdx = lastEventIdx + 1;
end

% Calculate the cumulative image
cumulativeImage = cumsum(windowedImage,3);

% Calculate colour range to include most pixels
P = 99; % Percentile of pixels to include in colour scale
Y = prctile(abs(windowedImage(:,:,end)),P,'all');
%% Display a few frames with the same colour scale
% Pick which 4 to show
tShow = 1 + linspace(10, 40, 4);

figure(2)
clf
for n = 1:length(tShow)
    subplot(2,2,n)
    imagesc(windowedImage(:,:,round(tShow(n))),[-1 1] .* Y)
    title(['Time = ' num2str(timesArr(round(tShow(n)))./1e6) ' s'])
end
%% Some stuff I've moved from eventCam_processing_v0.m
%%
%% Gaussian filtering
sigma = 500;

gaussSz = 2*ceil(2*sigma)+1;
gauss = 1:gaussSz;
gauss = exp(-(gauss-ceil(gaussSz/1)).^2./sigma.^2);

% pGauss = conv(cdEvents.p, gauss, 'valid')./sum(gauss);
idx = 1:numel(cdEvents.x);
tmp = cdEvents.x(idx);
xGauss = conv(tmp(cdEvents.p(idx) > 0), gauss, 'valid')./sum(gauss);
% xGauss = conv(tmp, gauss, 'valid')./sum(gauss);
tmp = cdEvents.y(idx);
yGauss = conv(tmp(cdEvents.p(idx) > 0), gauss, 'valid')./sum(gauss);
% yGauss = conv(tmp, gauss, 'valid')./sum(gauss);
tmp = cdEvents.ts(idx);
tGauss = conv(tmp(cdEvents.p(idx) > 0), gauss, 'valid')./sum(gauss);
% tGauss = conv(tmp, gauss, 'valid')./sum(gauss);
tGauss = tGauss - min(tGauss);
clear tmp
% xGauss = conv(cdEvents.x(1:4:end), gauss, 'valid')./sum(gauss);
% yGauss = conv(cdEvents.y(1:4:end), gauss, 'valid')./sum(gauss);
% tGauss = conv(cdEvents.ts(1:4:end), gauss, 'valid')./sum(gauss);
%
figure(3)
clf
hold on
scatter(xGauss, yGauss, [], tGauss./1e3)
% scatter(200+allPos(:,2)./3, 225-allPos(:,3)./3, [], allPos(:,1))
title('Bead location from event camera')
%% Animate
step = 1e3;
overlap = 5e3;
nSteps = 200;
figure(4)
clf
hold on
ax = gca;
% scatter(200+allPos(:,2)./3, 225-allPos(:,3)./3)%, [], tGauss(end)*allPos(:,1)./allPos(1))
for idx = 1:nSteps
    scatter(xGauss((idx-1)*step + (1:(step+2*overlap))), ...
        yGauss((idx-1)*step + (1:(step+2*overlap))), ...
        [], ...30*nGauss((idx-1)*step + (1:(step+2*overlap)))/max(nGauss) , ...
        tGauss((idx-1)*step + (1:(step+2*overlap))))
    title(sprintf('%.1fs to %.1fs',1e-6*tGauss((idx-1)*step + 1),1e-6*tGauss((idx-1)*step + step)))
    xlim([100 400])
    ylim([100 350])
    pause(0.15)
    if idx < nSteps
        delete(ax.Children(1))
    end
end
%%
figure(99) % I got 99 figures but...
clf
ax2 = subplot(2,1,2);
ax1 = subplot(2,1,1);
for idx = 1:10
%     title(ax1, sprintf('%i events', (step+2*overlap)))
    histogram(ax1, xGauss((idx-1)*step + (1:(step+2*overlap))))
    histogram(ax2, yGauss((idx-1)*step + (1:(step+2*overlap))))
    drawnow
    pause(0.25)
end
%% Integrate some
map = [0 10 0
    0 9 0
    0 8 0
    0 7 0 
    0 6 0
    0 5 0
    5 5 5
    5 0 0
    6 0 0
    7 0 0 
    8 0 0
    9 0 0
    10 0 0]/10;
idxs = 3896:72712;
figure(10)
plot(idxs, cdEvents.ts(idxs))
xlabel('Event number')
ylabel('Time stamp (μs)')
title('Time distribution of events')

im = zeros(max(cdEvents.y)+1,max(cdEvents.x)+1);
for idx = idxs
    x = cdEvents.x(idx)+1;
    y = cdEvents.y(idx)+1;
    im(y, x) = im(y, x) + cdEvents.p(idx);
end
% im = 1.2.^im;
% im = log(im);
figure(11)    
imagesc(im)
colorbar
caxis([-1 1] * 20)
xlim([160 240])
ylim([190 260])
title(sprintf('Sum of events from %.2gs to %.2gs',cdEvents.ts(idxs(1))*1e-6,cdEvents.ts(idxs(end))*1e-6))
%% Gaussian weighted Gaussian filtering
sigma = 500;

gaussSz = 2*ceil(2*sigma)+1;
gauss = 1:gaussSz;
gauss = exp(-(gauss-ceil(gaussSz/1)).^2./sigma.^2)/(sigma*sqrt(2*pi));

% pGauss = conv(cdEvents.p, gauss, 'valid')./sum(gauss);
idx = 1:numel(cdEvents.x)/1024;
tmp = cdEvents.x(idx);
% xGauss = conv(tmp(cdEvents.p(idx) > 0), gauss, 'valid')./sum(gauss);
xGauss = conv(tmp, gauss, 'valid')./sum(gauss);
tmp = cdEvents.y(idx);
% yGauss = conv(tmp(cdEvents.p(idx) > 0), gauss, 'valid')./sum(gauss);
yGauss = conv(tmp, gauss, 'valid')./sum(gauss);
tmp = cdEvents.ts(idx);
% tGauss = conv(tmp(cdEvents.p(idx) > 0), gauss, 'valid')./sum(gauss);
tGauss = conv(tmp, gauss, 'valid')./sum(gauss);



%% MSD???
% tCrop = cdEvents.ts(ceil(gaussSz/2):ceil(444024-gaussSz/2));
tracks = {[tGauss, xGauss], [tGauss, yGauss]};

msd = msdanalyzer(1, 'px', 'us','original');
msd = msd.addAll(tracks);
msd = msd.computeMSD;

figure(4)
clf
msd.plotMSD;
ax = gca;
ax.XAxis.Scale = 'log';
ax.YAxis.Scale = 'log';
%% Load mat file with MSDs
dirName = '~/Data/phd/EventCam/';
fName = 'TrappedBead_3_td_msd.mat'; % MSD for set 1 doesn't exist
% fName = 'Bead_Stage_SineWave_td_msd.mat';
% fName = 'Fluorescence_Stage_SineWave_smallsteps_td_msd.mat';
load([dirName fName],'msd')
whos('msd')
% xGauss = msd.tracks{1}(:,2);
% yGauss = msd.tracks{2}(:,2);
% tGauss = msd.tracks{1}(:,1);
%% Plot loaded MSDs
pow = 1;
mult = 3e4;
% Xrange = [-4.2 -2];
Xrange = [-4.2 -1];

figure(7)
clf
ax = gca;
hold on
msd.plotMSD
ax.XAxis.Scale = 'log';
ax.YAxis.Scale = 'log';

X = logspace(Xrange(1), Xrange(2));
Y = mult.*X.^pow;

plot(X,Y,'--','LineWidth',2)

legend('X', 'Y','α τ^1')
%% % ↓ This was when I was trying to load the raw files with the dat reader... 
% %% Sort
% eventsSorted = struct('ts',[]);
% [eventsSorted.ts, I] = sort(cdEvents.ts);
% eventsSorted.x = cdEvents.x(I);
% eventsSorted.y = cdEvents.y(I);
% eventsSorted.gray = cdEvents.gray(I);
% %% Plot sorted
% figure(1)
% subplot(4,2,2)
% plot(eventsSorted.ts,'.','MarkerSize', 1);
% ylabel('Time (μs)')
% 
% subplot(4,2,4)
% plot(eventsSorted.x,'.','MarkerSize', 1);
% ylabel('X (px)')
% 
% subplot(4,2,6)
% plot(eventsSorted.y,'.','MarkerSize', 1);
% ylabel('Y (px)')
% 
% subplot(4,2,8)
% plot(eventsSorted.gray,'.','MarkerSize', 1);
% ylabel('Gray (arb. U)')

%%
%% Regularize time data for reasons
nTimes = 1e5;

sigma = 5*range(cdEvents.ts)./nTimes;

gaussWdth = 2*ceil(2*sigma)+1;

nEvents = numel(cdEvents.ts);
dT = floor(range(cdEvents.ts)/nTimes);

times = round(0:dT:(max(cdEvents.ts)-dT-1))+min(cdEvents.ts);
times2 = round(0:dT:(max(cdEvents.ts)-dT-1))+min(cdEvents.ts)+gaussWdth;

tIdxs = zeros(length(times),1);
tIdxs2 = zeros(length(times),1);

residuals = zeros(length(times),1);
residuals2 = zeros(length(times),1);
%% Get indexes for raw data
%     % create matrix of abs differences, find minimum  to choose as index
%     % take those indexes out of xGauss yGauss (tGauss??) and MSD with log
%     % spacing

idx = 0;
idxs = 1:2e4;

while idx < length(times)
    idx = idx + 1;

    % You can go your own waaaayyy
    [residuals(idx), tIdxs(idx)] = min(abs(cdEvents.ts(idxs) - times(idx)));

    % Needs better check for idxs going beyond max
    while idxs(end) < nEvents  && residuals(idx) == abs(cdEvents.ts(idxs(end)) - times(idx))
        idxs = idxs + min(15e3,nEvents-idxs(end));
        [residuals(idx), tIdxs(idx)] = min(abs(cdEvents.ts(idxs) - times(idx)));
    end
    tIdxs(idx) = tIdxs(idx) + idxs(1) - 1;
end
disp('timed')
idx = 0;
idxs = 1:2e4;
while idx < length(times2)
    idx = idx + 1;
    
    % You can go your own waaaayyy
    [residuals2(idx), tIdxs2(idx)] = min(abs(cdEvents.ts(idxs) - times2(idx)));

    while idxs(end) < nEvents && residuals2(idx) == abs(cdEvents.ts(idxs(end)) - times2(idx))
        idxs = idxs + min(15e3,nEvents-idxs(end));
        [residuals2(idx), tIdxs2(idx)] = min(abs(cdEvents.ts(idxs) - times2(idx)));
    end
    tIdxs2(idx) = tIdxs2(idx) + idxs(1) - 1;
end
disp('time2''d')
%% Group events and apply Gaussian
eventsGrouped = cell(length(tIdxs),1);
for idx = 1:length(tIdxs)
    idxs = tIdxs(idx):tIdxs2(idx);
%     eventsGrouped{idx} = [cdEvents.ts(idxs) cdEvents.x(idxs) cdEvents.y(idxs) cdEvents.p(idxs)];
    eventsGrouped{idx} = [cdEvents.ts(idxs) cdEvents.x(idxs) cdEvents.y(idxs)];
end
disp('grouped')
% Apply filter in parallel
xGauss = zeros(length(tIdxs),1);
yGauss = zeros(length(tIdxs),1);
tGauss = zeros(length(tIdxs),1);
nGauss = tIdxs2 - tIdxs;

parfor idx = 1:length(tIdxs)
    data = eventsGrouped{idx};
    
    muT = (times2(idx) + times(idx))*0.5;
    gaussT = (1/sqrt(2*pi*sigma.^2)) * exp(- (data(:,1) - muT).^2/sigma.^2);
    normFactor = sum(gaussT);
    
    tGauss(idx) = sum(data(:,1).*gaussT)./normFactor;
    xGauss(idx) = sum(data(:,2).*gaussT)./normFactor;
    yGauss(idx) = sum(data(:,3).*gaussT)./normFactor;
end
clear eventsGrouped
disp('done!')
%%
xSF = (337 - 141)/range(allPos(:,2));
ySF = (324 - 115)/range(allPos(:,3));
xOS = 210;
yOS = 225;
%
figure(6)
clf
hold on
scatter(xGauss, yGauss, 30*nGauss/max(nGauss), 0.5*(times+times2))
% scatter(xGauss, yGauss, [], tGauss)
scatter(xOS+allPos(:,2)*xSF, yOS-allPos(:,3)*ySF,[],'k','.')
legend('Event data','Stage position')
xlabel('X (px)')
ylabel('Y (px)')
title(sprintf('Regularized reconstructed data, gaussian σ = %0.1fms',1e-3*sigma))
%% Event distribution in time
idx = 3e5;
figure(7)
histogram(1e-6*cdEvents.ts(1e5:idx),100)
xlabel('Event time (s)')
ylabel('Count')
title(sprintf('Time distribution of first %.G events',idx))

