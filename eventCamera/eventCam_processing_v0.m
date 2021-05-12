% First look at processing event camera data
%

%% Load data
dirList = dir('~/mnt2');
if length(dirList) == 2 ;system('sudo mount /dev/sda1 ~/mnt2'); end

addpath('~/junk/prophesee-matlab-scripts');

evtCamDir = '~/mnt2/Users/billh/Event_Camera/';
% subDir = 'Event Camera - OpticallyTrappedBead/';
subDir = 'Event Camera - Stage_SineWave/';

dirList = dir([evtCamDir subDir]);
ls([evtCamDir subDir])

% fileBase = 'TrappedBead_3'; % 1, 2, or 3
% fileBase = 'Bead_Stage_SineWave'; % Bead or Fluorescence or Fluorescence..._smallsteps
fileBase = 'Fluorescence_Stage_SineWave_smallsteps'; 

cdEvents = load_cd_events([evtCamDir subDir fileBase '_td.dat'])
%%
if strcmp(subDir(end-4:end),'Bead/')
    imStack = bfopen([evtCamDir subDir fileBase '_sCMOS/' fileBase '_MMStack_Pos0.ome.tif']);
    % This is at 0.216 um/px
    metadata = fileread([evtCamDir subDir fileBase '_sCMOS/' fileBase '_MMStack_Pos0_metadata.txt']);
    metadata = jsondecode(metadata);
    
    f = fields(metadata);
    ts = zeros(length(f)-1,1);
    for frIdx = 1:length(f)-1
        fr = num2str(frIdx-1);
        ts(frIdx) = metadata.(['FrameKey_' fr '_0_0']).ElapsedTime_ms;
    end
    ts = ts - ts(1);
    range(ts)./1e3
    metadata.Summary.Prefix
    
    xyPos = fileread([evtCamDir subDir 'sCMOS_' fileBase '/XYPositions.txt']);
    xyPos = strsplit(xyPos, {'\t'});
    allPos = zeros(size(xyPos));
    for idx = 1:numel(xyPos)
        allPos(idx) = str2double(xyPos{idx});
    end
    allPos = reshape(allPos(1:end-1),floor(numel(allPos)/5),5);
elseif strcmp(subDir(end-4:end),'Wave/')
    dirList = dir([evtCamDir subDir 'sCMOS_' fileBase '/']);
    imStack = {{}};
    nums = {};
    fIdx = 3; % Skip . and ..
    fCount = 1;
    while fIdx <= length(dirList)
        if length(dirList(fIdx).name) > 4
            if strcmp(dirList(fIdx).name(end-3:end), '.tif')
                % just assume the data will always be sorted correctly
%                 imStack{1}{fCount,1} = imread([evtCamDir subDir 'sCMOS_' fileBase '/' dirList(fIdx).name]);
                % This is at 0.108 um/px
                fCount = fCount + 1;
            elseif strcmp(dirList(fIdx).name(end-3:end),'.txt')
                %%
%                 fid = fopen([dirList(fIdx).folder '/' dirList(fIdx).name]);
                xyPos = fileread([dirList(fIdx).folder '/' dirList(fIdx).name]);
%                 fid = fclose(fid);
                xyPos = strsplit(xyPos, {'\t'});
                allPos = zeros(size(xyPos));
                for idx = 1:numel(xyPos)
                    allPos(idx) = str2double(xyPos{idx});
                end
                allPos = reshape(allPos(1:end-1),499,5);
                %{
                'XYPositions.txt' file: time and position data. The top row
                of this file is the time elapsed in ms at each stage
                position, the second and third rows are the detected X and
                Y stage position and the fourth and fifth rows are the set
                stage X and Y positions.
                [t; Xd; Yd; Xs; Ys]
                %}
                fCount = fCount + 1;
            end
        end
        fIdx = fIdx + 1;
    end
end

ims = cat(3, imStack{1}{:});
%%
centres = imCentreOfMass(ims .* uint16(ims > 2.2e4));
size(centres)
%% Load mat file with MSDs
dirName = '~/Data/phd/EventCam/';
% fName = 'TrappedBead_3_td_msd.mat';
% fName = 'Bead_Stage_SineWave_td_msd.mat';
fName = 'Fluorescence_Stage_SineWave_smallsteps_td_msd.mat';
load([dirName fName],'msd')
whos('msd')
xGauss = msd.tracks{1}(:,2);
yGauss = msd.tracks{2}(:,2);
tGauss = msd.tracks{1}(:,1);
%% Plot loaded MSDs
pow = 1;
mult = 3e4;
% Xrange = [-4.2 -2];
Xrange = [-3 -1];

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
%% Plot raw
figure(1)
clf 

subplot(4,2,1)
plot(cdEvents.ts,'.','MarkerSize', 1);
ylabel('Time (μs)')

subplot(4,2,3)
plot(cdEvents.x,'.','MarkerSize', 1);
ylabel('X (px)')

subplot(4,2,5)
plot(cdEvents.y,'.','MarkerSize', 1);
ylabel('Y (px)')

subplot(4,2,7)
if isfield(cdEvents,'gray')
    plot(cdEvents.gray,'.','MarkerSize', 1);
    ylabel('Gray (arb. U)')
elseif isfield(cdEvents,'p')
    plot(cdEvents.p,'.','MarkerSize', 1);
    ylabel('Polarity (+/-)')
end

subplot(2,2,2)
imagesc(imStack{1}{1,1})
axis image

subplot(2,2,4)
imagesc(imStack{1}{end,1})
axis image
%% Histograms!
subplot(4,2,2)
histogram(cdEvents.ts./60e6)
xlabel('Time (mins)')

subplot(4,2,4)
histogram(cdEvents.x);
xlabel('X (px)')

subplot(4,2,6)
histogram(cdEvents.y);
xlabel('Y (px)')

subplot(4,2,8)
if isfield(cdEvents,'gray')
    histogram(cdEvents.gray);
    xlabel('Gray (arb. U)')
elseif isfield(cdEvents,'p')
    histogram(cdEvents.p);
    xlabel('Polarity (+/-)')
end

% % ↓ This was when I was trying to load the raw files with the dat reader... 
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
%% Square windowing
kernSz = 1000;

kern = ones(kernSz,1);
pWindowed = conv(cdEvents.p, kern, 'valid')./sum(kern);
xWindowed = conv(cdEvents.x, kern, 'valid')./sum(kern);
yWindowed = conv(cdEvents.y, kern, 'valid')./sum(kern);
tWindowed = conv(cdEvents.ts, kern, 'valid')./sum(kern);

figure(2)
clf
scatter(xWindowed, yWindowed, [], tWindowed)
%% Gaussian filtering
sigma = 500;

gaussSz = 2*ceil(2*sigma)+1;
gauss = 1:gaussSz;
gauss = exp(-(gauss-ceil(gaussSz/1)).^2./sigma.^2);

% pGauss = conv(cdEvents.p, gauss, 'valid')./sum(gauss);
idx = 1:4:numel(cdEvents.x);
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
scatter(200+allPos(:,2)./3, 225-allPos(:,3)./3, [], allPos(:,1))
title('Bead location from event camera')
%%
step = 5e3;
overlap = 10e3;
nSteps = 200;
figure(4)
clf
hold on
ax = gca;
scatter(200+allPos(:,2)./3, 225-allPos(:,3)./3)%, [], tGauss(end)*allPos(:,1)./allPos(1))
for idx = 1:nSteps
    scatter(xGauss((idx-1)*step + (1:(step+2*overlap))), yGauss((idx-1)*step + (1:(step+2*overlap))), [] , tGauss((idx-1)*step + (1:(step+2*overlap))))
    title(sprintf('%.1fs to %.1fs',1e-6*tGauss((idx-1)*step + 1),1e-6*tGauss((idx-1)*step + step)))
    xlim([100 400])
    ylim([100 350])
    pause(0.15)
    if idx < nSteps
        delete(ax.Children(1))
    end
end
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
%% 
% save(msd,'eventcam_msd.mat')
%% SCMOS data
thresh = 700;

ims = cat(3, imStack{1}{:,1});
ims = ims(31:90, 31:90, :);
centres = imCentreOfMass(double(ims).*(ims > thresh), 'simple') + [30; 30];

figure(5)
clf
subplot(2,1,1)
scatter(centres(1,:), centres(2,:), [], ts)
axis image

subplot(2,2,3)
histogram(centres(1,:))
xlabel('X (px)')

subplot(2,2,4)
histogram(centres(2,:))
xlabel('Y (px)')
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
%%
eventsGrouped = cell(length(tIdxs),1);
for idx = 1:length(tIdxs)
    idxs = tIdxs(idx):tIdxs2(idx);
%     eventsGrouped{idx} = [cdEvents.ts(idxs) cdEvents.x(idxs) cdEvents.y(idxs) cdEvents.p(idxs)];
    eventsGrouped{idx} = [cdEvents.ts(idxs) cdEvents.x(idxs) cdEvents.y(idxs)];
end
disp('grouped')
%% Apply filter in parallel
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
% scatter(xGauss, yGauss, 30*nGauss/max(nGauss), 0.5*(times+times2))
scatter(xGauss, yGauss, [], tGauss)
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


