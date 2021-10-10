% First look at processing event camera data
%

%% Load data
hn = system('hostname');
if strcmp(hn, 'bowtie')
    evtCamDir = '~/Documents/data/EventCam/';
    addpath('~/Documents/Analysis/prophesee-matlab-scripts');

else
    dirList = dir('~/mnt2');
    if length(dirList) == 2 ;system('sudo mount /dev/sda1 ~/mnt2'); end
    evtCamDir = '~/mnt2/Users/billh/Event_Camera/';
    addpath('~/junk/prophesee-matlab-scripts');

end
subDir = 'Event Camera - OpticallyTrappedBead/';
% subDir = 'Event Camera - Stage_SineWave/';
dirList = dir([evtCamDir subDir]);
ls([evtCamDir subDir])

fileBase = 'TrappedBead_1'; % 1, 2, or 3
% fileBase = 'Bead_Stage_SineWave'; % Bead or Fluorescence or Fluorescence..._smallsteps
% fileBase = 'Fluorescence_Stage_SineWave_smallsteps'; 

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
    
    try
        xyPos = fileread([evtCamDir subDir 'sCMOS_' fileBase '/XYPositions.txt']);
        xyPos = strsplit(xyPos, {'\t'});
        allPos = zeros(size(xyPos));
        for idx = 1:numel(xyPos)
            allPos(idx) = str2double(xyPos{idx});
        end
        allPos = reshape(allPos(1:end-1),floor(numel(allPos)/5),5);
    catch
        warning('woops no XYPositions.txt')
    end
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
                imStack{1}{fCount,1} = imread([evtCamDir subDir 'sCMOS_' fileBase '/' dirList(fIdx).name]);
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
            end
        end
        fIdx = fIdx + 1;
    end
end

ims = cat(3, imStack{1}{:,1});
imTs = (1:size(ims,3)) .*(metadata.FrameKey_49999_0_0.ElapsedTime_ms - metadata.FrameKey_0_0_0.ElapsedTime_ms) ... 
    ./ 5e4; % I sure hope there are this many frames in the imstack!
%%
vidObj = VideoWriter([evtCamDir subDir fileBase '_sCMOS.avi'],'Grayscale AVI')
open(vidObj)
fh = figure(1);
clf
colormap gray
for fr = 1:size(imStack{1},1)
    imagesc(imStack{1}{fr,1})
    axis off
    axis equal
    thisFrame = getframe(fh);
    thisFrame.cdata = thisFrame.cdata(:,:,1);
    writeVideo(vidObj, thisFrame);
end
vidObj
close(vidObj)

%%
centres = imCentreOfMass(ims(10:end-10, 10:end-10, :) .* uint16(ims(10:end-10, 10:end-10, :) > 500)); %2.2e4));
size(centres)

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
%%
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
histogram(cdEvents.x, 'BinLimits', [315 328]);
xlabel('X (px)')

subplot(4,2,6)
histogram(cdEvents.y, 'BinLimits', [243 257]);
xlabel('Y (px)')

subplot(4,2,8)
if isfield(cdEvents,'gray')
    histogram(cdEvents.gray);
    xlabel('Gray (arb. U)')
elseif isfield(cdEvents,'p')
    histogram(cdEvents.p);
    xlabel('Polarity (+/-)')
end

%% And interpolate
% How many events do you need per data point?
%% SCMOS data
thresh = 700;%2e4; % For trapped, 700? For sinewave, 2e4.
crop = [31 90 31 90]; %[400 1200 700 1400]; % For sinewave. For trapped, [31 90 31 90].
    
if exist('imStack','var')
    ims = cat(3, imStack{1}{:,1});
    ims = ims(crop(3):crop(4), crop(1):crop(2), :);
    %clear imStack

elseif ~exist('ims','var')
    error('Neither imStack nor ims exist as variables')
end
% % This is useful for determining the crop area
% ims = max(ims, [], 3);
% figure(1), clf, imagesc(ims)

% %%
% ims = ims(31:90, 31:90, :); % For trapped
% centres = imCentreOfMass(double(ims).*(ims > thresh), 'simple') + [30; 30];

centres = imCentreOfMass(double(ims).*(ims > thresh), 'simple') + crop(1:2)';
stiffIms = calcStiffness(centres, 0.216e-6) % μm/px is different for trapped bead (2x2 binning)
% stiffIms = 1e-7 * [8.37; 7.21]; % For TrappedBead_1. Ratio X/Y = 1.16.

% figure(5)
% clf
% subplot(2,1,1)
% scatter(centres(1,:), centres(2,:))
% axis image
% 
% subplot(2,2,3)
% histogram(centres(1,:))
% xlabel('X (px)')
% 
% subplot(2,2,4)
% histogram(centres(2,:))
% xlabel('Y (px)')

% Use y Range of centres to determine "height" (amplitude) of sine wave
amp = 0.108.*range(centres(2,:)); % For sineWave sets, 0.108μm/px.
% For sineWave, amp = 569.5284px=61.51μm, but from allPos, amp should be 600[μm??].

%% Buffer events over window (refinement ROI)
% Run the integration twice: Once with coarse time spacing and no ROI, and
% once with fine time spacing and an ROI based on the coarse spacing.


% nt_coarse = 1e4;
% [tsC, centresC] = bufferEvents(cdEvents.ts, cdEvents.x, cdEvents.y, nt_coarse);
% centresNoNaN = centresC(:,~isnan(centresC(1,:)));
% m = mean(centresNoNaN,2)';
% r = range(centresNoNaN,2)';
% ROI = [m-r/2, r]
% % lol that was a waste of time

% nt_fine = 1e3;
% ROI = [300 230 50 50];
% [tsF, centresF, ns] = bufferEvents(cdEvents.ts, cdEvents.x, cdEvents.y, nt_fine, ROI);
% tsNoNaN = tsF(:,~isnan(centresF(1,:)));
% centresNoNaN = centresF(:,~isnan(centresF(1,:)));
% % stiffEVFB = calcStiffness(centresNoNaN, 2.64e-6); % Dunno where I got this pixel calibration from... 
% stiffEVFB2 = calcStiffness(centresNoNaN, 0.164e-6); % Maybe should be 74e-6?

nt_fine = 1e3;
ROI = [300 230 50 50];
[tsF, centresF, ns] = bufferEvents_2(cdEvents.ts, cdEvents.x, cdEvents.y, nt_fine, ROI);
stiffBuffEv2 = calcStiffness(centresF, 2.64e-6); % stiffBuffEv2 = 1e-9 * [1.330; 1.127]; % For TrappedBead_1. Ratio X/Y = 1.1804

% plot(ROI([1 1 1 1 1])+[0 1 1 0 0].*ROI(3), ROI([2 2 2 2 2])+[0 0 1 1 0].*ROI(4), 'g-')

% With ROI width&height = 2*range
% stiffEVFB = 1e-11 * [0.831; 2.891]; % For TrappedBead_1. Ratio X/Y = 0.287.
% With ROI w&h = range
% stiffEVFB = 1e-11 * [3.368; 0.472]; % For TrappedBead_1. Ratio X/Y = 0.140.
% With ROI = [300 230 50 50]
% stiffEVFB = 1e-9 * [2.455; 2.101]; % For TrappedBead_1. Ratio X/Y = 1.168.

%% Interpolate to regularise
nT = 1e-3*range(tsF);

xin = centresF(1, :);
yin = centresF(2, :);

ti = linspace(min(tsF), max(tsF), nT);
xi = interp1(tsF, xin, ti, 'pchip');
yi = interp1(tsF, yin, ti, 'pchip');

msd = msdanalyzer(1, 'px', 'us', 'log');
tracks = {[ti' xi'], [ti' yi']};
msd = msd.addAll(tracks);
msd = msd.computeMSD;
msd.plotMSD
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
%% Animate 2.0
step = 0.5e3;
overlap = 250;
nSteps = 100;
figure(4)
clf
% hold on
subplot(10, 10, [1 10])
plot(cdEvents.ts(1:nSteps*step + 2*overlap)./1e6,'k', 'LineWidth', 2)
hold on
xlabel('#')
ylabel('t (s)')
subplot(10, 10, [21 100])
hold on
% ax = gca;
% scatter(200+allPos(:,2)./3, 225-allPos(:,3)./3)%, [], tGauss(end)*allPos(:,1)./allPos(1))
% axis equal
for idx = 1:nSteps
    im = zeros(480, 640);
    for jdx = (idx-1)*step + (1:(step+2*overlap))
        y = 1+cdEvents.y(jdx); 
        x = 1+cdEvents.x(jdx);
        im(y, x) = im(y,x) + cdEvents.p(jdx);
    end
    [t, c] = bufferEvents(...
        cdEvents.ts((idx-1)*step + (1:(step+2*overlap))), ...
        cdEvents.x((idx-1)*step + (1:(step+2*overlap))), ...
        cdEvents.y((idx-1)*step + (1:(step+2*overlap))), ...
        step + 2*overlap, ROI);
    t = t(step + 2*overlap);
    c = c(:,step + 2*overlap);
    subplot(10, 10, [1 10])
    cla
    i1 = (idx-1)*step + 1;
    i2 = (idx-1)*step + (step+2*overlap);
    plot(cdEvents.ts(1:nSteps*step + 2*overlap)./1e6, 1:nSteps*step + 2*overlap,'k', 'LineWidth', 2)
    plot(cdEvents.ts([i1 i2])./1e6, [i1 i2], 'y', 'LineWidth', 2)
    plot(t./1e6, (i1+i2)/2, 'rx','MarkerSize',6,'LineWidth',3)
    subplot(10, 10, [21 100])
    imagesc(im, [-1 1]*3);
    hold on
    plot(ROI([1 1 1 1 1])+[0 1 1 0 0].*ROI(3), ROI([2 2 2 2 2])+[0 0 1 1 0].*ROI(4), 'r-')
    plot(c(1), c(2), 'ks', 'LineWidth', 3, 'MarkerSize', 9)
    
    colorbar
    axis equal
    axis tight
    title(sprintf('%.1fs to %.1fs, %g events',1e-6*cdEvents.ts((idx-1)*step + 1),1e-6*cdEvents.ts((idx-1)*step + step+2*overlap), step+2*overlap))
%     fprintf('%.1fs to %.1fs, %g events\n',1e-6*cdEvents.ts((idx-1)*step + 1),1e-6*cdEvents.ts((idx-1)*step + step+2*overlap), step+2*overlap)
    pause(0.15)
end
%% Interpolate to regularise
nT = 1e-3*range(tsNoNaN);

[tIn, ia, ic]= unique(tsNoNaN);
xin = centresNoNaN(1, ia);
yin = centresNoNaN(2, ia);

ti = linspace(min(tsNoNaN), max(tsNoNaN), nT);
xi = interp1(tIn, xin, ti, 'pchip');
yi = interp1(tIn, yin, ti, 'pchip');

msd = msdanalyzer(1, 'px', 'us', 'log');
tracks = {[ti' xi'], [ti' yi']};
msd = msd.addAll(tracks);
msd = msd.computeMSD;
msd.plotMSD
%%
figure(73)
clf

for plt = 1:2
% subplot(2,1,plt)
loglog(msd.msd{plt}(:,1)./1e6, msd.msd{plt}(:,2))
hold on
end
xlabel('Delay time τ (s)')
ylabel('MSD (px.^2)')
title({'Mean-squared displacement' 'for trapped bead using event cam'})
grid on
%% Integrate events over window (manual ROI)
% Use an a priori ROI. You could use a single integration time to see where
% to place this
nt = 1e4;
ROI = [300 230 40 60];
[ts, centres] = integrateEvents(cdEvents.ts, cdEvents.x, cdEvents.y, nt, ROI);
figure(2)
clf
histogram(centres(1,:))
centresNoNaN = centres(:,~isnan(centres(1,:)));
stiffEVROI = calcStiffness(centresNoNaN, 2.64e-6);
% stiffEVROI = 1e-11 * [3.35; 2.36]; % For TrappedBead_1. Ratio X/Y = 1.422.

%% MSD for refinement method
msd = msdanalyzer(1, 'px', 's', 'log');
tracks = {[tsF', centresF(1,:)'], [tsF', centresF(2,:)']};
msd = msd.addAll(tracks);
msd = msd.computeMSD;

figure(8)
clf
hold on
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';

plot(msd.msd{1}(:,1), msd.msd{1}(:,2))
plot(msd.msd{1}(:,1), msd.msd{2}(:,2))