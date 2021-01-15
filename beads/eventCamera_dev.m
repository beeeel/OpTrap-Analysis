%% Have a look at the event camera data
% You need to set the path to the DAT file
filePath = '~/Data/phd/Bead_Stage_SineWave_td.dat';

% This data loading function comes from the Prophesee MATLAB scripts and
% needs to be on your MATLAB path
td_data = load_atis_data(filePath)

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