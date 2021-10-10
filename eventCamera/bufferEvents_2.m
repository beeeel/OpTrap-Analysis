function [ts, centres, ns] = bufferEvents_2(tIn, xIn, yIn, n_t, ROI)
%% [ts, centres, ns] = bufferEvents_2(t, x, y, n_t, [ROI])
%% Buffer events into sliding window with n_t events in at all times
% Optional: Discard events outside of ROI defined as [x y w h] (or [] for
% default). 

if ~exist('ROI', 'var') 
    ROI = [0 0 320 240];
elseif isempty(ROI)
    ROI = [0 0 320 240];
end

% Half-width and -height
w2 = floor(ROI(3)/2);
h2 = floor(ROI(4)/2);
% ROI centres
xc = ROI(1) + w2;
yc = ROI(2) + h2;
% Vectorize the in ROI calculation
inROI = abs(xIn - xc) < w2 & abs(yIn - yc) < h2;
tIdxs = inROI .* (1:length(tIn))';
tIdxs = tIdxs(tIdxs ~= 0);
xIn = xIn .* inROI;
yIn = yIn .* inROI;
tIn = tIn .* inROI;

% Prepare arrays for integration window and centres
xsum = 0;
ysum = 0;
nsum = 0;
tsum = 0;
centres(2, length(tIdxs) + n_t - 1) = 0;% Faster preallocation by implicit sizing
ts(1, length(tIdxs) + n_t - 1) = 0;
ns(1, length(tIdxs) + n_t - 1) = 0;

% Maybe not the most efficient way to reconsruct as it is not fully
% vectorised - could it be done with a rolling window now? Also, now it's
% harder to integrate an adaptive ROI into this, unless this function is
% called with only part of the event stream.
tic

% For each event
% (First do enough that the buffer is full)
for i = 1:n_t
    tIdx = tIdxs(i);
    xsum = xsum + xIn(tIdx);
    ysum = ysum + yIn(tIdx);
    nsum = nsum + inROI(tIdx);
    tsum = tsum + tIn(tIdx);
    
    ts(i) = tsum./nsum;
    ns(i) = nsum;
    centres(:,i) = [xsum; ysum] ./ nsum;
end
% (Then remove an event before adding another)
for i = n_t + 1:length(tIdxs)
    tIdx = tIdxs(i);
    bIdx = tIdxs(i - n_t);
    xsum = xsum - xIn(bIdx);
    ysum = ysum - yIn(bIdx);
    nsum = nsum - inROI(bIdx);
    tsum = tsum - tIn(bIdx);
    
    xsum = xsum + xIn(tIdx);
    ysum = ysum + yIn(tIdx);
    nsum = nsum + inROI(tIdx);
    tsum = tsum + tIn(tIdx);
    
    ts(i) = tsum./nsum;
    ns(i) = nsum;
    centres(:,i) = [xsum; ysum] ./ nsum;
end
% (Then just empty the buffer)
for i = length(tIdxs) + 1: length(tIdxs) + n_t - 1
    bIdx = tIdxs(i - n_t);
    xsum = xsum - xIn(bIdx);
    ysum = ysum - yIn(bIdx);
    nsum = nsum - inROI(bIdx);
    tsum = tsum - tIn(bIdx);
    
    ts(i) = tsum./nsum;
    ns(i) = nsum;
    centres(:,i) = [xsum; ysum] ./ nsum;
end
% Replace infs with nans
infs = any( isinf( centres ), 1);
centres(:, infs) = nan(2,sum(infs));
% Celebrate!
% fprintf('Integrated %g events into window length %g in %gs\n', length(tIn), n_t, toc)
