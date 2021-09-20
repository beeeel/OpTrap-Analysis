% It will be easy to implement the Hough transform if wanted.
function [ts, centres, ns] = bufferEvents(tIn, xIn, yIn, n_t, ROI)
%% centres = bufferEvents(t, x, y, n_t, [ROI])
%% Buffer events into sliding window of width n_t
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
xIn = xIn .* inROI;
yIn = yIn .* inROI;
tIn = tIn .* inROI;

% Prepare arrays for integration window and centres
xsum = 0;
ysum = 0;
nsum = 0;
tsum = 0;
centres(2, length(tIn) + n_t) = 0;% Faster preallocation by implicit sizing
ts(1, length(tIn) + n_t) = 0;
ns(1, length(tIn) + n_t) = 0;

% Maybe not the most efficient way to reconsruct as it is not fully
% vectorised - could it be done with a rolling window now? Also, now it's
% harder to integrate an adaptive ROI into this, unless this function is
% called with only part of the event stream.
tic

% For each event
% (First do enough that the buffer is full)
for tIdx = 1:n_t
    xsum = xsum + xIn(tIdx);
    ysum = ysum + yIn(tIdx);
    nsum = nsum + inROI(tIdx);
    tsum = tsum + tIn(tIdx);
    
    ts(tIdx) = tsum./nsum;
    ns(tIdx) = nsum;
    centres(:,tIdx) = [xsum; ysum] ./ nsum;
end
% (Then remove an event before adding another)
for tIdx = n_t + 1:length(tIn)
    bIdx = tIdx - n_t;
    xsum = xsum - xIn(bIdx);
    ysum = ysum - yIn(bIdx);
    nsum = nsum - inROI(bIdx);
    tsum = tsum - tIn(bIdx);
    
    xsum = xsum + xIn(tIdx);
    ysum = ysum + yIn(tIdx);
    nsum = nsum + inROI(tIdx);
    tsum = tsum + tIn(tIdx);
    
    ts(tIdx) = tsum./nsum;
    ns(tIdx) = nsum;
    centres(:,tIdx) = [xsum; ysum] ./ nsum;
end
% (Then just empty the buffer)
for tIdx = length(tIn) + 1: length(tIn) + n_t
    bIdx = tIdx - n_t;
    xsum = xsum - xIn(bIdx);
    ysum = ysum - yIn(bIdx);
    nsum = nsum - inROI(bIdx);
    tsum = tsum - tIn(bIdx);
    
    ts(tIdx) = tsum./nsum;
    ns(tIdx) = nsum;
    centres(:,tIdx) = [xsum; ysum] ./ nsum;
end
% Replace infs with nans
infs = any( isinf( centres ), 1);
centres(:, infs) = nan(2,sum(infs));
% Celebrate!
% fprintf('Integrated %g events into window length %g in %gs\n', length(tIn), n_t, toc)


%%  Something I didn't want to delete
% % Width of each integration time
% d_t = range(tIn)./n_t;  
%
% if ~eqTime
%     % This way: integrate with sliding window of width n_t
%     
%     % Preallocate maximum times array
%     ts(length(tIn)) = 0;
%     % 
%     firstIdxFun = 'firstEventIdx + 1';
%     lastIdxFun = 'lastEventIdx + 1';
%     tsFun = 'something';
%     % Function evaluated for each event processed
%     perEvtFun = '';
% else
%     % Centre of each integration time
%     ts = linspace(min(tIn)+d_t/2, max(tIn)-d_t/2, n_t);
%     % End of each integration time
%     tEs = ts + d_t/2;
%     firstIdxFun = 'lastEventIdx + 1';
%     % Function to calculate last index
%     lastIdxFun = 'find(tIn > tEs(timeIdx),1) - 1';
%     % Function evaluated for each event processed - 
%     perEvtFun = '';
% end
% 
% 
% 
% % 0th time window ends at 0th event
% lastEventIdx = 0 + ~eqTime * n_t;
%
% % For each time window
% for timeIdx = 1:length(ts)
%     % Start the next window after the previous one
%     firstEventIdx = eval(firstIdxFun);
%     % Find the end of the window
%     lastEventIdx = eval(lastIdxFun);
%     % For each event in the window
%     for pIdx = firstEventIdx:lastEventIdx
%         % If it's within the ROI
%         if inROI(pIdx)
%             % Get the array indices for that pixel
%             y = yIn(pIdx) + 1;
%             x = xIn(pIdx) + 1;
%             % Add polarity to the integrated image. If user did not supply it,
%             % it will have defaulted to 1 instead
%             integrated(y, x) = integrated(y, x) + eval(pFun);
%             % Do a per-event time update
%             eval(perEvtFun);
%         end
%     end
%     % Leave NaNs method
%     centres(:,timeIdx) = imCentreOfMass(integrated, 'simple'); 
%     % This could be done much better by tracking the centre of mass sums as we roll
% 
%     %     % Remove NaNs method
%     %     c = imCentreOfMass(integrated, 'simple');
%     %     if any(isnan(c))
%     %         c = [xc; yc];
%     %     end
%     %     centres(:,timeIdx) = c;
%     
%     % Reset integration array
%     integrated = zeros(480, 640);
% end