%% This needs to do what is says: Keep a record of which events are in the buffer, 
% and the state of the position/image space and time. Then it will be easy
% to implement the Hough transform if wanted.
function [ts, centres] = bufferEvents(tIn, xIn, yIn, n_t, ROI, pIn, eqTime)
%% centres = integrateEvents(t, x, y, n_t, [ROI, p])
%% Integrate events into sliding window of width n_t OR integrate events into n_t windows of equal time
% Optional: Discard events outside of ROI defined as [x y w h] (or [] for
% default). 
% Optional: Use polarity (you need to provide ROI also)
% Optional: Actually win

if ~exist('ROI', 'var') 
    ROI = [0 0 320 240];
elseif isempty(ROI)
    ROI = [0 0 320 240];
elseif length(ROI) == length(xIn)
    error('I think you gave me polarity but not ROI')
end
if ~exist('pIn', 'var') 
    pFun = '1';
elseif any(size(pIn) ~= size(xIn))
    pFun = '1';
else
    pFun = 'pIn(pIdx)';
end
if ~exist('eqTime', 'var') 
    eqTime = false;
elseif isempty(eqTime)
    eqTime = false;
end

% Half-width and -height
w2 = floor(ROI(3)/2);
h2 = floor(ROI(4)/2);
% ROI centres
xc = ROI(1) + w2;
yc = ROI(2) + h2;

% Width of each integration time
d_t = range(tIn)./n_t;     

if ~eqTime
    % This way: integrate with sliding window of width n_t
    
    % Preallocate maximum times array
    ts(length(tIn)) = 0;
    % 
    firstIdxFun = 'firstEventIdx + 1';
    lastIdxFun = 'lastEventIdx + 1';
    tsFun = 'something';
    % Function evaluated for each event processed
    perEvtFun = '';
else
    % Centre of each integration time
    ts = linspace(min(tIn)+d_t/2, max(tIn)-d_t/2, n_t);
    % End of each integration time
    tEs = ts + d_t/2;
    firstIdxFun = 'lastEventIdx + 1';
    % Function to calculate last index
    lastIdxFun = 'find(tIn > tEs(timeIdx),1) - 1';
    % Function evaluated for each event processed - 
    perEvtFun = '';
end

% Prepare arrays for integration window and centres
integrated(480, 640) = 0; % Faster preallocation by implicit sizing
centres(2,n_t) = 0;

% 0th time window ends at 0th event
lastEventIdx = 0 + ~eqTime * n_t;

% Probably not the most efficient way to reconsruct as it is not vectorised 
tic

% For each time window
for timeIdx = 1:length(ts)
    % Start the next window after the previous one
    firstEventIdx = eval(firstIdxFun);
    % Find the end of the window
    lastEventIdx = eval(lastIdxFun);
    % For each event in the window
    for pIdx = firstEventIdx:lastEventIdx
        % If it's within the ROI
        if abs(xIn(pIdx) - xc) < w2 && abs(yIn(pIdx) - yc) < h2
            % Get the array indices for that pixel
            y = yIn(pIdx) + 1;
            x = xIn(pIdx) + 1;
            % Add polarity to the integrated image. If user did not supply it,
            % it will have defaulted to 1 instead
            integrated(y, x) = integrated(y, x) + eval(pFun);
            % Do a per-event time update
            eval(perEvtFun);
        end
    end
    % Leave NaNs method
    centres(:,timeIdx) = imCentreOfMass(integrated, 'simple');

    %     % Remove NaNs method
    %     c = imCentreOfMass(integrated, 'simple');
    %     if any(isnan(c))
    %         c = [xc; yc];
    %     end
    %     centres(:,timeIdx) = c;
    
    % Reset integration array
    integrated = zeros(480, 640);
end
fprintf('Integrated %g events into %g windows in %gs\n', length(tIn), n_t, toc)