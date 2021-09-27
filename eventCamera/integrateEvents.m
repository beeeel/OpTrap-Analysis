function [ts, centres] = integrateEvents(tIn, xIn, yIn, n_t, ROI, pIn)
%% centres = integrateEvents(t, x, y, n_t, [ROI, p])
%% Integrate events from cdEvents struct into n_t windows of equal time
% Optional: Discard events outside of ROI defined as [x y w h] (or [] for
% default). 
% Optional: Use polarity (you need to provide ROI also)

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

% Half-width and -height
w2 = floor(ROI(3)/2);
h2 = floor(ROI(4)/2);
% ROI centres
xc = ROI(1) + w2;
yc = ROI(2) + h2;

% Width of each integration time
d_t = range(tIn)./n_t;     

% Centre of each integration time
ts = linspace(min(tIn)+d_t/2, max(tIn)-d_t/2, n_t); 
% End of each integration time
tEs = ts + d_t/2;

% Prepare arrays for integration window and centres
integrated(480, 640) = 0; % Faster preallocation by implicit sizing
centres(2,n_t) = 0;

% 0th time window ends at 0th event
lastEventIdx = 0;

% Probably not the most efficient way to reconsruct as it is not vectorised 
tic

% For each time window
for timeIdx = 1:length(ts)
    % Start the next window after the previous one
    firstEventIdx = lastEventIdx + 1;
    % Find the end of the window
    lastEventIdx = find(tIn > tEs(timeIdx),1) - 1;
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