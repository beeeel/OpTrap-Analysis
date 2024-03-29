function fh = bead_plotDCData(data, varargin)
%% figureHandle = plotDCData(data, [fftLims, figNum])
% Plot trace and FFT for DC

%% Setup
% Only plot if DC is here for every point in timeVec
if ~(isfield(data.raw, 'dcAvg') && (min(size(data.raw.dcAvg)==size(data.raw.timeVecMs)) ~= 0))
    warning('Expected data.raw.dcAvg, instead failing calmly')
    fh = [];
    return
end

% Parse inputs
if nargin >= 2
    % help give better errors when misused!
    argN = 1;
    fftLims = varargin{argN};
    validateattributes(fftLims,{'numeric'},{'numel',2,'increasing'},...
        'bead_plotDCData','fftLims',nargin+argN-length(varargin));
else
    fftLims = [];
end

if nargin < 3
    fh = figure('Name',data.fName);
elseif nargin == 3
    fh = figure(varargin{2}); %#ok<*UNRCH>
else
    error('How many nargins did you use? Should be 1 to 3!')
end

timeVec = data.raw.timeVecMs;
cropT = data.opts.cropT;

timeVec = timeVec(cropT(1):cropT(2));
DC = data.raw.dcAvg(cropT(1):cropT(2));

%% Plotting
subplot(2,1,1)
plot(1e-3*timeVec, DC,'.')

xlabel('Time (s)')
ylabel('Intensity (Arb. U.)')
title('Mean brightness of ROI')

ax = subplot(2,1,2);
fft_scaled(1e-3*timeVec, DC, true, gca);
ylabel('Amplitude (Arb. U.)')
if ~isempty(fftLims)
    ylim(fftLims)
else
    xlim([10 ax.XLim(2)])
end

end
