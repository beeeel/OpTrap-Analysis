function fh = bead_plotDCData(data, varargin)
%% figureHandle = plotDCData(data, [fftLims, figNum])
% Plot trace and FFT for DC

if isfield(data.raw, 'dcAvg') && (min(size(data.raw.dcAvg)==size(data.raw.timeVecMs)) ~= 0)
    
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
    
    cropT = data.opts.cropT;
    timeVec = data.raw.timeVecMs(cropT(1):cropT(2));
    DC = data.raw.dcAvg(cropT(1):cropT(2));
    
    subplot(2,1,1)
    plot(1e-3*timeVec, DC,'.')
    
    xlabel('Time (s)')
    ylabel('Intensity (Arb. U.)')
    title('Mean brightness of ROI')
    
    subplot(2,1,2)
    fft_scaled(1e-3*timeVec, DC, true, gca);
    ylabel('Amplitude (Arb. U.)')
    if ~isempty(fftLims)
        ylim(fftLims)
    end
else
    warning('Expected data.raw.dcAvg, instead failing calmly')
    fh = [];
end