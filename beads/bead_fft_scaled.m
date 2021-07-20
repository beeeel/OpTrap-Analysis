function data = bead_fft_scaled(data, doPlots, varargin)
%% data = bead_fft_scaled(data, doPlots, [setLims, fh])
% Calculate Fourier transform of centres data in both directions, store
% one-sided spectra in data struct along with frequency vector in Hz.
%
% Optional: Set axis y limits, provide figure handle

%% Setup
if nargin == 1
    doPlots = false;
end

% Don't ask.
argN = 1;
if nargin > 2 && ~ isempty(varargin{argN})
    % help give better errors when misused!
    setLims = varargin{argN};
    validateattributes(setLims,{'numeric'},{'numel',2},...
        'bead_fft_scaled','setLims',nargin+argN-length(varargin));
else
    setLims = [];
end
    

if length(data.opts.cropT) == 2
    cropT = data.opts.cropT;      
else
    cropT = [1 length(data.raw.timeVecMs)];
end

n_points = diff(cropT)+1;

% Choose centres row
if isfield(data.opts, 'centresRow')
    cRow = data.opts.centresRow;
elseif isfield(data.raw, 'suffixes')
    cRow = 1:length(data.raw.suffixes);
else
    cRow = 1:size(data.raw.xCentresPx,1);
end

%% Calculation
for direction = 'xy'
    % Get the FFT. Copied from MATLAB's example of scaling FFT to physical
    % units
    X = fft(data.mPerPx * data.raw.([ direction 'CentresPx'])(cRow,cropT(1):cropT(2)), [], 2);
    P = abs(X/n_points);
    P = P(:,1:floor(end/2)+1);
    P(:,2:end-1) = 2*P(:,2:end-1);
    
    data.pro.([direction 'fftM']) = P;
end

% Frequency in index (k+1) is k cycles per whole dataset, so convert to Hz
data.pro.fftFreqHz = (0:ceil((n_points-1)/2))./diff(data.raw.timeVecMs([cropT(1) cropT(2)])*1e-3);

%% Plot
if doPlots
    argN = 2;
    if nargin >= 4 && ~ isempty(varargin{argN})
        fh = varargin{argN};
    else
        fh = figure;
    end
    fh.Name = data.fName;
    
    if size(data.pro.xfftM,1) == 1
        subplot(2,1,1)
        plot(data.pro.fftFreqHz, 1e9 * data.pro.xfftM)
        if ~isempty(setLims)
            ylim(setLims)
        else
            xlim([10 fh.Children(1).XLim(2)])
        end
        
        title('X centres frequency spectrum')
        xlabel('Frequency (Hz)')
        ylabel('Amplitude (nm)')
        
        subplot(2,1,2)
        plot(data.pro.fftFreqHz, 1e9 * data.pro.yfftM)
        if ~isempty(setLims)
            ylim(setLims)
        end
        
        title('Y centres frequency spectrum')
        xlabel('Frequency (Hz)')
        ylabel('Amplitude (nm)')
    else
        n_obj = size(data.pro.xfftM,1);
        LR = data.raw.suffixes;
        for obj = 1:n_obj
            subplot(2,n_obj,obj)
            hold on
            plot(data.pro.fftFreqHz, 1e9 * data.pro.xfftM(obj,:))
            if ~isempty(setLims)
                ylim(setLims)
            else
                xlim([10 fh.Children(1).XLim(2)])
            end
            
            
            title(['X centres frequency spectrum for ' LR{obj} ' bead'])
            xlabel('Frequency (Hz)')
            ylabel('Amplitude (nm)')
            
            subplot(2,n_obj,n_obj + obj)
            hold on
            plot(data.pro.fftFreqHz, 1e9 * data.pro.yfftM(obj,:))
            if ~isempty(setLims)
                ylim(setLims)
            else
                xlim([10 fh.Children(1).XLim(2)])
            end
            
            title(['Y centres frequency spectrum for ' LR{obj} ' bead'])
            xlabel('Frequency (Hz)')
            ylabel('Amplitude (nm)')
        end
    end
    drawnow
end