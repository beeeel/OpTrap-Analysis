function cTau = msd_cornerator(msdObj, obsT, tRanges, varargin)
%% TO DO: Write header comment
% Check that stuff actually works
% Incorporate into bead_processing_accumulator_v1
% Consider whether it's worth making bead_processor_v2 with detailed
%  analysis (should mostly be copypasta from accumulator_v1)
% Sleep

%% Parse inputs
nMSDs = size(msdObj.msd, 1);

p = inputParser;

p.addRequired('msdObj',@(x)isa(x,'msdanalyzer')&&isscalar(x))
p.addRequired('obsT',@(x)validateattributes(x,{'numeric'},{'scalar'}))
p.addRequired('tRanges',@(x)validateattributes(x,{'cell'},{'ncols',nMSDs}))

p.addParameter('nSkip', 40, @(x)validateattributes(x, {'numeric'},{'positive','<',length(msdObj.msd{1})}))
p.addParameter('dims', 1:nMSDs, @(x)validateattributes(x, {'numeric'},{'positive','nonzero','<=',nMSDs}))
p.addParameter('yLims', [1e-6 5e1], @(x)validateattributes(x, {'numeric'},{'increasing','positive','nonzero','numel',2}))
p.addParameter('figHand', [], @(x)isa(x,'matlab.ui.Figure'))
p.addParameter('lineColour', 'k', @(x)(isa(x,'char') && isscalar(x)) || (isa(x,'numeric') && all(x < 1) && length(x) == 3))
p.addParameter('lineStyle', '-', @(x) any(strcmp(x,{'-',':','-.','--','none'})))
p.addParameter('marker','none', @(x) any(strcmp(x,{'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram', 'none'})))
p.addParameter('interpM', 'pchip', @(x) any(strcmp(x, {'linear', 'nearest', 'next', 'previous', 'spline', 'pchip', 'cubic', 'v5cubic', 'makima'})))
p.addParameter('interpF', 1e2, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}))

p.parse(msdObj, obsT, tRanges, varargin{:});

% Intercept
tRanges = p.Results.tRange;
interpM = p.Results.interpM;
interpF = p.Results.interpF;
% MSD options
nSkip = p.Results.nSkip;
dims = p.Results.dims;
% Plot options
yl = p.Results.yLims;
fh = p.Results.figHand;
colour = p.Results.lineColour;
lS = p.Results.lineStyle;
mS = p.Results.marker;
%% Setup

% Titles for X and Y directions
if length(dims) == 2
    tits = {'Radial','Tangential'};
elseif length(dims) == 4
    tits = {'left bead Radial','left bead Tangential', 'right bead Radial', 'right bead Tangential'};
else 
    error('Length of dimensions needs to be 2 or 4 because I''m lazy')
end
% Font size
fSz = 16;

if isempty(fh)
    fh = figure(203);
    clf
else
    figure(fh)
end

if obsT < 0
    legs = sprintf('%i min before drug',abs(round(obsT)));
else
    legs = sprintf('%i min after drug',round(obsT));
end

cTau = zeros(numel(tRanges{:}),length(dims));

% For each dimension, get the tau and msd, interpolate and then fit all the
% tRanges
for dIdx = dims
    plt = dims(dIdx);
    ax = subplot(nTracks/2,2,plt);
    hold on
    clear h
    try
        tau = msdObj.msd{plt}(2:end-nSkip,1);
        msd = msdObj.msd{plt}(2:end-nSkip,2);
        
        taui = logspace(log10(tau(1)), log10(tau(end)), length(tau).*interpF)';
        msdi = interp1(tau, msd, taui, interpM);
        
        h = plot(taui, msdi, ...
            'Color',colour, 'LineWidth', 2, 'LineStyle', lS, 'Marker', mS);
        
        % Fit to interpolated data
        fps = zeros(2, length(tRanges));
        for fIdx = 1:length(tRanges)
            % Get the indexes for the specified time range
            msdIdx = find(taui > tRanges{fIdx}(1), 1) ...
                : find(taui < tRanges{fIdx}(2), 1, 'last');
            tauData = taui(msdIdx);
            msdData = msdi(msdIdx);
            % Fit to log data and show the fits
            fo = fit(log(tauData), log(msdData), 'poly1');
            plot(tauData, exp(fo.p1 * log(tauData) + fo.p2), 'k:', 'LineWidth', 1.5)
            % Store the data for later
            fps(:,fIdx) = [fo.p1; fo.p2];
        end
        % Find corners by intercept of fits
        dfps = diff(fps,1,2);
        cTau(:, dIdx) = exp(-dfps(2,:)./dfps(1,:));
    catch ME
        warning(['caught error: ' ME.message])
    end
    
    try
        legend(h,legs, 'Location', 'northwest');
    catch ME
        warning(['Couldn''t set legend: ' ME.message])
    end
    
    % This is horrible
    ax.XScale = 'log';
    ax.YScale = 'log';
    
    ax.FontSize = fSz;
    xl = [1 2];
    [xl(1), xl(2)] = bounds(msdObj.msd{1}(:,1));
    
    xlim(xl)
    ylim(yl)
    
    xlabel('Delay time \tau (s)')
    ylabel('Mean-squared displacement (\mu m^2)')
    title(sprintf('%s direction',tits{plt}))
    grid on
end

end
