function [FT, varargout] = msd_fourier_transformator(msdObj, obsT, varargin)
%% [FT, [oCs]] = msd_fourier_transformator(msdObj, obsT, varargin)
% Do Fourier transform of MSD and find intercept frequency(s) if given.
% Takes parameters in name-value pairs. Possible options:
% wRange        - cell row vector with 1 element per dimension to be FT'd
% trunc         - truncation mode ('none', or 'minima')
% extrap        - extrapolation mode ('none', or 'linear')
% norm          - which corner to normalise time to ('none','low', or 'high')
% show_int      - show plots with tan(δ) used for intercept finding
% nSkip         - number of MSD points to skip (default 40)
% dims          - which dimensions to do FT on (indexes to msdObj.msd)
% yLims         - y limits on MSD plot
% figHand       - figure handle to plot upon (automatically holds existing plots if correct axes are present)
% lineColour    - line colour to plot
% lineStyle     - line style to plot MSD with (FT is always storage '-' and loss '--'
% marker        - line marker to plot MSD and FT with
%
% Could I make a "hold on" version of this to redraw on same axis?

%% Parse inputs
% If you're troubleshooting this bit, don't bother. Much easier to give up.
nMSDs = size(msdObj.msd, 1);
nargoutchk(1,2);

p = inputParser;

p.addRequired('msdObj',@(x)isa(x,'msdanalyzer')&&isscalar(x))
p.addRequired('obsT',@(x)validateattributes(x,{'numeric'},{'scalar'}))

p.addParameter('wRange',{{}, {}},@(x)validateattributes(x,{'cell'},{'ncols',nMSDs}))
p.addParameter('trunc','none',@(x)any(strcmp(x,{'none','minima'})))%,'FF'
p.addParameter('truncFF',0, @(x)validateattributes(x, {'numeric'},{'scalar','positive','nonzero','<',length(msdObj.msd{1})}))
p.addParameter('extrap','none',@(x)any(strcmp(x,{'none','linear'})))
p.addParameter('norm','none',@(x)any(strcmp(x,{'none','low','high'})))
p.addParameter('show_int',false,@(x)islogical(x))
p.addParameter('nSkip', 40, @(x)validateattributes(x, {'numeric'},{'positive','<',length(msdObj.msd{1})}))
p.addParameter('dims', 1:nMSDs, @(x)validateattributes(x, {'numeric'},{'positive','nonzero','<=',nMSDs}))
p.addParameter('yLims', [2e-6 1e-1], @(x)validateattributes(x, {'numeric'},{'increasing','positive','nonzero','numel',2}))
p.addParameter('fh', [], @(x)isa(x,'matlab.ui.Figure'))
p.addParameter('lineColour', 'k', @(x)(isa(x,'char') && isscalar(x)) || (isa(x,'numeric') && all(x < 1) && length(x) == 3))
p.addParameter('lineStyle', '-', @(x) any(strcmp(x,{'-',':','-.','--','none'})))
p.addParameter('marker','none', @(x) any(strcmp(x,{'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram', 'none'})))

p.parse(msdObj, obsT, varargin{:});

% FT options
wR = p.Results.wRange;
trunc_mode = p.Results.trunc;
FF = p.Results.truncFF;
extrap_mode = p.Results.extrap;
norm_mode = p.Results.norm;
show_ints = p.Results.show_int;
% MSD options
nSkip = p.Results.nSkip;
dims = p.Results.dims;
% Plot options
yLs = p.Results.yLims;
fh = p.Results.fh;
colour = p.Results.lineColour;
lS = p.Results.lineStyle;
mS = p.Results.marker;

%% Preparatory
% hee hee
if length(dims) == 2
    tits = {'Radial', 'Tangential'};    
elseif length(dims) == 4
    tits = {'Left Radial', 'Left Tangential', 'Right Radial', 'Right Tangential'};
else
    error('huh, how many beads?');
end

if obsT < 0
    legs = sprintf('%i min before drug',abs(round(obsT)));
else
    legs = sprintf('%i min after drug',round(obsT));
end

fSz = 16;

% silence warnings
warns = {'curvefit:fit:complexYusingOnlyReal','MATLAB:Axes:NegativeDataInLogAxis', 'MATLAB:legend:IgnoringExtraEntries'};
for wI = 1:length(warns)
    st(wI) = warning('query', warns{wI});
    warning('off', warns{wI});
end

prep_figure(fh, tits, fSz, yLs, length(dims), show_ints, norm_mode);
clear h

allOCs = {};
FT = cell(length(dims),1);

for dimI = 1:length(dims)
    dim = dims(dimI);
    
    % Get the MSD
    msdV = msdObj.msd{dim};
    tau = msdV(2:end-nSkip,1);
    msd = msdV(2:end-nSkip,2);
    
    
    % Determine truncation point
    [idx, eta] = msd_truncator(tau, msd, trunc_mode, FF);
    
    % Extrapolate if necessary
    [tau, msd, eta, idx] = msd_extrapolator(tau, msd, idx, eta, extrap_mode);
    
    % Do the interpolation and rheoFT
    [omega, G1, G2] = msd_interp_FT(tau(1:idx), msd(1:idx), 0, eta, idx, 1e3);
    
    FT{dimI} = [omega, G1, G2];
    
    % Find and show intercepts
    subplot(2+show_ints, length(dims),length(dims)* (1 + show_ints) + dimI)
    oC = gstar_interceptor(omega, G1, G2, wR{dimI}, show_ints, colour);
    
    allOCs{dim} = oC;
    
    % Plot MSD
    subplot(2+show_ints, length(dims),dimI)
    
    switch norm_mode
        case 'intercept'
            nF = oC(end);
            xlim auto
        otherwise
            nF = 1;
    end
    
    loglog(nF.*msdV(2:end,1), msdV(2:end,2), 'LineWidth', 2, ...
        'Color', [1 1 1 0.25], 'LineStyle', lS, 'Marker', mS);
    h = loglog(nF.*tau, msd, 'LineWidth', 2, ...
        'Color', colour, 'LineStyle', lS, 'Marker', mS);
    h(end+1) = plot(nF.*tau(idx), msd(idx), ...
        'rx', 'LineWidth', 3, 'MarkerSize', 12);
    
    % Plot inverse of intercept frequencies
    yl = ylim;
    tmp = plot(nF.*[1; 1]./oC, yl, '--',...
        'Color',0.7*[1 1 1 0.8], 'LineWidth', 2);
    h(end+1) = tmp(1);
    ylim(yl);
    
    %     legend('MSD','1 ÷ Intercept frequency','Location','best')
    
    
    % Plot the FT with intercept frequencies
    subplot(2+show_ints, length(dims),length(dims)+dimI)
    loglog(omega, ...
        G1, 'LineWidth', 2, ...
        'Color', colour, 'Marker', mS, 'LineStyle', '-', 'LineWidth', 2);
    loglog(omega, ...
        G2, 'LineWidth', 2, ...
        'Color', colour, 'Marker', mS, 'LineStyle', '--', 'LineWidth', 2);
    
    % Show Intercept frequency without changing YLims
    yl = ylim;
    plot(oC.*[1; 1], ylim, '--','Color',0.7*[1 1 1 0.8], 'LineWidth', 2)
    ylim(yl);
    
    legend('"Storage"','"Loss"','Intercept frequency','Location','best')
    
    legend(h, legs,'Final point used in FT','1 ÷ Intercept frequency','Location','best')
end

if nargout == 2
    varargout{1} = allOCs;
end

for wI = 1:length(warns)
    warning(st(wI).state, warns{wI});
end
% End of main function
end
%% Sub-function definitions

function prep_figure(fh, tits, fSz, yLs, n_dim, show_ints, norm_mode)
% Prepare figure window
if isempty(fh)
    figure(20)
    clf
else
    figure(fh.Number)
end
for plt = 1:n_dim
    subplot(2+show_ints, n_dim,plt)
    hold on
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    if strcmp( norm_mode, 'none' )
        xlabel('Lag time \tau (s)')
    else
        xlabel('Normalized lag time \tau')
    end
    ylabel('MSD (\mum^2)')
    title([tits{plt} ' MSD'])
    set(gca,'FontSize',fSz)
    xlim([1e-4 2e2])
    ylim(yLs)
    
    subplot(2+show_ints, n_dim,plt+n_dim)
    hold on
    set(gca,'XScale','log')
    set(gca,'YScale','log')
%     axis equal
%     ylim([1e-6 1e1])
    xlim([8e-3 1e4])
    xlabel('Frequency (Hz)')
    ylabel('"G'', G''''"')
    title([tits{plt} ' Fourier transform'])
    set(gca,'FontSize',fSz)
    
    if show_ints
        subplot(3,n_dim,plt+n_dim*(1+show_ints))
        hold on
%         title('Ratio of storage to loss moduli')
        xlabel('Frequency (Hz)')
        ylabel('"tan(\delta) = G''''÷G''"')
        set(gca,'FontSize',fSz)
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        grid on
    end
end

end

function [idx, eta] = msd_truncator(tau, msd, mode, FF)
% Tells you which idx to truncate MSD at (helpfully reuses gradient
% calculation)
dydxfilt = msd_gradientor(tau, msd);

% Take minima
[~, idx] = min(dydxfilt);

% Fudge factor extra points to include in FT
if ~exist('FF','var')
    FF = 40;
end

switch mode
    % Truncate at plateau/gradient minima + fudge factor
    case 'minima'
        idx = idx + FF;
    % Don't truncate
    case 'FF'
        idx = length(tau) - FF;
    otherwise
        idx = length(tau);
end
if idx > size(dydxfilt,1)
    eta = 1./ dydxfilt(end);
else
    eta = 1./ dydxfilt(idx);
end
end

function [tau, msde, eta, idx] = msd_extrapolator(tau, msd, idx, eta, mode)
% Needs to extrapolate MSD, replacing some points with (log) linear fit for
% the same tau values
switch mode
    case 'linear'
        nP = 30;
        fo = fit(log(tau(idx-nP:idx)), log(msd(idx-nP:idx)), 'Poly1');
        msde = msd;
        msde(idx-nP:end) = exp(fo.p1 * log(tau(idx-nP:end)) + fo.p2);
        
        eta = 1./fo.p1;
        idx = length(tau);
    otherwise
        msde = msd;
end

% figure(3)
% clf
% loglog(tau, msd,'-','LineWidth',2)
% hold on
% loglog(tau, msde,'-')
end
% 
% function [dydx, varargout] = msd_gradientor(tau, msd, varargin)
% % Calculates MSD gradient. I've made this so I can replace the pointwise
% % numerical with a fitting method
% 
% if nargin == 2
%     method = 'piecewise';
% else
%     method = varargin{1};
% end
% if nargin == 4
%     nP = varargin{2};
% end
% 
% switch method
%     case 'piecewise'
%         % Calculate pointwise gradient in log space
%         tRs = tau(2:end)./tau(1:end-1);
%         mRs = msd(2:end,:)./msd(1:end-1,:);
%         
%         dydx = log(mRs) ./ log(tRs);
%         
%         % Rolling average
%         if ~exist('nP','var')
%             nP = 15;
%         end
%         kern = ones(nP,1);
%         dydx = conv(dydx, kern, 'valid')./nP;
%         tout = conv(tau(1:end-1), kern, 'valid') ./ nP;
%     case 'fitting'
%         if ~exist('nP','var')
%             nP = 40;
%         end
%         dydx = zeros(size(msd,1)-nP, size(msd,2));
%         tout = nan(length(tau)-nP, 1);
%         tau = tau(msd > 0); % Careful, if MSD is a matrix this will shit itself
%         msd = msd(msd > 0);
%         for jdx = 1:size(msd,2)
%             for idx = 1:length(tau)-nP
%                 tdata = tau(idx:idx+nP);
%                 mdata = msd(idx:idx+nP,jdx);
%                 
%                 tout(idx) = mean(tdata);
%                 
%                 fo = fit(log(tdata), log(mdata), 'Poly1');
%                 dydx(idx, jdx) = fo.p1;
%             end
%         end
%     otherwise
%         error('Choose a correct method')
% end
% 
% if nargout == 2
%     varargout = {tout};
% elseif nargout > 2
%     error('Too many nargouts in msd_gradientor')
% end
% 
% end

function oC = gstar_interceptor(omega, G1, G2, wRange, doPlot, colour)
% Find the intercepts between G1 and G2 as functions of frequency omega.
% Uses linear fit over ranges given in cell array wRange.

oC = nan(1,max(1,length(wRange)));

if doPlot
    h(2) = loglog(omega, G1.\G2, 'Color', colour, 'LineWidth', 2);
    plot(xlim, [1 1], '--','Color',[1 1 1]*0.8, 'LineWidth', 3)
end

for wRIdx = 1:length(wRange)
    wdx = [0 0];
    wdx(2) = max(sum(omega > wRange{wRIdx}(1)), 1);
    wdx(1) = max(sum(omega > wRange{wRIdx}(2)), 1);
    
    odata = omega(wdx(1):wdx(2));
    Gdata = G1(wdx(1):wdx(2)).\G2(wdx(1):wdx(2));
    
    try
        fo = fit(log(odata), log(Gdata), 'Poly1');
        
        oC(wRIdx) = exp(- fo.p1 \ fo.p2);
        
        if oC(wRIdx) > odata(1) || oC(wRIdx) < odata(end) % Frequency omega is high to low
            oC(wRIdx) = nan;
        end
        
        if doPlot
            h = plot(odata, exp(fo.p1 * log(odata) + fo.p2), 'k:', 'LineWidth', 2.5);
        end
    catch ME
        warning(ME.identifier, 'Fitting failed with error %s\n',ME.message)
        oC(wRIdx) = nan;
    end
end
if doPlot
    try
        legend(h, 'Fit result', 'Ratio G''''÷G''','Location','best')
    catch
        try
            legend(h(2), 'Ratio G''''÷G''','Location','best')
        catch ME
            error(ME.identifier, 'No idea what happened in gstar_interceptor')
        end
    end
end
end
