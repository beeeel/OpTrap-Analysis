function [FT] = msd_fourier_transformator(msdObj, obsT, varargin)

%% TO DO:
% Check through everything once
% Allow for drawing on extant figure and clearing/holding
% Check wRange behaviour
% Replace references to accumulated in main loop

%% Parse inputs

nMSDs = size(msdObj.msd, 1);

p = inputParser;

p.addRequired('msdObj',@(x)isa(x,'msdanalyzer')&&isscalar(x))
p.addRequired('obsT',@(x)validateattributes(x,{'numeric'},{'scalar'}))

p.addParameter('wRange',{},@(x)validateattributes(x,{'cell'},{'ncols',nMSDs}))
p.addParameter('trunc','none',@(x)any(strcmp(x,{'none','minima'})))
p.addParameter('extrap','none',@(x)any(strcmp(x,{'none','linear'})))
p.addParameter('norm','none',@(x)any(strcmp(x,{'none','low','high'})))
p.addParameter('show_int',false,@(x)islogical(x))

p.parse(msdObj, obsT, varargin{:});

%% Preparatory
% hee hee
if nMSDs == 2
    tits = {'Radial', 'Tangential'};    
elseif nMSDs == 4
    tits = {'Left Radial', 'Left Tangential', 'Right Radial', 'Right Tangential'};
else
    error('huh, how many beads?');
end
cols = {'k', [0.8 0 0], 'b', [0 0.75 0], 'k', 'k', 'k', 'k', 'k'};
lS = {'-','-','-','-','-','-','-','-','-'};
mS = {'none','none','none','none','none','none','none','none','none'}; {'s','o'};
if obsT < 0
    legs = sprintf('%i min before drug',abs(round(obsT)));
else
    legs = sprintf('%i min after drug',round(obsT));
end

yLs = [2e-6 1e-1];
fSz = 14;


figure(20)
prep_figure(tits, fSz, yLs, length(dims), show_ints, norm_mode);
clear h


for dimI = 1:length(dims)
    dim = dims(dimI);
    for rep = 1:length(rIdxs)
        rIdx = rIdxs(rep);
        
        % Number of points to skip
        nSkip = 40;
        
        % Get the MSD
        track = accumulated{dIdx}{1,cIdx}(rIdx).msd.msd{dim};
        tau = track(2:end-nSkip,1);
        msd = track(2:end-nSkip,2);
        
        % Determine truncation point
        [idx, eta] = msd_truncator(tau, msd, trunc_mode, -30*(dimI-1));
        
        % Extrapolate if necessary
        [tau, msd, eta, idx] = msd_extrapolator(tau, msd, idx, eta, extrap_mode);
        
        switch wRange_mode
            case 'auto'
                % Get plateau index for auto wRanging
                [pIdx, ~] = msd_truncator(tau, msd, 'minima', 0);
                edx = min(pIdx + 50, length(tau));
                sdx = max(pIdx-100, 1);
                wR = {1./tau([pIdx, sdx]), 1./tau([edx, pIdx])};
            case 'r-theta'
                wR = wRanges{dim};
            case 'all'
                wR = wRanges{dim}{rep};
            case 'struct'
                if ~isstruct(wRanges)
                    tmp = whos('wRanges');
                    error('wRanges needs to be a struct, instead got: %s\n', tmp.class)
                end
                day = dayDirs{dIdx};
                if isfield(wRanges, ['d' day])
                    c = num2str(cIdx);
                    if isfield(wRanges.(['d' day]), ['c' c])
                        if size(wRanges.(['d' day]).(['c' c]), 1) >= rIdx
                            wR = wRanges.(['d' day]).(['c' c]){rIdx,dim};
                        else 
                            wR = {};
                        end
                    else
                        wR = {};
                    end
                else
                    wR = {};
                end
            otherwise
                wR = wRanges{dIdx, cIdx};
        end
        
        % Do the interpolation and rheoFT
        [omega, G1, G2] = msd_interp_FT(tau(1:idx), msd(1:idx), 0, eta, idx, 1e3);
        
        % Find and show intercepts
        subplot(2+show_ints, length(dims),length(dims)* (1 + show_ints) + dimI)
        oC = gstar_interceptor(omega, G1, G2, wR, show_ints, cols{rep});
        
        allOCs{dim, rIdx} = oC;
        
        % Plot MSD
        subplot(2+show_ints, length(dims),dimI)
        
        switch norm_mode
            case 'intercept'
                nF = oC(end);
                xlim auto
            otherwise
                nF = 1;
        end
        
        h(rep) = loglog(nF.*tau, msd, 'LineWidth', 2, ...
            'Color', cols{rep}, 'LineStyle', lS{rep}, 'Marker', mS{rep});
        h(length(rIdxs)+1) = plot(nF.*tau(idx), msd(idx), ...
            'rx', 'LineWidth', 3, 'MarkerSize', 12);
        
        % Plot inverse of intercept frequencies
        yl = ylim;
        tmp = plot(nF.*[1; 1]./oC, yl, '--',...
            'Color',0.7*[1 1 1 0.8], 'LineWidth', 2);
        h(length(rIdxs)+2) = tmp(1);
        ylim(yl);
        
        %     legend('MSD','1 ÷ Intercept frequency','Location','best')
        
        
        % Plot the FT with intercept frequencies
        subplot(2+show_ints, length(dims),length(dims)+dimI)
        loglog(omega, ...
            G1, 'LineWidth', 2, ...
            'Color', cols{rep}, 'Marker', mS{rep}, 'LineStyle', '-', 'LineWidth', 2);
        loglog(omega, ...
            G2, 'LineWidth', 2, ...
            'Color', cols{rep}, 'Marker', mS{rep}, 'LineStyle', '--', 'LineWidth', 2);
        
        % Show Intercept frequency without changing YLims
        yl = ylim;
        plot(oC.*[1; 1], ylim, '--','Color',0.7*[1 1 1 0.8], 'LineWidth', 2)
        ylim(yl);
        
        legend('Storage','Loss','Intercept frequency','Location','best')
        
    end
    legend(h, legs{:},'Final point used in FT','1 ÷ Intercept frequency','Location','best')
end

%% Function definitions

function prep_figure(tits, fSz, yLs, n_dim, show_ints, norm_mode)
% Prepare figure window
clf
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
    ylabel('G'', G''''')
    title([tits{plt} ' Moduli'])
    set(gca,'FontSize',fSz)
    
    if show_ints
        subplot(3,n_dim,plt+n_dim*(1+show_ints))
        hold on
        title('Ratio of storage to loss moduli')
        xlabel('Frequency (Hz)')
        ylabel('Ratio G''/G''''')
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
        nP = 15;
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

function [dydx, varargout] = msd_gradientor(tau, msd, varargin)
% Calculates MSD gradient. I've made this so I can replace the pointwise
% numerical with a fitting method

if nargin == 2
    method = 'piecewise';
else
    method = varargin{1};
end
if nargin == 4
    nP = varargin{2};
end

switch method
    case 'piecewise'
        % Calculate pointwise gradient in log space
        tRs = tau(2:end)./tau(1:end-1);
        mRs = msd(2:end,:)./msd(1:end-1,:);
        
        dydx = log(mRs) ./ log(tRs);
        
        % Rolling average
        if ~exist('nP','var')
            nP = 15;
        end
        kern = ones(nP,1);
        dydx = conv(dydx, kern, 'valid')./nP;
        tout = conv(tau(1:end-1), kern, 'valid') ./ nP;
    case 'fitting'
        if ~exist('nP','var')
            nP = 40;
        end
        dydx = zeros(size(msd,1)-nP, size(msd,2));
        tout = nan(length(tau)-nP, 1);
        tau = tau(msd > 0); % Careful, if MSD is a matrix this will shit itself
        msd = msd(msd > 0);
        for jdx = 1:size(msd,2)
            for idx = 1:length(tau)-nP
                tdata = tau(idx:idx+nP);
                mdata = msd(idx:idx+nP,jdx);
                
                tout(idx) = mean(tdata);
                
                fo = fit(log(tdata), log(mdata), 'Poly1');
                dydx(idx, jdx) = fo.p1;
            end
        end
    otherwise
        error('Choose a correct method')
end

if nargout == 2
    varargout = {tout};
elseif nargout > 2
    error('Too many nargouts in msd_gradientor')
end

end

function oC = gstar_interceptor(omega, G1, G2, wRange, doPlot, colour)
% Find the intercepts between G1 and G2 as functions of frequency omega.
% Uses linear fit over ranges given in cell array wRange.

oC = zeros(1,length(wRange));

for wRIdx = 1:length(wRange)
    wdx = [0 0];
    wdx(2) = max(sum(omega > wRange{wRIdx}(1)), 1);
    wdx(1) = max(sum(omega > wRange{wRIdx}(2)), 1);
    
    odata = omega(wdx(1):wdx(2));
    Gdata = G1(wdx(1):wdx(2))./G2(wdx(1):wdx(2));
    
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
    h(2) = loglog(omega, G1./G2, 'Color', colour, 'LineWidth', 2);
    plot(xlim, [1 1], '--','Color',[1 1 1]*0.8, 'LineWidth', 3)
    try
        legend(h, 'Fit result', 'Ratio G''÷G''''','Location','best')
    catch ME
        error(ME.identifier, 'No wRanges supplied?')
    end
end
end
