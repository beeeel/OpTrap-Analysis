%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Compare different FT things           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First I did FT without any nSkip at the end - this gave poor results at
% low frequencies. Then I used nSkip = 40 (200 points MSDs), which gave
% good results when there was a plateau. Then I calculated the gradient at
% the end, and passed that to the FT as eta = 1/ginf, but it had little effect on the
% shape/noise.
%
%   DONE: Try finding plateau using gradient, truncate MSD at plateau and
%   extrapolate the plateau, giving the inverse gradient.
%
%   DONE: Better intercept fitting and ω range finding - use plateau as
%   midpoint to determine upper and lower ranges
%   TODO: Apply smoothing filter to G'/G'' to try and auto find crossings
%   TODO: Use initial quadratic fit, weighted by gaussian about plateau ω,
%   to find crossings and then improve with linear fits, no weighting.
%
%   TODO: Better gradient measurement through fitting or 5-point stencil or
%   similar.
%
%   TODO: Compare FTs with: truncate at plateau w/ginf = 0, all data w/ginf
%   from fit, replace data with fitted data
%
%   TODO: Use intercept times to normalize τ and MSD


%% Fourier transform shenanigans
% First make sure you've loaded accumulated 
if ~exist('accumulated','var')
    error('Need to have variable "accumulated"')
end

%% Main block - do FT and find crossover points

trunc_mode      = 'none';       % none or minima
extrap_mode     = 'none';       % none or linear
wRange_mode     = 'struct';      % none (manual), r-theta, all, struct or auto
norm_mode       = 'none';    % none or intercept
show_ints       = true;       % Show third row for intercepts

% Day set rep and dimension (x/y)
dIdx = 2;
cIdx = 1;
rIdxs = [1];
dims = 1:2;

% Frequency range over which to consider finding intercept, sorted by
% method

% for 'none'
% wRanges = {{[10 2e3]}, {[1 50]}}; 

% for 'all':
% wRanges = {{{[1e-2 0.9] [4e2 2e3]}, {[1e-2 0.2] [50 5e2]}},... % This for radial
%     { { [3 30] }, { [3 30] } }}; % This for tangential
% wRanges = {{{[1e-2 1.5] [10 300]}, {[1e-2 0.1] [10 3e2]}},... % This for radial
%     { { [1e-2 1.5] [8 80] }, { [3 30] } }}; % This for tangential

% for 'r-theta'
% wRanges = {{[0.01 0.5] [20 1000]}, {[0.01 0.5] [10 200]}, {[0.01 0.1] [50 1000]}, {[2 20]}};
% wRanges = {{[0.01 1] [3 100]}, {[0.01 2] [10 1000]}}; % for 'r-theta'

% For 'struct':
[wRanges, ~] = get_wRanges_struct(dayDirs);


nB = length(accumulated{dIdx}{1,cIdx}(rIdxs(1)).msd.msd);
% hee hee
if nB == 2
    tits = {'Radial', 'Tangential'};    
elseif nB == 4
    tits = {'Left Radial', 'Left Tangential', 'Right Radial', 'Right Tangential'};
else
    error('huh, how many beads?');
end
cols = {'k', [0.8 0 0], 'b', [0 0.75 0], 'k', 'k', 'k', 'k', 'k'};
lS = {'-','-','-','-','-','-','-','-','-'};
mS = {'none','none','none','none','none','none','none','none','none'}; {'s','o'};
legs = cell(size(rIdxs));
for idx = 1:length(rIdxs)
    t = accumulated{dIdx}{2,cIdx}(rIdxs(idx));
    if t < 0
        legs{idx} = sprintf('%i min before drug',abs(round(t)));
    else
        legs{idx} = sprintf('%i min after drug',round(t));
    end
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

% %% G' G'' gradients
% [a1, t1] = msd_gradientor(omega, G1);
% [a2, t2] = msd_gradientor(omega, G2);
% figure(1)
% % clf
% semilogx(t1, a1, 'bs')
% hold on
% semilogx(t2, a2, 'bx')
% 
% xlabel('Frequency f (Hz)')
% ylabel('Power law scaling of moduli')
% legend('Storage (piecewise)', 'Loss (piecewise)', ...
%     'Storage (fitting)', 'Loss (fitting)')
% %% MSD gradient
% % t = logspace(-2, 2, 100)';
% % maxMode = 11;
% % g = zeros(size(t));
% % for n = 1:maxMode
% %     g = g + (1/n.^4) * (1 - exp(-t * n.^4));
% % end
% % 
% % % This doesn't seem to work
% % gg = zeros(size(t));
% % for n = 1:maxMode
% %     gg = gg + (1/n.^4 + exp(-t * n.^4));
% % end
% t = tau;
% g = msd;
% 
% [a1, t1] = msd_gradientor(t, g, 'piecewise', 15);
% [a2, t2] = msd_gradientor(t, g, 'fitting', 40);
% [a3, t3] = msd_gradientor(t, g, 'fitting', 3);
% 
% % Plot it
% figure(2)
% clf
% % semilogx(t, gg, 'b^')
% yyaxis left
% semilogx(t1, a1, 'ks')
% hold on
% semilogx(t2, a2, 'kx')
% semilogx(t3, a3, 'ko')
% 
% ax = gca;
% ax.XScale = 'log';
% 
% title('Simulated MSD gradient measured with two methods')
% xlabel('Delay  time τ (s)')
% ylabel('MSD gradient')
% 
% yyaxis right
% ax = gca;
% ax.YScale = 'log';
% loglog(t, g, 'r.')
% ylabel('MSD (arb. units)')
% 
% legend(...'Analytical soluion', 
%     'Piecewise numerical method', 'Linear fit method, n=40', ...
%     'Linear fit method, n=5', 'Simulated MSD', 'Location', 'east')
% grid on


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

function [wRanges, rIdxs] = get_wRanges_struct(dayDirs)
rIdxs = 1;
wRanges = struct(strcat('d',dayDirs{end}), struct('c1', {{}}));

wRanges.d2021_07_27.c2 = {...
    {[2e-2 1] [5 1e2]} {[2e-2 1] [5 1e2]};
    {[2e-2 1] [5 1e2]} {[2e-2 1] [5 1e2]};
    {[2e-2 0.8] [8e2 2e3]} {[8 80]};
    {[2e-2 0.8] [8e2 2e3]} {[8 80]};
    {[2 10]} {[2 10]};
    {[2 10]} {[2 10]};
    {[1 7]} {[1 10]};
    {[1e-2 1.5e-1] [30 1e4]} {[2 10]};};
wRanges.d2021_07_27.c1 = {...
    {[1e-2 2] [200 2e3]} {[1e-2 5] [60 8e2]};
    {[1e-2 5] [500 5e4]} {[1e-2 6] [200 1e4]};
    {[1e-2 2] [8e2 2e3]} {[1e-2 1] [200 1e4]};
    {[1e-2 2] [8e2 2e3]} {[1e-1 2] [200 1e4]};
    {[1e-2 1] [8e2 2e3]} {[1e-2 1] [200 1e4]};
    {[1e-1 3] [20 1e2]} {[1e-2 1] [200 1e4]};
    {[1e-2 1] [6e2 2e3]} {[1e-2 1] [200 1e4]};
    {[1e-2 1] [6e2 2e3]} {[1e-2 1] [200 1e4]};
    {[1e-2 1] [6e2 2e3]} {[1e-2 1] [200 1e4]}};
wRanges.d2021_07_26.c1 = {...
    {[1e-2 10] [50 3e2]} {[1e-2 10] [20 8e2]};
    {[1e-2 5] [50 5e2]} {[1e-2 6] [20 8e2]};
    {[1e-2 10] [50 3e2]} {[1e-2 10] [20 8e2]};
    {[6e-2 2] [50 3e2]} {[6e-2 2] [20 8e2]};
    {[1e-2 1] [50 3e2]} {[1e-2 1] [20 8e2]};
    {[1e-1 3] [20 1e2]} {[1e-2 3] [20 2e2]};
    {[1e-1 3] [50 3e2]} {[1e-2 3] [50 3e2]};
    {[1e-1 3] [10 2e2]} {[1e-2 3] [10 2e2]}};
wRanges.d2021_07_23.c3 = {...
    {[100 1e4]} {[20 8e2]};
    {[100 1e4]} {[20 8e2]};
    {[100 1e4]} {[20 8e2]};
    {[100 1e4]} {[20 8e2]};
    {[100 1e4]} {[20 8e2]};
    {[100 1e4]} {[20 8e2]};
    {[1e-2 0.8] [100 1e4]} {[100 1e4]};
    {[100 1e4]} {[20 8e2]};
    {[100 1e4]} {[20 1e3]}};
wRanges.d2021_07_23.c2 = {...
    {[1e-2 1] [70 1e4]} {[5 80]};
    {[3 50]} {[2 30]};
    {[3 50]} {[0.5 5]};
    {[3 50]} {[0.5 5]};
    {[0.5 4]} {[0.5 4]};
    {[3e-1 5]} {[3e-1 3]};
    {[1 10]} {[1 10]}; %
    {[2e-1 10]} {[2e-1 8]};
    {[4e-1 3]} {[1e-1 3]};
    {[1 5]} {[1 5]};};
wRanges.d2021_07_23.c1 = {...
    {[1e-1 5]} {[1e-1 5]};
    {[5e-1 500]} {[1e-1 5]};
    {[1e-1 5]} {[1e-1 5]};
    {[1e-1 3]} {[1e-1 3]};
    {[1e-1 3]} {[1e-1 3]};
    {[1e-1 5]} {[1e-1 5]};
    {[1e-2 2]} {[1e-2 2]}; %
    {[1e-2 30]} {[1e-2 8]};
    {[1e-2 1] [20 200] [300 1e3]} {[1e-2 2]};
    {[1e-2 1] [20 200] [300 1e3]} {[1e-2 2]}};

end