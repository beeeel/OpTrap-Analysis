function [FT, varargout] = msd_fourier_transformator(msdObj, obsT, varargin)
%% [FT, [oCs]] = msd_fourier_transformator(msdObj, obsT, varargin)
% Do Fourier transform of MSD and find intercept frequency(s) if given.
% Takes parameters in name-value pairs. Possible options:
% wRange        - cell row vector with 1 element per dimension to be FT'd
% nBead         - numerical scalar for number of beads in dataset
% trunc         - truncation mode ('none', 'minima', 'FF', or 'timeFF')
% truncFF       - fudge factor for truncation (points skipped)
% truncT        - time to truncate at (when using trunc = 'timeFF')
% extrap        - extrapolation mode ('none', or 'linear')
% norm          - which corner to normalise time to ('none','low', or 'high')
% show_int      - show plots with tan(δ) used for intercept finding
% nSkip         - number of MSD points to skip (default 40)
% dims          - which dimensions to do FT on (indexes to msdObj.msd)
% yLims         - y limits on MSD plot
% fh            - figure handle to plot upon (automatically holds existing plots if correct axes are present)
% lineColour    - line colour to plot
% lineStyle     - line style to plot MSD with (FT is always storage '-' and loss '--'
% marker        - line marker to plot MSD and FT with
% lowPassFreq   - Apply low pass to FT (WIP)
% interpF       - Variable interpolation factor (default 1000)
% msdNorm       - Normalize MSD before FT.
% showLeg       - Show legend on plots
% doPlot        - What it says

%% Parse inputs
% If you're troubleshooting this bit, don't bother. Much easier to give up.
nMSDs = size(msdObj.msd, 1);
nargoutchk(1,2);

p = inputParser;

p.addRequired('msdObj',@(x)(isa(x,'msdanalyzer')&&isscalar(x)) || (false));
p.addRequired('obsT',@(x)validateattributes(x,{'numeric'},{'scalar'}))

p.addParameter('nBead',[],@(x)validateattributes(x,{'numeric'},{'scalar','positive'}))
p.addParameter('wRange',{{}, {}},@(x)validateattributes(x,{'cell'},{'ncols',nMSDs}))
p.addParameter('trunc','none',@(x)any(strcmp(x,{'none','minima','FF', 'timeFF'})))%
p.addParameter('truncFF',0, @(x)validateattributes(x, {'numeric'},{'scalar','positive','nonzero','<',length(msdObj.msd{1}), 'integer'}))
p.addParameter('truncT',[], @(x)validateattributes(x, {'numeric'},{'positive','nonzero','<',msdObj.msd{1}(end,1)}))
p.addParameter('extrap','none',@(x)any(strcmp(x,{'none','linear'})))
p.addParameter('eta',[], @(x)validateattributes(x, {'numeric'},{'scalar','positive','nonzero'}))
p.addParameter('norm','none',@(x)any(strcmp(x,{'none','low','high'})))
p.addParameter('show_int',false,@(x)islogical(x))
p.addParameter('nSkip', 40, @(x)validateattributes(x, {'numeric'},{'positive','<',length(msdObj.msd{1})}))
p.addParameter('dims', 1:nMSDs, @(x)validateattributes(x, {'numeric'},{'positive','nonzero','<=',nMSDs}))
p.addParameter('yLims', [], @(x)validateattributes(x, {'numeric'},{'increasing','positive','nonzero','numel',2}))
p.addParameter('fh', [], @(x)isa(x,'matlab.ui.Figure'))
p.addParameter('lineColour', 'k', @(x)(isa(x,'char') && isscalar(x)) || (isa(x,'numeric') && all(x <= 1) && length(x) == 3))
p.addParameter('lineStyle', '-', @(x) any(strcmp(x,{'-',':','-.','--','none'})))
p.addParameter('marker','none', @(x) any(strcmp(x,{'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram', 'none'})))
p.addParameter('lowPassFreq',[],@(x) validateattributes(x, {'numeric'},{'positive','scalar','nonzero'}))
p.addParameter('interpF', 1e3, @(x)validateattributes(x, {'numeric'},{'positive','scalar'}))
p.addParameter('msdNorm', ones(1,nMSDs), @(x)isa(x,'double') && length(x) == nMSDs && all(x > 0))
p.addParameter('showLeg', true, @(x) isa(x,'logical'))
p.addParameter('doPlot', true, @(x) isa(x,'logical'))

p.parse(msdObj, obsT, varargin{:});

nB = p.Results.nBead;
% FT options
wR = p.Results.wRange;
trunc_mode = p.Results.trunc;
FF = p.Results.truncFF;
truncT = p.Results.truncT;
extrap_mode = p.Results.extrap;
etaIn = p.Results.eta;
norm_mode = p.Results.norm;
show_ints = p.Results.show_int;
lpFrq = p.Results.lowPassFreq;
% MSD options
nSkip = p.Results.nSkip;
dims = p.Results.dims;
msdNorm = p.Results.msdNorm;
if numel(msdNorm) ~= numel(dims)
    try 
        msdNorm = msdNorm(dims);
    catch
        error('Normalisation and dimensions mismatch that you thought wouldn''t happen. Well guess what. It has happened');
    end
end
% Plot options
yLs = p.Results.yLims;
fh = p.Results.fh;
colour = p.Results.lineColour;
lS = p.Results.lineStyle;
mS = p.Results.marker;
interpF = p.Results.interpF;
showLeg = p.Results.showLeg;
doPlot = p.Results.doPlot;

%% Preparatory
% hee hee
if length(dims) == 2
    tits = {'Radial', 'Tangential'};    
elseif length(dims) == 4
    tits = {'Left Radial', 'Left Tangential', 'Right Radial', 'Right Tangential'};
elseif length(dims) == nB
    warning('New code: Use average of MSDs for FT')
    tits = {'Thing'};
    dims = 1;
elseif length(dims) == 1
    tits = {'Radial', 'Tangential'};    
    tits = tits(dims);
else
    error('huh, how many beads?');
end

if obsT < 0
    legs = sprintf('%i min before',abs(round(obsT)));
else
    legs = sprintf('%i min after',round(obsT));
end

fSz = 16;

% silence warnings
warns = {'curvefit:fit:complexYusingOnlyReal','MATLAB:Axes:NegativeDataInLogAxis', 'MATLAB:legend:IgnoringExtraEntries'};
for wI = 1:length(warns)
    st(wI) = warning('query', warns{wI});
    warning('off', warns{wI});
end

if doPlot
    prep_figure(fh, tits, fSz, yLs, length(dims), show_ints, norm_mode);
    clear h
end

allOCs = {};
FT = cell(length(dims),1);

for dimI = 1:length(dims)
    dim = dims(dimI);
    
    % Get the MSD
    if isscalar(dims)
        msdV = msdObj.getMeanMSD;
    else
        msdV = msdObj.msd{dim};
    end
    tau = msdV(2:end-nSkip,1);
    msd = msdV(2:end-nSkip,2)./msdNorm(dimI);
    
    % Do the lowpass
    if ~isempty(lpFrq)
        msd = lowpass_logspace(tau, msd,lpFrq);
        if lpFrq > 1
            idx1 = round(lpFrq/2);
        else
            idx1 = 1;
        end
    else
        idx1 = 1;
    end
    
    % Determine truncation point
    [idx, eta] = msd_truncator(tau, msd, trunc_mode, FF, truncT);
    
    % Extrapolate if necessary
    [tau, msd, eta, idx] = msd_extrapolator(tau, msd, idx, eta, extrap_mode, etaIn);
    
    if ~isempty(etaIn)
        warning('Manually overriding eta (1/long time MSD gradient) for FT')
        eta = etaIn;
    end
    % Do the interpolation and rheoFT
    [omega, G1, G2] = msd_interp_FT(tau(idx1:idx), msd(idx1:idx), 0, eta, idx, interpF); % Need a smart way of changing nPoints (idx)
    
    
    % Lowpass if there's a frequency to use
    %     if ~isempty(lpFrq)
    %         % (I've actually not written this yet)
    %         G1 = lowpass_logspace(omega, G1, lpFrq);
    %         G2 = lowpass_logspace(omega, G2, lpFrq);
    %     end
    
    FT{dimI} = [omega, G1, G2];
    
    % Find and show intercepts
    if doPlot && show_ints
        subplot(2+show_ints, length(dims),length(dims)* (1 + show_ints) + dimI)
    end
    oC = gstar_interceptor(omega, G1, G2, wR{dimI}, show_ints && doPlot, colour);
    
    allOCs{dim} = oC;
    
    % Plot MSD
    if doPlot
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
        h(end+1) = plot(nF.*tau([idx1 idx]), msd([idx1 idx]), ...
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
            'Color', 'b', 'Marker', 'x', 'LineStyle', 'none', 'LineWidth', 2);
        loglog(omega, ...
            G2, 'LineWidth', 2, ...
            'Color', 'r', 'Marker', 's', 'LineStyle', 'none', 'LineWidth', 2);
        
        % % If you wanted to do them custom
        %     loglog(omega, ...
        %         G1, 'LineWidth', 2, ...
        %         'Color', colour, 'Marker', mS, 'LineStyle', '-', 'LineWidth', 2);
        %     loglog(omega, ...
        %         G2, 'LineWidth', 2, ...
        %         'Color', colour, 'Marker', mS, 'LineStyle', '--', 'LineWidth', 2);
        
        % Show Intercept frequency without changing YLims
        yl = ylim;
        plot(oC.*[1; 1], ylim, '--','Color',0.7*[1 1 1 0.8], 'LineWidth', 2)
        ylim(yl);
        
        if showLeg
            legend('"Storage"','"Loss"','Intercept frequency','Location','best')
            legend(h, legs,'ω (min, max)','1 ÷ Intercept frequency','Location','best')
        end
    end
end

if nargout == 2
    varargout{1} = allOCs;
end

for wI = 1:length(warns)
    warning(st(wI).state, warns{wI});
end
% End of main function

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
            if ~isempty(yLs)
                xlim([1e-4 2e2])
                ylim(yLs)
            end
            
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

    function [idx, eta] = msd_truncator(tau, msd, mode, FF, truncT)
        % Tells you which idx to truncate MSD at (helpfully reuses gradient
        % calculation)
        dydxfilt = msd_gradientor(tau, msd);
        
        % Take minima
        [~, idx] = min(dydxfilt);
        
        % Fudge factor extra points to include in FT
        if ~exist('FF','var')
            FF = 40;
        elseif isempty(FF)
            FF = 40;
        end
        
        switch mode
            % Truncate at plateau/gradient minima + fudge factor
            case 'minima'
                idx = idx + FF;
                % Don't truncate
            case 'FF'
                idx = length(tau) - FF;
            case 'timeFF'
                idx = find(tau > truncT, 1);
            otherwise
                idx = length(tau);
        end
        if idx > size(dydxfilt,1)
            eta = 1./ dydxfilt(end);
        else
            eta = 1./ dydxfilt(idx);
        end
    end

    function [tau, msde, eta, idx] = msd_extrapolator(tau, msd, idx, eta, mode, etaIn)
        % Needs to extrapolate MSD, replacing some points with (log) linear fit for
        % the same tau values
        switch mode
            case 'linear'
                nP = 30;
                %         fo = fit(log(tau(idx-nP:idx)), log(msd(idx-nP:idx)), 'Poly1');
                
                [a, b] = leastSq(log(tau(idx-nP:idx)), log(msd(idx-nP:idx)));
                %
                if ~isempty(etaIn)
                    a = 1/etaIn;
                end
                fo = struct('p1',a,'p2',log(2*b));
                
                msde = msd;
                msde(idx-nP:end) = exp(fo.p1 * log(tau(idx-nP:end)) + fo.p2);
                
                eta = 1./fo.p1;
                idx = length(tau);
            otherwise
                msde = msd;
        end
        
    end

    function oC = gstar_interceptor(omega, G1, G2, wRange, doPlot, colour)
        % Find the intercepts between G1 and G2 as functions of frequency omega.
        % Uses linear fit over ranges given in cell array wRange.
        
        oC = nan(1,max(1,length(wRange)));
        
        if doPlot
            % Do it as points
            h(2) = loglog(omega, G1.\G2, 'Color', colour, 'LineWidth', 2, 'LineStyle', 'none', 'Marker', 'o');
            % % Do it as a line
            %     h(2) = loglog(omega, G1.\G2, 'Color', colour, 'LineWidth', 2);
            plot(xlim, [1 1], '--','Color',[1 1 1]*0.8, 'LineWidth', 3)
        end
        
        for wRIdx = 1:length(wRange)
            wdx = [0 0];
            wdx(2) = max(sum(omega > wRange{wRIdx}(1)), 1);
            wdx(1) = max(sum(omega > wRange{wRIdx}(2)), 1);
            
            odata = omega(wdx(1):wdx(2));
            Gdata = G1(wdx(1):wdx(2)).\G2(wdx(1):wdx(2));
            
            try
                %         fo = fit(log(odata), log(Gdata), 'Poly1');
                
                [a, b] = leastSq(real(log(odata)), real(log(Gdata)));
                %
                fo = struct('p1',a,'p2',log(2*b));
                
                oC(wRIdx) = exp(- fo.p1 \ fo.p2);
                
                if oC(wRIdx) > odata(1) || oC(wRIdx) < odata(end) % Frequency omega is high to low
                    oC(wRIdx) = nan;
                end
                
                if doPlot
                    h = plot(odata, exp(fo.p1 * log(odata) + fo.p2), 'g:', 'LineWidth', 2.5);
                end
            catch ME
                warning(ME.identifier, 'Fitting failed with error %s\n',ME.message)
                oC(wRIdx) = nan;
            end
        end
        if doPlot
            try
                if showLeg
                    legend(h, 'Fit result', 'Ratio G''''÷G''','Location','best')
                end
            catch
                try
                    if showLeg
                        legend(h(2), 'Ratio G''''÷G''','Location','best')
                    end
                catch ME
                    error(ME.identifier, 'No idea what happened in gstar_interceptor')
                end
            end
        end
    end
end