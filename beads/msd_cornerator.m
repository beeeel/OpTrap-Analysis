function varargout = msd_cornerator(msdObj, obsT, tRanges, varargin)
%% [cTau, fitParams, fitErr] = msd_cornerator(msdObj, obsT, tRanges, varargin)
% Find forner times between tRanges using linear fits to loglog data.
% Accepts additional parameters in name-value pairs. Possible options:
% nSkip         - Number of points from MSD to skip
% dims          - Indices within msdObj.msd to use
% yLims         - Y limits when plotting MSDs
% figHand       - Figure handle to plot upon
% lineColour    - Line colour to plot MSD
% lineStyle     - Line style to plot MSD
% marker        - Marker to plot MSD
% interpM       - Interpolation method for fitting
% interpF       - Interpolation factor for fitting
% doPlot        - Exactly what it says.

%% TO DO:
% Consider whether it's worth making bead_processor_v2 with detailed
%  analysis (should mostly be copypasta from accumulator_v1)... probably
%  not
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
p.addParameter('doPlot', true, @(x) validateattributes(logical(x), {'logical'},{'scalar'}))
p.addParameter('figHand', [], @(x)isa(x,'matlab.ui.Figure'))
p.addParameter('lineColour', 'k', @(x)(isa(x,'char') && isscalar(x)) || (isa(x,'numeric') && all(x < 1) && length(x) == 3))
p.addParameter('lineStyle', '-', @(x) any(strcmp(x,{'-',':','-.','--','none'})))
p.addParameter('marker','none', @(x) any(strcmp(x,{'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram', 'none'})))
p.addParameter('interpM', 'pchip', @(x) any(strcmp(x, {'linear', 'nearest', 'next', 'previous', 'spline', 'pchip', 'cubic', 'v5cubic', 'makima'})))
p.addParameter('interpF', 1e2, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}))
p.addParameter('estimator', 'lsq', @(x) any(strcmp(x,{'lsq', 'fit'})))
p.addParameter('normT', [1 1], @(x)isa(x,'double') && length(x) == nMSDs && all(x < msdObj.msd{1}(1,end)) && all(x > 0))

p.parse(msdObj, obsT, tRanges, varargin{:});

% Intercept
tRanges = p.Results.tRanges;
interpM = p.Results.interpM;
interpF = p.Results.interpF;
% MSD options
nSkip = p.Results.nSkip;
dims = p.Results.dims;
% Estimator
est = p.Results.estimator;
% Plot options
doPlot = p.Results.doPlot;
yl = p.Results.yLims;
fh = p.Results.figHand;
colour = p.Results.lineColour;
lS = p.Results.lineStyle;
mS = p.Results.marker;
normT = p.Results.normT;
%% Setup

cTau = nan(round(max(length(tRanges{1}),length(tRanges{2}))/2),length(dims));

% % Just give up if there's no tRanges to work on
% if ~any([length(tRanges{1}), length(tRanges{2})] > 1)
%     varargout{1} = cTau;
%     if nargout > 1
%         varargout{2} = nan(size(cTau));
%         return
%     end
% end

legs = {};
if doPlot
    N_setup_fig;
end

fps = zeros(2, length(tRanges), length(dims));
fitErr = zeros(2, length(tRanges), length(dims));

% For each dimension 
for dIdx = dims
    d = dims(dIdx);
    subplot(length(dims)/2,2,dIdx)
    try
        % get the tau and msd, interpolate 
        tau = msdObj.msd{d}(2:end-nSkip,1);
        msd = msdObj.msd{d}(2:end-nSkip,2);
        
        taui = logspace(log10(tau(1)), log10(tau(end)), length(tau).*interpF)';
        msdi = interp1(tau, msd, taui, interpM);
        
        if doPlot
            h = plot(taui./normT(d), msdi, ...
                'Color',colour, 'LineWidth', 2, 'LineStyle', lS, 'Marker', mS);
        end
        
        % Fit to interpolated data from all the tRanges
        
        for fIdx = 1:length(tRanges{d})
            % Get the indexes for the specified time range
            msdIdx = find(taui > tRanges{d}{fIdx}(1), 1) ...
                : find(taui < tRanges{d}{fIdx}(2), 1, 'last');
            tauData = taui(msdIdx)./normT(d);
            msdData = msdi(msdIdx);
            
            % Fit to log data (and get errors)
            [fps(:,fIdx, dIdx), fitErr(:,fIdx, dIdx)] = N_get_fits;
            % % Errors now come from fit function above
            %             % Calculate fitErr
            %             fitErr(:,fIdx, dIdx) = N_get_RMSE;
            % Show fits
            if doPlot
                plot(tauData, exp(fps(1,fIdx, dIdx) * log(tauData) + log(fps(2,fIdx, dIdx))) , ...
                    'r:', 'LineWidth', 2.5)
            end
        end
        % Find corners by intercept of fits. I fucking hope I wrote this in
        % my lab book.... (I probably didn't)
        dfps = diff(fps,1,2);
        dfps(2,:,:) = diff(log(fps(2,:,:)),1,2);
        cTau(:, dIdx) = exp(-dfps(2,:,dIdx)./dfps(1,:,dIdx));
        
        if doPlot
            for corn = 1:size(cTau,1)
                h(2) = plot(cTau(corn,dIdx) * [1 1], yl, 'color', [0 0 0 0.25], 'LineWidth', 3);
            end
        end
    catch ME
        warning(['caught error: ' ME.message])
    end
    
    % Set legend
    if doPlot
        try
%             legend(h,{legs 'Fit intercept time'}, 'Location', 'northwest');
            legend('MSD',sprintf('α = %.2g, D = %.2gμm^2/s',fps(:,1,dIdx)),sprintf('α = %.2g, D = %.2gμm^2',fps(:,2,dIdx)))
        catch ME
            warning(['Couldn''t set legend: ' ME.message])
        end
    end
end

varargout{1} = cTau;
if nargout > 1
    varargout{2} = fps;
end
if nargout > 2
    varargout{3} = fitErr;
end

%% Function definitions

    function fitErr = N_get_RMSE
        % hahaha they'll never suspect that this isn't reall the RMSE
        err = msdData - exp(fps(1,fIdx, dIdx) * log(tauData) + fps(2,fIdx, dIdx));
%         fitErr = sqrt(mean(err.^2,'all'));
        fitErr = std(err);
    end

    function [fps, fitErr] = N_get_fits
        % Either use linear fit or the least squares estimator from [1]Ling
        % 2019 eq. 2.4
        switch est
            case 'lsq'
                [alpha, D] = leastSq(log(tauData), log(msdData));
                % This needs to output the gradient (alpha) and y-intercept
                % of the linear fit, whereas 2D is the MSD at τ=1. [1]
                fps = [alpha; (2*D)];
                % Truth is I couldn't tell you how this works, but it does.
                
                % hahaha they'll never suspect that this isn't reall the RMSE
                err = msdData - exp(fps(1) * log(tauData) + fps(2));
                %         fitErr = sqrt(mean(err.^2,'all'));
                fitErr = std(err); % So this is actually the standard devi
            case 'fit'
                % Fit to log data
                fo = fit(log(tauData), log(msdData), 'poly1');
                % Store the data for later
                fps = [fo.p1; fo.p2];
                fitErr = diff(confint(fo))';
        end
        
    end


    function N_setup_fig
        % Titles for X and Y directions
        if length(dims) == 2
            if obsT ~= 0
                tits = {'Radial','Tangential'};
            else
                tits = {'X','Y'};
            end
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
        
        for pIdx = dims
            plt = dims(pIdx);
            
            ax = subplot(length(dims)/2,2,plt);
            
            hold on
            clear h
            
            % This is horrible
            % Less so now
            ax.XScale = 'log';
            ax.YScale = 'log';
            
            ax.FontSize = fSz;
            
            [xl1, xl2] = bounds(msdObj.msd{1}(:,1));
            
            xlim([xl1 xl2])
            ylim(yl)
            
            xlabel('Delay time \tau (s)')
            ylabel('MSD (\mu m^2)')
            
            title(sprintf('%s direction',tits{plt}))
            
            grid on
            
        end
    end
end