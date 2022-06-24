function msd_colorfulplot(accumulated, dIs, cIs, rIs, tRanges, varargin) %#ok<INUSL>
%% msd_colorfulplot(accumulated, dIs, cIs, rIs, ... )
% Plot MSDs using nice colourmap according to time. I will add support for
% name-value pair arguments to help this be good, but not yet.
% You can either supply two axis for the plots to go on (radial and
% tangential), or you can supply normX and normY
%% TO DO
% Name-value pair inputs
% Fix normalisation, including supply of normalization times/frequencies
if exist('rIs','var') && ~isempty(rIs)
elseif (~exist('rIs','var') || isempty(rIs)) && isscalar(dIs) && isscalar(cIs)
    rIs = 1:size(accumulated{1,dIs}{1,cIs},2);
else
    error('Either dIs and cIs must be scalar, or you must provide rIs')
end

if nargin >= 7
    if isa(varargin{1}, 'matlab.graphics.axis.Axes') && isa(varargin{2}, 'matlab.graphics.axis.Axes')
        axs = [varargin{1}, varargin{2}];
    end
    if isnumeric(varargin{1})
        validateattributes(varargin{1}, {'numeric'}, {'nrows',2,'ncols',max(rIs)},'msd_colorfulplot','time domain normalisation',5)
        lambda = varargin{1};
        normX = true;
    else
        normX = false;
    end
    if islogical(varargin{2})
        normY = varargin{2};
    end
end

% Assume sample temp = room temp = 22 C
T = 22;
kbT = (272+T).*1.38064852e-23;
% Viscosity of water at 22 C
eta =  0.9544e-3; % need to copy equation for eta(T)
% Bead radius is ~2.5um (maybe as much as 3)
a = 2.5e-6;

% Number of points for gradient calculation
nP = 25;

% Ignore the last n points of MSD
nSkip = 40;
% LineStyles - default to solid line
lsty = repmat({'-'},1,20);
% MarkerStyles - default to solid line
msty = repmat({'.'},1,20);
% Color scale factors - multiply RGB triplet by this
cF = ones(1,20);
% Normalize the MSDs? (use MSDp to normalize)
if ~exist('normY','var')
    normY = false;
end
if ~exist('normX','var')
    normX = false;
end

if ~exist('axs','var')
    % Setup figure window
    tits = {'Radial', 'Tangential'};
    figure(2+normY+2*normX)
    clf
    for plt = 1:2
        axs(plt) = subplot(1,2,plt);
        hold on
        set(gca,'XScale','log')
        set(gca,'YScale','log')
        %     axis equal
        if normX
            xlabel('Normalized time, λτ')
            %         xlabel('Normalized time, τω_c')
        else
            xlabel('Lag time \tau (s)')
        end
        if normY
            ylabel('normalized MSD')
            ylim([5e-4 5e2])
        else
            ylabel('MSD (μm^2)')
            %         ylim([8e-6 1e-1])
            ylim([5e-6 3e-1])
        end
        title(tits{plt})
        set(gca,'FontSize',15)
    end
end
clear h
pC = 2;
legCell = {'Gradient minima'};
% legCell = {'1 ÷ ω_c'};
for dIdx = dIs
    % msdps =
    allTs = sort([accumulated{1,dIdx}{2,:}]);
    
    for cIdx = cIs
        
        MSDs = [accumulated{1,dIdx}{1,cIdx}.msd];
        ts = accumulated{1,dIdx}{2,cIdx};
        
        %         colormap cool
        %         colourmap = colormap;
        %     min(ceil(size(colourmap,1)*(ts-ts(1)+1e-6)/(allTs(end) - allTs(1)+1e-6)))
        %     max(ceil(size(colourmap,1)*(ts-ts(1)+1e-6)/(allTs(end) - allTs(1)+1e-6)))
        
        nNeg = sum(ts < 0);
        colormap winter
        colourmap = colormap;
        colour = colourmap(16+ceil(0.75*size(colourmap,1)*(ts(1:nNeg)-ts(1)+1e-6)/(abs(ts(1))+1e-6)),:);
        colormap autumn
        colourmap = flipud(colormap);
        colour = [colour; colourmap(ceil(size(colourmap,1)*(ts(nNeg+1:end)+1e-6)/(ts(end)+1e-6)),:)];
        
        %     legCell = {};
        for rIdx = rIs
            legCell{pC} = sprintf('%.2g min', ts(rIdx));
            % Extract all data from msd field into 3D array, then into 2D array
            MSDp = cat(3, MSDs(rIdx).msd{:});
            tau = MSDp(1:end-nSkip,1,1);
            msd = squeeze(MSDp(1:end-nSkip,2,:));
            % Calculate derivative
            [dydx, tout] = msd_gradientor(tau, msd, 'lsq', nP);
            % Take minima
            [~, idx] = min(dydx);
            tmin = tout(idx);
            tmp = tau > tmin';
            tmp = tmp(2:end,:) > tmp(1:end-1,:);
            idx = find(tmp, 2);
            ind = idx;
            idx = ind - [0; size(tau,1)];
            p = msd(ind);
            cTau(:,rIdx) = tau(idx);
            if normX && ~exist('lambda','var')
                lambda(:, rIdx) = [1;1];
                %             lambda(:,rIdx, cIdx) = kbT./ (p * 3 * pi * eta * a);
                %             lambda(1,rIdx, cIdx) = allOCs{1,rIdx}(1);
                %             lambda(2,rIdx, cIdx) = allOCs{2,rIdx}(1);
                lambda(:,1:2) = 1./[0.0122, 0.0805; 0.0089, 0.0334];
            elseif ~exist('lambda','var')
                lambda = repmat([1;1],1,max(rIs));
            end
            if normY
                normF = p;
            else
                normF = [1 1];
            end
            
            
            for plt = 1:2
                h(pC) = loglog(axs(plt), ...
                    lambda(plt,rIdx).*MSDs(rIdx).msd{plt}(1:end-nSkip,1), ...
                    MSDs(rIdx).msd{plt}(1:end-nSkip,2)./normF(plt), 'LineWidth', 2, ...
                    'Color', colour(rIdx,:).*cF(dIdx), 'LineStyle', lsty{dIdx}, ...
                    'Marker', msty{dIdx});
                h(1) = plot(axs(plt), lambda(plt,rIdx).*tau(idx(plt)), p(plt)./normF(plt), 'kx','LineWidth', 2);
                %             h(1) = plot([1 1], ylim, '--', 'Color', 0.8 * [1 1 1 0.8],...
                %                 'LineWidth', 3);
                if exist('tRanges','var') && size(tRanges,1) >= rIdx
                    taui = logspace(log10(tau(2)), log10(tau(end)), length(tau).*1e3)';
                    msdi = interp1(tau, msd(:,plt), taui, 'pchip');
                    for fIdx = 1:length(tRanges{rIdx,plt})
                        % Get the indexes for the specified time range
                        msdIdx = find(taui > tRanges{rIdx,plt}{fIdx}(1), 1) ...
                            : find(taui < tRanges{rIdx,plt}{fIdx}(2), 1, 'last');
                        tauData = taui(msdIdx).*lambda(plt,rIdx);
                        msdData = msdi(msdIdx)./normF(plt);
                        
                        % Fit to log data (and get errors)
                        [fps(:,fIdx, dIdx), fitErr(:,fIdx, dIdx)] = N_get_fits;
                        % % Errors now come from fit function above
                        %             % Calculate fitErr
                        %             fitErr(:,fIdx, dIdx) = N_get_RMSE;
                        % Show fits
                        plot(axs(plt), tauData, exp(fps(1,fIdx, dIdx) * log(tauData) + log(fps(2,fIdx, dIdx))) , ...
                            'k:', 'LineWidth', 2.5)
                    end
                end
            end
            pC = pC + 1;
        end
    end
end
legend(h, legCell, 'Location','best')

    function [fps, fitErr] = N_get_fits
        % Either use linear fit or the least squares estimator from [1]Ling
        % 2019 eq. 2.4
        est = 'lsq'; % not at all hacky
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
end