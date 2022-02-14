function msd_colorfulplot(accumulated, dIs, cIs, rIs, varargin) %#ok<INUSL>
%% msd_colorfulplot(accumulated, dIs, cIs, rIs, ... )
% Plot MSDs using nice colourmap according to time. I will add support for
% name-value pair arguments to help this be good, but not yet.
%% TO DO
% Name-value pair inputs
% Fix normalisation, including supply of normalization times/frequencies
if nargin >= 6
    if isa(varargin{1}, 'matlab.graphics.axis.Axes') && isa(varargin{2}, 'matlab.graphics.axis.Axes')
        axs = [varargin{1}, varargin{2}];
    end
end

if exist('rIs','var')
elseif ~exist('rIs','var') && isscalar(dIs) && isscalar(cIs)
    rIs = 1:size(accumulated{dIs}{1,cIs},2);
else
    error('Either dIs and cIs must be scalar, or you must provide rIs')
end
% Assume sample temp = room temp = 22 C
T = 22;
kbT = (272+T).*1.38064852e-23;
% Viscosity of water at 22 C
eta =  0.9544e-3; % need to copy equation for eta(T)
% Bead radius is ~2.5um (maybe as much as 3)
a = 2.5e-6;

% Ignore the last n points of MSD
nSkip = 40;
% LineStyles - default to solid line
lsty = {'-','-','-','-','-','-'}; 
% MarkerStyles - default to solid line
msty = {'.','.','.','.','.','.'}; 
% Color scale factors - multiply RGB triplet by this
cF = [1 1 1 1 1 1];
% Normalize the MSDs? (use MSDp to normalize)
normY = false;
normX = false;

if ~exist('axs','var')
    % Setup figure window
    tits = {'Radial', 'Tangential'};
    figure(2)
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
allTs = sort([accumulated{dIdx}{2,:}]);

for cIdx = cIs
    
   MSDs = [accumulated{dIdx}{1,cIdx}.msd];
    ts = accumulated{dIdx}{2,cIdx};
    
    colormap cool
    colourmap = colormap;
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
        [dydx, tout] = msd_gradientor(tau, msd, 'lsq');
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
        if normX
            lambda(:,rIdx, cIdx) = kbT./ (p * 3 * pi * eta * a);
            lambda(1,rIdx, cIdx) = allOCs{1,rIdx}(1);
            lambda(2,rIdx, cIdx) = allOCs{2,rIdx}(1);
        else
            lambda(:, rIdx, cIdx) = [1;1];
        end
        if normY
            normF = p;
        else
            normF = [1 1];
        end
        
        
        for plt = 1:2
            h(pC) = loglog(axs(plt), lambda(plt,rIdx, cIdx).*MSDs(rIdx).msd{plt}(1:end-nSkip,1), ...
                MSDs(rIdx).msd{plt}(1:end-nSkip,2)./normF(plt), 'LineWidth', 2, ...
                'Color', colour(rIdx,:).*cF(dIdx), 'LineStyle', lsty{dIdx}, ...
                'Marker', msty{dIdx});
            h(1) = plot(axs(plt), lambda(plt,rIdx, cIdx).*tau(idx(plt)), p(plt)./normF(plt), 'kx','LineWidth', 2);
%             h(1) = plot([1 1], ylim, '--', 'Color', 0.8 * [1 1 1 0.8],...
%                 'LineWidth', 3);
        end
        pC = pC + 1;
    end
end
end
legend(h, legCell, 'Location','best')