function [varargout] = plot_EPsTime(datOut, pIdxs, varargin)
%% Plot...
% [a] = plot_EPs2D(datOut, pIdxs, [varargin])
% varargin in name-value pairs:
%   mult        - multiplicative factor, 1 element per dimension/EP
%   NF          - additative factor, as above
%   xls         - xlimits for plot
%   yls         - ylimits for plot
%   holdOn      - Draw on (assumably) extant axes
%   mark        - Marker to plot
%   colour      - Colour to plot
%   dofit       - Plot a least-squares fit to data

p = inputParser;

p.addRequired('datOut',@(x)isa(x,'cell'));
p.addRequired('pIdxs',@(x)validateattributes(x,{'numeric'},{'numel',1,'<=',size(datOut,2)}))

p.addParameter('mult',1,@(x)validateattributes(x,{'numeric'},{'numel',1,'positive'}))
p.addParameter('NF',0,@(x)validateattributes(x,{'numeric'},{'numel',1}))
p.addParameter('xls',[],@(x)validateattributes(x,{'numeric'},{'numel',2,'increasing'}))
p.addParameter('yls',[],@(x)validateattributes(x,{'numeric'},{'numel',2,'increasing'}))
p.addParameter('dofit', true, @(x) isa(x,'logical'))
p.addParameter('holdOn', false, @(x) isa(x,'logical'))
p.addParameter('plotFit', true, @(x) isa(x,'logical'))

p.addParameter('colour', 'k', @(x)(isa(x,'char') && isscalar(x)) || (isa(x,'numeric') && all(x <= 1) && length(x) == 3) || any(strcmp(x,{'scatter','scat2'})))
p.addParameter('mark','o', @(x) any(strcmp(x,{'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram', 'none'})))

p.parse(datOut, pIdxs, varargin{:});

mult = p.Results.mult;
NF = p.Results.NF;
    
RT = {'Radial','Tangential'};
if ~p.Results.holdOn
    fh = figure(99); % I got 99 problems and figure numbering ain't one
    clf
else
    fh = gcf;
end
for dim = 1:2
    tDat = datOut{1,pIdxs}(:,1);
    allDat = mult * datOut{1,pIdxs}(:,dim+1) + NF;
    
    ax = subplot(1,2,dim);
    hold on
    
    % It's not pretty but it works
    % Number of cells is equal to number of times the observation time goes down
    nC = sum(tDat(1:end-1) > tDat(2:end))+1;
    % Use the cumulative sum to get a list of which observation is which cell
    Cidx = [1; cumsum(tDat(1:end-1) > tDat(2:end))+1];
    
    if strcmp(p.Results.colour, 'scatter')
        colourmap = colormap;
        colourmap = colourmap(round(linspace(1,size(colourmap,1),nC)),:);
        colormap(colourmap);
        scatter(tDat, allDat, [], Cidx, 'LineWidth',2, 'Marker',p.Results.mark)
%         if dim == 2
%             colorbar
%         end
    elseif strcmp(p.Results.colour, 'scat2')
        Ms = {'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram'};
        colourmap = colormap;
        colourmap = colourmap(round(linspace(1, size(colourmap,1), 5)), :);
        l = 0;
        n = 1;
        while n <= nC
            l = l + 1;
            m = 1;
            while m <= 5
                idxs = Cidx == n;
                plot(tDat(idxs), allDat(idxs), 'Color', colourmap(m, :), 'Marker', Ms{l},'LineStyle','none')
                m = m + 1;
                n = n + 1;
            end
        end
    else
        Ms = {'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram'};
        Cs = 'rgbcm';
        l = 0;
        n = 1;
        while n <= nC
            l = l + 1;
            m = 1;
            while m <= 5
                idxs = Cidx == n;
                plot(tDat(idxs), allDat(idxs), 'Color', Cs(m), 'Marker', Ms{l},'LineStyle','-')
                m = m + 1;
                n = n + 1;
            end
        end
    end
    
%     if p.Results.dofit
%         xData= tDat;
%         xPlt = [min(xData) max(xData)];
%         yData = log(allDat);
%         N = isnan(xData) | isnan(yData);
%         if any(strcmp(datOut(2,pIdxs), 'βG_0'))
%             fo = fit(xData(~N), yData(~N), @(a, x) a-x, 'Start', 1);
%             if p.Results.plotFit
%                 plot(exp(xPlt), exp(fo.a)./exp(xPlt)', 'k--')
%             end
%         elseif any(strcmp(datOut(2,pIdxs), '(βG_0)^{-1}'))
%             fo = fit(xData(~N), yData(~N), @(a, x) a+x, 'Start', 1);
%             if p.Results.plotFit
%                 plot(exp(xPlt), exp(fo.a).*exp(xPlt)', 'k--')
%             end
%         end
%         if nargout > 0
%             varargout{1}(dim) = exp(fo.a);
%         end
%     end
    
    if ~isempty(p.Results.xls)
        xlim(p.Results.xls);
    end
    if ~isempty(p.Results.yls)
        ylim(p.Results.yls);
    end
    
    if any(datOut{1,2} == 1)
        s = sprintf('%s/%s_1',datOut{2,pIdxs},datOut{2,pIdxs});
    else
        s = datOut{2,pIdxs};
    end
    xlabel('Time after bead adhesion (min)')
    ylabel(s)
    title(RT{dim})
    
%     ax.XTick = unique([xt(1) 0 xt(end)]);
%     ax.XTickLabel = unique([xtl(1) '0' xtl(end)]);
%     ax.YTick = unique([yt(1) 0 yt(end)]);
%     ax.YTickLabel = unique([ytl(1) '0' ytl(end)]);
    ax.FontSize = 18;
    
%     xline(0);
%     yline(0);
    
%     lh = legend(h, 'Control','Lat B');
%     lh.String = lh.String(1:2);
%     lh.Location = 'west';

    ax.YScale = 'log';
%     ax.XScale = 'log';
end