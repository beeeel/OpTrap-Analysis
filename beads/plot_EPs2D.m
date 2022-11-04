function [varargout] = plot_EPs2D(datOut, pIdxs, varargin)
%% Plot 2D scatters of 2 parameters, can output 1 parameter fit for viscosity
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
p.addRequired('pIdxs',@(x)validateattributes(x,{'numeric'},{'numel',2}))

p.addParameter('mult',[1 1],@(x)validateattributes(x,{'numeric'},{'numel',2,'positive'}))
p.addParameter('NF',[0 0],@(x)validateattributes(x,{'numeric'},{'numel',2}))
p.addParameter('xls',[],@(x)validateattributes(x,{'numeric'},{'numel',2,'increasing'}))
p.addParameter('yls',[],@(x)validateattributes(x,{'numeric'},{'numel',2,'increasing'}))
p.addParameter('dofit', true, @(x) isa(x,'logical'))
p.addParameter('holdOn', false, @(x) isa(x,'logical'))
p.addParameter('plotFit', true, @(x) isa(x,'logical'))

p.addParameter('colour', 'k', @(x)(isa(x,'char') && isscalar(x)) || (isa(x,'numeric') && all(x <= 1) && length(x) == 3) || strcmp(x,'scatter') || strcmp(x,'rainbow'))
p.addParameter('mark','o', @(x) any(strcmp(x,{'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram', 'none'})))

p.parse(datOut, pIdxs, varargin{:});

% Which endpoints to plot? (Column number in datCtl)
if length(pIdxs) ~= 2 || any(pIdxs > size(datOut,2))
    error('invalid pIdxs and lazy error checking')
end

mult = p.Results.mult;
NF = p.Results.NF;
    
RT = {'Radial','Tangential'};
if ~p.Results.holdOn
    fh = figure(99); % I got 99 problems and figure numbering ain't one
    clf
end
for dim = 1:2
    allDat = {datOut{1,pIdxs(1)}(:,dim+1), datOut{1,pIdxs(2)}(:,dim+1)};
    
    ax = subplot(1,2,dim);
    hold on

%     for idx = 1:2
%         if pIdxs(idx) == 1
%             mult(idx) = 1e17;
% %         elseif any(pIdxs(idx) == [6 7 9 11])
% %             mult(idx) = 1e-17;
%         elseif pIdxs(idx) == 3 && pIdxs(3-idx) == 2 % If we're plotting α against αH
%             mult(idx) = 0.5;
%             NF(idx) = -0.5;
%         end
%     end
    
    if ~strcmp(p.Results.colour, 'scatter') && ~strcmp(p.Results.colour, 'rainbow')
        h = plot(mult(1)*allDat{1}+NF(1),mult(2)*allDat{2}+NF(2), p.Results.mark, 'LineWidth', 2, 'Color', p.Results.colour, 'MarkerFaceColor', p.Results.colour);
    elseif strcmp(p.Results.colour, 'rainbow')
        Ms = {'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram'};
        Cs = 'rgbcm';
        l = 0;
        n = 1;
        hold on
        
        tDat = datOut{1,1}(:,1);
        % It's prettier but it might not work
        % Number of cells is last index
        nC = datOut{1,1}(end,4);
        % Get a list of which observation is which cell
        Cidx = datOut{1,1}(:,4);
        while n <= nC
            l = l + 1;
            m = 1;
            while m <= 5
                idxs = Cidx == n;
                plot(mult(1)*allDat{1}(idxs)+NF(1), mult(2)*allDat{2}(idxs)+NF(2), 'Color', Cs(m), 'Marker', Ms{l},'LineStyle','-')
                m = m + 1;
                n = n + 1;
            end
        end
    else
        h = scatter(mult(1)*allDat{1}+NF(1),mult(2)*allDat{2}+NF(2), [], datOut{1,1}(:,1), p.Results.mark, 'LineWidth', 2);
    end
    
    if p.Results.dofit
%         [a, b] = leastSq(log(mult(1)*allDat{1}+NF(1)),log(mult(2)*allDat{2}+NF(2)), ax);
        xData= log(mult(1)*allDat{1}+NF(1));
        xPlt = [min(xData) max(xData)];
        yData = log(mult(2)*allDat{2}+NF(2));
        N = isnan(xData) | isnan(yData);
        if any(strcmp(datOut(2,pIdxs), '(βG_0)^{-1}'))
            [fo, ~, out] = fit(xData(~N), yData(~N), @(a, x) a-x, 'Start', 1);
            error('You haven''t checked that the fitting equation is correct here')
            if p.Results.plotFit
                plot(exp(xPlt), exp(fo.a)./exp(xPlt)', 'k--')
            end
        elseif any(strcmp(datOut(2,pIdxs), 'βG_0'))
            [fo, ~, out] = fit(xData(~N), yData(~N), @(a, x) x-a, 'Start', 1);
            if p.Results.plotFit
                plot(exp(xPlt), exp(fo.a).*exp(xPlt)', 'k--')
            end
        end
        
        if nargout > 0
            as = [exp(fo.a+std(out.residuals)) exp(fo.a-std(out.residuals))] - exp(fo.a);
%         fprintf('\nConfidence interval for η: %g\t %g\tPa.s\n',as+exp(fo.a))
%         fprintf('\t\tError on η: %g Pa.s\n', mean(as))
            varargout{1}(dim) = exp(fo.a);
            if nargout > 1
                varargout{2}(dim) = mean(as);
            end
        end
    end
    
    if ~isempty(p.Results.xls)
        xlim(p.Results.xls);
    end
    if ~isempty(p.Results.yls)
        ylim(p.Results.yls);
    end
    
%     xl = xlim;
%     yl = ylim;
    xt = ax.XTick;
    yt = ax.YTick;
    xtl = ax.XTickLabel;
    ytl = ax.YTickLabel;
%     xlim(xl)
%     ylim(yl)


    
    if any(datOut{1,2} < 0)
        s = 'Δ';
    else
        s = '';
    end
    xlabel([s datOut{2,pIdxs(1)}])
    ylabel([s datOut{2,pIdxs(2)}])
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
    ax.XScale = 'log';
end