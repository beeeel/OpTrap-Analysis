function fh = bead_plotProData(data, varargin)
%% figureHandle = plotProBeadData(data, [setLims, figNum])
% Do histograms and scatterplots for 2D position/time data

xCentresM = data.pro.xCentresM;
yCentresM = data.pro.yCentresM;
timeVec = data.raw.timeVecMs;
cropT = data.opts.cropT;      

xStiff = calcStiffness(xCentresM);
yStiff = calcStiffness(yCentresM);
    
if nargin >= 2
    % help give better errors when misused!
    argN = 1;
    setLims = varargin{argN};
    validateattributes(setLims,{'numeric'},{'numel',2},...
        'bead_plotProData','setLims',nargin+argN-length(varargin));
    else
        setLims = [];
end

if nargin < 3
    fh = figure;
elseif nargin == 3
    fh = figure(varargin{2}); %#ok<*UNRCH>
else
    error('How many nargins did you use? Should be 1 to 3!')
end

fh.Name = data.fName;
clf

if size(xCentresM, 1) == 1
    % Histogram of xCentres and yCentres in units um
    subplot(3,1,1)
    hold on
    histogram(xCentresM.*1e6,'Normalization','probability')
    histogram(yCentresM.*1e6,'Normalization','probability')
    if ~isempty(setLims)
        xlim(setLims)
    end
    xlabel('Centre position (\mu m)')
    ylabel('Probability')
    order = num2str(data.opts.pOrder);
    title({['Histogram of centres, polynomial order ' order ' filtered '], [ 'trap stiffness kx = ' num2str(xStiff./1e-6) ' pN/\mu m, ky = ' num2str(yStiff./1e-6) ' pN/\mu m']})
    legend('X','Y')
    
    % Scatterplot of each centre in units um
    subplot(3,1,2)
    plot(xCentresM.*1e6,yCentresM.*1e6,'.')
    if ~isempty(setLims)
        xlim(setLims)
        ylim(setLims)
    end
    title(['Scatterplot of centres, ' num2str(length(xCentresM)) ' frames'])
    xlabel('X (\mu m)')
    ylabel('Y (\mu m)')
    axis equal
    
    % Time traces of X and Y in units um
    subplot(3,2,5)
    plot(1e-3.*timeVec(cropT(1):cropT(2)), xCentresM.*1e6,'.')
    xlabel('Time (s)')
    ylabel('X (\mu m)')
    title('X time trace')
    if ~isempty(setLims)
        ylim(setLims)
    end
    subplot(3,2,6)
    plot(1e-3.*timeVec(cropT(1):cropT(2)), yCentresM.*1e6,'.')
    xlabel('Time (s)')
    ylabel('Y (\mu m)')
    title('Y time trace')
    if ~isempty(setLims)
        ylim(setLims)
    end
else
    for obj = 1:size(xCentresM,1)
        % Histogram of xCentres and yCentres in units um
        subplot(3,1,1)
        hold on
        histogram(xCentresM(obj,:).*1e6,'Normalization','probability')
        histogram(yCentresM(obj,:).*1e6,'Normalization','probability')
        if ~isempty(setLims)
            xlim(setLims)
        end
        xlabel('Centre position (\mu m)')
        ylabel('Probability')
        order = num2str(data.opts.pOrder);
        title({['Histogram of centres, polynomial order ' order ' filtered ']});%, [ 'trap stiffness kx = ' num2str(xStiff./1e-6) ' pN/\mu m, ky = ' num2str(yStiff./1e-6) ' pN/\mu m']})
        legend('X left','Y left','X right','Y right')
        
        % Scatterplot of each centre in units um
        subplot(3,1,2)
        hold on
        plot(xCentresM(obj,:).*1e6,yCentresM(obj,:).*1e6,'.')
        if ~isempty(setLims)
            xlim(setLims)
            ylim(setLims)
        end
        title(['Scatterplot of centres, ' num2str(length(xCentresM)) ' frames'])
        xlabel('X (\mu m)')
        ylabel('Y (\mu m)')
        axis equal
        
        % Time traces of X and Y in units um
        subplot(3,2,5)
        hold on
        plot(1e-3.*timeVec(cropT(1):cropT(2)), xCentresM(obj,:).*1e6,'.')
        xlabel('Time (s)')
        ylabel('X (\mu m)')
        title('X time trace')
        if ~isempty(setLims)
            ylim(setLims)
        end
        
        subplot(3,2,6)
        hold on
        plot(1e-3.*timeVec(cropT(1):cropT(2)), yCentresM(obj,:).*1e6,'.')
        xlabel('Time (s)')
        ylabel('Y (\mu m)')
        title('Y time trace')
        if ~isempty(setLims)
            ylim(setLims)
        end
    end
end
drawnow    
