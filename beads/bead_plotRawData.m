function fh = bead_plotRawData(data, varargin)
%% figureHandle = plotRawBeadData(data, [setLims, figNum])
% Do histograms and scatterplots for 2D position/time data

%% Setup
if ~isfield(data, 'mPerPx')
    error('No pixel calibration in data struct')
end

% Get time and cropping
timeVec = data.raw.timeVecMs;
if length(data.opts.cropT) == 2
    cropT = data.opts.cropT;      
else
    cropT = [1 length(timeVec)];
end

% Parse varargins
if nargin >= 2
    setLims = varargin{1};
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

% Choose centres row
if isfield(data.opts, 'centresRow')
    cRow = data.opts.centresRow;
elseif isfield(data.raw, 'suffixes')
    cRow = 1:length(data.raw.suffixes);
else
    cRow = 1:size(data.raw.xCentresPx,1);
end

% Get centres and demean
xCentresM = data.raw.xCentresPx(cRow,cropT(1):cropT(2));
yCentresM = data.raw.yCentresPx(cRow,cropT(1):cropT(2));
zCentresM = data.raw.zCentresPx(cRow,cropT(1):cropT(2));

% Calculate stiffness using equipartition method
xStiff = calcStiffness(xCentresM);
yStiff = calcStiffness(yCentresM);
zStiff = calcStiffness(zCentresM);

%% Plotting
fh.Name = data.fName;
clf

% Histogram of xCentres and yCentres in units um
if size(xCentresM, 1) == 1
    subplot(3,1,1)
    hold on
    histogram(xCentresM.*1e6,'Normalization','probability','DisplayName','X')
    histogram(yCentresM.*1e6,'Normalization','probability','DisplayName','Y')
    if isfield(data.raw,'zCentresPx'); histogram(zCentresM.*1e6,'Normalization','probability','DisplayName','Z'); end
    if ~isempty(setLims)
            xlim(setLims)
    end
    xlabel('Centre position (\mu m)')
    ylabel('Probability')
    order = num2str(data.opts.pOrder);
    str = [ 'trap stiffness \kappa_x = ' num2str(xStiff(1)./1e-6,3) ' N/Arb.U., \kappa_y = ' num2str(yStiff(1)./1e-6,3) ' N/Arb.U.']; if isfield(data.raw, 'zCentresPx'); str = [str ', \kappa_z = ' num2str(zStiff(1),3) ' N/Arb.U.']; end
    title(['Histogram of centres, unfiltered'], str)
    legend()
    
    % Scatterplot of each centre in units um
    subplot(3,1,2)
    if ~isfield(data.raw,'zCentresPx')
        plot(xCentresM,yCentresM,'.')
    else
        scatter3(xCentresM,yCentresM,zCentresM,[],data.raw.timeVecMs*1e-3)
        zlabel('Z (Arb. U.)')
    end

    if ~isempty(setLims)
        xlim(setLims)
        ylim(setLims)
    end
    title(['Scatterplot of centres, ' num2str(length(xCentresM)) ' frames'])
    xlabel('X (Arb. U.)')
    ylabel('Y (Arb. U.)')
    axis equal
    
    % Time traces of X and Y in units um
    subplot(3,2,5)
    plot(1e-3.*timeVec(cropT(1):cropT(2)), xCentresM.*1e6,'.')
    xlabel('Time (s)')
    ylabel('X (Arb. U.)')
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
        subplot( 3, size(xCentresM,1), obj)
        hold on
        histogram(xCentresM(obj,:).*1e6,30,'Normalization','probability')
        histogram(yCentresM(obj,:).*1e6,30,'Normalization','probability')
        if setLims
            xlim(setLims)
        end
        xlabel('Centre position (\mu m)')
        ylabel('Probability')
        title({'Histogram of centres, unfiltered' [ 'trap stiffness kx = ' num2str(xStiff(1)./1e-6) ' pN/\mu m, ky = ' num2str(yStiff(1)./1e-6) ' pN/\mu m']})
        warning('off','MATLAB:legend:IgnoringExtraEntries') % Probably not needed
        legend('X','Y')
        warning('on','MATLAB:legend:IgnoringExtraEntries')% Probably not needed
        
        % Scatterplot of each centre in units um
        subplot( 3, size(xCentresM,1), size(xCentresM,1) + obj)
        hold on
        scatter(xCentresM(obj,:).*1e6,yCentresM(obj,:).*1e6,[],1e-3.*timeVec(cropT(1):cropT(2)),'.')
        if setLims
            xlim(setLims)
            ylim(setLims)
        end
        title(['Scatterplot of centres, ' num2str(length(xCentresM)) ' frames'])
        xlabel('X (\mu m)')
        ylabel('Y (\mu m)')
        axis equal
        if obj == 1; colorbar; end
        
        % Time traces of X and Y in units um
        subplot(3,2,5)
        hold on
        plot(1e-3.*timeVec(cropT(1):cropT(2)), xCentresM(obj,:).*1e6,'.','MarkerSize',3)
        xlabel('Time (s)')
        ylabel('X (\mu m)')
        title('X time trace')
        if setLims
            ylim(setLims)
        end
        
        subplot(3,2,6)
        hold on
        plot(1e-3.*timeVec(cropT(1):cropT(2)), yCentresM(obj,:).*1e6,'.','MarkerSize',3)
        xlabel('Time (s)')
        ylabel('Y (\mu m)')
        title('Y time trace')
        if setLims
            ylim(setLims)
        end
    end
end
drawnow
