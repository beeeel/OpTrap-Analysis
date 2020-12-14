%% Confocal bead analysis
% Does "basic" analysis for each file- CofM, draw histogram and graph,
% attempt to measure width. At end is attempt to do autocorrelation/MSD
% measure

%% Load some confocal bead data
fNames = {'laser-100-50kfr_40-50kfr_too_many_px.czi', 'laser-25-100kfr_few_px.czi',...
    'laser-25-100kfr_more_px_bead_leaves_trap.czi','laser-25-100kfr_more_px_reducing_laser_power_fast.czi'};
cropTs = {[1 134078],[1 84440],[1 65e3],[1 92e3]};
cropXs = {[40 60],[50 70],[15 50], [15 50]};
for fileNameIdx = 1
    fileName = fNames(4);
    start = tic;
    Imstack = bfopen(['~/Documents/data/Zeiss/2020-11-10-2um-beads-repaired-XY-files/' fileName{:}]);
    imTime = toc(start);
    Metastack = bfmeta(['~/Documents/data/Zeiss/2020-11-10-2um-beads-broken-XT-files/' fileName{:}]);
    metaTime = toc(start) - imTime;
    fprintf('Time to load images: %f.0\nTime to load metadata: %f.0\n',imTime,metaTime)
    %% Do the processing
    % Metadata handling
    omeMeta = Metastack{4};
    hashTable = Metastack{2};
    stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
    stackSizeT = omeMeta.getPixelsSizeT(0).getValue(); % stack height, lines
    
    voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.METER); % in m
    voxelSizeXdouble = voxelSizeX.doubleValue();                                  % The numeric value represented by this object after conversion to type double
    voxelSizeTdouble = str2double(hashTable.get('Global Information|Image|Channel|LaserScanInfo|LineTime #1'));
    
    dwellTimedouble = str2double(hashTable.get('Global HardwareSetting|ParameterCollection|PixelDwellTime #1'));
    
    directionality = hashTable.get('Global HardwareSetting|ParameterCollection|ScanDirection #1');
    nDirections = strcmp(directionality,'Bidirectional') + 1;
    
    scanSpeed = voxelSizeXdouble./(dwellTimedouble.*1e-6)
    beadTime = 2e-6./scanSpeed
    
    % Image processing
    if size(Imstack{1},1) == 1
        Ims = Imstack{1}{1,1};
    else
        warning('untested line of code ahead. Check Ims variable is what it should be');
        Ims = cat(2,Imstack{1}{:,1});
    end
    
    % Crop to just the bead
    cropX = cropXs{fileNameIdx};
    cropT = cropTs{fileNameIdx};
    ImsC = zeros((diff(cropT)+1)/nDirections, diff(cropX)+1,nDirections);
    if nDirections == 2
        ImsC(:,:,1) = Ims(cropT(1):2:cropT(2), cropX(1):cropX(2));
        ImsC(:,:,2) = Ims(cropT(1)+1:2:cropT(2), cropX(1):cropX(2));
    else
        ImsC = Ims(cropT(1):cropT(2), cropX(1):cropX(2));
    end
    % Calculate centres in units of m (determined by unit set getting
    % voxelSizeX)
    ImsT = double(permute(ImsC,[2 1 3]));
    X = voxelSizeXdouble.*(1:size(ImsT,1))';
    Centres = sum(X.*ImsT)./sum(ImsT);
    if nDirections == 2
        centresIL(1:2:diff(cropT)+1) = Centres(:,:,1);
        centresIL(2:2:diff(cropT)+1) = Centres(:,:,2);
    else
        centresIL = Centres;
    end
    
    % Calculate bead widths in units of um
    thresholded = abs(ImsT) > 50;
    % If D(i+1) > D(i) you're in a new block of changing position
    shiftedD = zeros(size(thresholded),'logical');
    shiftedD(2:end,:) = thresholded(1:end-1,:);
    beadStart = shiftedD < thresholded;
    % If D(i-1) > D(i) you're in a new block of equilibrium
    shiftedU = zeros(size(thresholded),'logical');
    shiftedU(1:end-1,:) = thresholded(2:end,:);
    beadEnd = shiftedU < thresholded;
    % max returns the first true from each column
    [~, beadsStartX]= max(beadStart);
    [~, beadsEndX] = max(flipud(beadEnd));
    % Blocks ready for output
    beadsX = [beadsStartX; beadsEndX];
    beadWidthsX = diff(beadsX);
    
    % Calculate stiffness in units N/m (or uN/um)
    Kb = 1.38064852e-23; % Boltzmann constant
    T = 273 + 20; % Assume trap at room temperature - maybe I should model this?
    % From Sarshar et al 2014 eq 4. Need to convert CofM from px to m, so this
    % gives stiffness in N/m.
    varCentres = var(Centres,0,2);
    varStart = var(beadsStartX,0,2);
    varEnd = var(beadsEndX,0,2);
    if nDirections == 2
        Stiffness = Kb .* T ./ sqrt(0.5 * (varCentres(1).^2 + varCentres(2).^2));
        StiffnessL = Kb .* T ./ sqrt(0.5 * (varStart(1).^2 + varStart(2).^2));
        StiffnessR = Kb .* T ./ sqrt(0.5 * (varEnd(1).^2 + varEnd(2).^2));
    else
        Stiffness = Kb .* T ./ varCentres;
        StiffnessL = Kb .* T ./ varStart;
        StiffnessR = Kb .* T ./ varEnd;
    end
    %% Plot cropped data, histogram bead position and time series bead position
    fh = figure;
    fh.Name = fileName{:};
    clf
    m = 3;
    
    subplot(m,1,1) % Cropped bead images
    imagesc([0 voxelSizeTdouble*diff(cropT)],1e6.*[X(1) X(end)],Ims(cropT(1):cropT(2), cropX(1):cropX(2))')
    set(gca,'YDir','normal')
    xlabel('Time(s)')
    ylabel('Scan Axis (\mu m)')
    title(['Cropped data (' num2str(diff(cropT)+1) ' frames)'])
    
    subplot(m,1,2) 
    if strcmp(directionality,'Bidirectional')
        % Histogram of centres from two scan directions
        histogram(1e6.*reshape(Centres(:,:,1),1,[],1))
        hold on
        histogram(1e6.*reshape(Centres(:,:,2),1,[],1))
        ylabel('Count')
        xlabel('Position (\mu m)')
        title('Bead position distribution (bidirectional scan)')
        legend('Scan direction 1','Scan direction 2','Location','best')
    else
        histogram(Centres)
        ylabel('Count')
        xlabel('Position (\mu m)')
        title('Bead position distribution (unidirectional scan)')
    end
    
    subplot(m,1,3) % Centre position against time
    plot((cropT(1)-1:cropT(2)-1)*voxelSizeTdouble,1e6.*centresIL(cropT(1):cropT(2)),'k.')
    xlabel('Time(s)')
    ylabel('Position (\mu m)')
    xlim(cropT*voxelSizeTdouble)
    title(['Bead positions (Trap stiffness ' num2str(mean(Stiffness)*1e6) 'pN/um)'])
    %% Plot position of centre, LHS and RHS, and histogram width of bead
    m = 2;
    fh = figure;
    fh.Name = fileName{:};
    clf
    subplot(m,1,1)
    hold on
    if nDirections == 2
        xData = (cropT(1)-1:cropT(2)-1) * voxelSizeTdouble;
        yData = zeros(diff(cropT)+1,1);
        yData(1:2:end) = 1e6.*voxelSizeXdouble.*beadsStartX(:,:,1);
        yData(2:2:end) = 1e6.*voxelSizeXdouble.*beadsStartX(:,:,2);
        plot(xData, yData)
        yData(1:2:end) = 1e6.*voxelSizeXdouble.*beadsEndX(:,:,1);
        yData(2:2:end) = 1e6.*voxelSizeXdouble.*beadsEndX(:,:,2);
        plot(xData, yData)
    else
        plot((cropT(1)-1:cropT(2)-1)*voxelSizeTdouble,1e6.*voxelSizeXdouble.* beadsStartX)
        plot((cropT(1)-1:cropT(2)-1)*voxelSizeTdouble,1e6.*voxelSizeXdouble.* beadsEndX)
    end
    plot((cropT(1)-1:cropT(2)-1)*voxelSizeTdouble,1e6.*centresIL(cropT(1):cropT(2)),'k.')
    xlabel('Time(s)')
    ylabel('Position (\mu m)')
    xlim(cropT*voxelSizeTdouble)
    title(['Bead positions (Trap stiffness ' num2str(mean(Stiffness)*1e6) 'pN/um)'])
    legend('LHS','RHS','Centre','Location','southeast')
    subplot(m,1,2)
    histogram(voxelSizeXdouble*1e6*beadWidthsX,'BinWidth',voxelSizeXdouble*0.9e6)
    xlim([0, 5])
    title(['Measured width of 2 \mu m bead (pixel size ' num2str(voxelSizeXdouble*1e6) '\mu m)'])
end
%% Autocorrelation
eta = 0.89e-4;
rBead = 2e-6;
diffCoeff = Kb*T/(6*eta*rBead);
figure(3)
clf
fs = 1/voxelSizeTdouble; % Because interleaved and normalised
centresNorm = zeros(1,size(Centres,2)*nDirections);
centresNorm(1:2:end) = Centres(:,:,1) - mean(Centres(:,:,1));
centresNorm(2:2:end) = Centres(:,:,2) - mean(Centres(:,:,2));
[autocor, lags] = xcorr(centresNorm, ceil(1*fs),'coeff');
autocorNorm = autocor./(2*diffCoeff*abs(lags/fs));
autocorNorm(ceil(length(autocorNorm)/2)) = [];
lags(ceil(length(lags)/2)) = [];
a = fit(log(abs(lags/fs))' ,log(autocorNorm)','poly1');
plot(a,log(abs(lags/fs)),log(autocorNorm))
xlabel('Log Lag (seconds)')
ylabel('Log \lt','Interpreter','Tex')