%% Process multiple sets of bead data sequentially
% Experiment parameters
mPerPx = 0.07e-6;
laserPowers = 30:5:60;
cropTs = {[1 5e5], [1 5e5], [1 5e5], [1e5 3.2e5], [1.4e5 4.4e5], [1.7e5 5e5], [1 0.9e5], [1.6e5 2.8e5] };

% Processing parameters
fitPoly = [0, 0, 0, 0, 0, 0, 0, 1]; % Fit a polynomial to remove drift. Only do this for calibration sets!
fitPolyOrder = 1;                   % Order of polynomial to be fitted

% Plotting parameters
showStack = false;   % Open the image data in ImageJ
doPlots = true;      % Plot the centres data
compCentres = false; % Show the Imstack with live calculated and offline calculated centres

% Get all the children directories in a struct
dirList = dir;
dirList = dirList([dirList.isdir]);
dirList = dirList(3:end);
% Check cropTs has been made with the right number of elements. Throws an
% error and gives an empty list if not.
checkCropTs(cropTs, dirList);

% Preallocate 
stiffXY = zeros(2,length(laserPowers));

for fileIdx = 8%1:length(dirList)
    % Load all the data and the metadata
    dirPath = [dirList(fileIdx).folder '/' dirList(fileIdx).name];
    fName = dirList(fileIdx).name;
    xCentres = byteStreamToDouble([dirPath '/XCentres.dat']);
    yCentres = byteStreamToDouble([dirPath '/YCentres.dat']);
    timeVec  = byteStreamToDouble([dirPath '/Times.dat']);
    Imstack  = bfopen([dirPath '/images_and_metadata/images_and_metadata_MMStack_Default.ome.tif']);
    metadata = fileread([dirPath '/images_and_metadata/images_and_metadata_MMStack_Default_metadata.txt']);
    metadata = jsondecode(metadata);
    
    % Apply calibration and calculate stiffnesses
    cropT = cropTs{fileIdx};
    xCentresM = xCentres(cropT(1):cropT(2)) .* mPerPx;
    yCentresM = yCentres(cropT(1):cropT(2)) .* mPerPx;
    
    % Conditional drift removal (only use for calibration sets!)
    if fitPoly(fileIdx)
        dims = [1, 3, 2];
        [~, xCentresM, ~] = func_thermal_rm(1:length(xCentresM), ...
            permute(xCentresM, dims), fitPolyOrder, 1, length(xCentresM));
        [~, yCentresM, ~] = func_thermal_rm(1:length(yCentresM), ...
            permute(yCentresM, dims), fitPolyOrder, 1, length(yCentresM));
        xCentresM = ipermute(xCentresM, dims);
        yCentresM = ipermute(yCentresM, dims);
    end
    
    % Calculate the stiffnesses and put into output array
    xStiff = calcStiffness(xCentresM);
    yStiff = calcStiffness(yCentresM);
    stiffXY(:, fileIdx) = [xStiff, yStiff];
    
    % Compare MATLAB calculated centres with live (Java) calculated centres
    if compCentres
        % Get two centres arrays
        imCentres = imCentreOfMass(cat(3,Imstack{1}{:,1}));
        xyCentres = [xCentres(1:1e4:end); yCentres(1:1e4:end)];
        % Create a UIFigure with axes holding plots and plot the first
        % image
        fh = uifigure('Name','Image with both calculated centres',...
            'Position',[680 160 980 800]);
        ax = axes(fh);
        hold(ax,'on')
        changeIm(struct('Value',1), Imstack, ax, imCentres, xyCentres);
        % Create a slider for user control of which image is showing
        sld = uislider(fh, 'Position', [100, 50, 600, 40], ...
            'ValueChangedFcn', @(sld, event) changeIm(sld, Imstack, ax, imCentres, xyCentres),...
            'Limits', [1 length(Imstack{1})], 'MinorTicks', 1:length(Imstack{1}));
    end
    
    % Plot the processed data
    if doPlots
        fh = figure(fileIdx); %#ok<*UNRCH>
        fh.Name = fName;
        clf
        
        % Histogram of xCentres and yCentres in units um
        subplot(3,1,1)
        hold on
        histogram(xCentresM.*1e6)
        histogram(yCentresM.*1e6)
        xlabel('Centre position (\mu m)')
        ylabel('Bin count')
        title(['Histogram of centres, trap stiffness kx = ' num2str(xStiff./1e-6) ' pN/\mu m, ky = ' num2str(yStiff./1e-6) ' pN/\mu m'])
        legend('X','Y')
        
        % Scatterplot of each centre in units um
        subplot(3,1,2)
        plot(xCentresM.*1e6,yCentresM.*1e6,'.')
        title(['Scatterplot of centres, ' num2str(diff(cropT)+1) ' frames'])
        xlabel('X Centre position (\mu m)')
        ylabel('Y Centre position (\mu m)')
        axis equal
        
        % Time traces of X and Y in units um 
        subplot(3,2,5)
        plot(1e-3.*timeVec(cropT(1):cropT(2)), xCentresM.*1e6,'.')
        xlabel('Time (s)')
        ylabel('X centre (\mu m)')
        title('X centre position time trace')
        subplot(3,2,6)
        plot(1e-3.*timeVec(cropT(1):cropT(2)), yCentresM.*1e6,'.')
        xlabel('Time (s)')
        ylabel('Y centre (\mu m)')
        title('Y centre position time trace')
        drawnow
    end
    
    % Open the Imstack file using ImageJ (kinda redundant since I have
    % compCentres)
    if showStack
        disp('MATLAB is locked until you close imagej')
        system(['imagej ' dirPath '/images_and_metadata/images_and_metadata_MMStack_Default.ome.tif']);
    end
end
%% Plot stiffnesses of all sets
figure(fileIdx+1)
clf
subplot(2,1,1)
hold on
plot(laserPowers,1e6*stiffXY(1,1:length(laserPowers)))
plot(laserPowers,1e6*stiffXY(2,1:length(laserPowers)))
title('Stiffness varying with laser power setting')
legend('X stiffness', 'Y stiffness')
xlabel('Laser power setting (%)')
ylabel('Trap stiffness (pN/\mu m)')
subplot(2,1,2)
bar(laserPowers, stiffXY(1,1:length(laserPowers))./stiffXY(2,1:length(laserPowers)))
title('Ratio $\frac{k_x}{k_y}$','Interpreter','latex','Fontsize',20)
xlabel('Laser power setting (%)')
ylabel('Ratio')
%% Show calculated centres with image
frames = 1:15;

imCentres = imCentreOfMass(cat(3,Imstack{1}{frames,1}));
figure(4)
fh.Name = fName;
clf
hold on
plot(0.07.*xCentres(1:1e4:max(frames)*1e4),0.07.*yCentres(1:1e4:max(frames)*1e4),'k.')
plot(0.07.*imCentres(1,:), 0.07.*imCentres(2,:), 'rx')
axis equal

figure(5)
fh.Name = fName;
for frame = 1%:16
clf
hold on
imagesc(Imstack{1}{frame,1})
plot(imCentres(1,frame), imCentres(2,frame),'x')
plot(xCentres((frame-1).*1e4+1)+0.5,yCentres((frame-1).*1e4+1)+0.5,'.')
legend('MATLAB measured centre','ImageJ measured centre')
title(num2str(frame))
pause(0.5)
end

function checkCropTs(cell, struct)
if length(cell) ~= length(struct)
    str = 'cropTs = {';
    str = [str repmat('[1 5e5], ', 1, length(struct))];
    str = [str(1:end-2) ' };'];
    disp(str);
    error('CropTs is empty, copy the above line into the script');
end
end

function changeIm(sld, ims, ax, imCentres, xyCentres) 
fr = round(sld.Value);
cla(ax);
imagesc(ax, ims{1}{fr,1});
plot(ax, imCentres(1, fr), imCentres(2, fr), 'kx');
plot(ax, xyCentres(1, fr), xyCentres(2, fr), 'k.');
legend(ax, 'MATLAB calculated centre','Live calculated centre')
diffCentres = imCentres - xyCentres(:, 1:length(imCentres));
title(ax,['Frame ' num2str(fr) ' x_{diff} = ' num2str(diffCentres(1,fr)) ...
    ' y_{diff} = ' num2str(diffCentres(2,fr))]);
end