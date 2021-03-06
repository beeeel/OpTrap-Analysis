%% Process multiple sets of bead data sequentially
% Experiment parameters
mPerPx = 0.07e-6;           % Camera pixel size calibration
laserPowers = 30:5:60;      % Laser power in % for the datasets used
ignoreDirs = {'hela_s3_w_bead_1','hela_s3_w_bead_2',...
    'hela_s3_w_bead_3','hela_s3_w_bead_z_3'}; % Directories to ignore

% Processing parameters
cropTs = {[1 10e5], [1 20e5], [1 20e5], [1 20e5], [1 10e5] };
fitPoly = [0 0 0 0 0 ]; % Fit a polynomial to remove drift. Only do this for calibration sets!
fitPolyOrder = 1;                   % Order of polynomial to be fitted

% Plotting parameters
showStack = false;   % Open the image data in ImageJ
doPlots = true;      % Plot the centres data
compCentres = false; % Show the Imstack with live calculated and offline calculated centres
setLims = false;     % Set axis limits on plots

% Get all the children directories in a struct
dirList = dir;
dirList = dirList([dirList.isdir]);
dirList = dirList(3:end);
for d = 1:length(ignoreDirs)
    dirList = dirList(~strcmp({dirList.name},ignoreDirs{d}));
end
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
    
    % Conditional drift removal only demeans when pOrder = 0
    dims = [1, 3, 2];
    pOrder = fitPolyOrder*fitPoly(fileIdx);
    [~, xCentresM, ~] = func_thermal_rm(1:length(xCentresM), ...
        permute(xCentresM, dims), pOrder, 1, length(xCentresM));
    [~, yCentresM, ~] = func_thermal_rm(1:length(yCentresM), ...
        permute(yCentresM, dims), pOrder, 1, length(yCentresM));
    xCentresM = ipermute(xCentresM, dims);
    yCentresM = ipermute(yCentresM, dims);
    
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
        input('Enter to close figure window')
        close(fh)
    end
    
    % Plot the processed data
    if doPlots
        fh = figure(fileIdx); %#ok<*UNRCH>
        fh.Name = fName;
        clf
        
        % Histogram of xCentres and yCentres in units um
        subplot(3,1,1)
        hold on
        histogram(xCentresM.*1e6,'Normalization','probability')
        histogram(yCentresM.*1e6,'Normalization','probability')
        if setLims
            xlim([-1 1] * 0.05)
        end
        xlabel('Centre position (\mu m)')
        ylabel('Bin probability')
        title(['Histogram of centres, trap stiffness kx = ' num2str(xStiff./1e-6) ' pN/\mu m, ky = ' num2str(yStiff./1e-6) ' pN/\mu m'])
        legend('X','Y')
        
        % Scatterplot of each centre in units um
        subplot(3,1,2)
        plot(xCentresM.*1e6,yCentresM.*1e6,'.')
        if setLims
            xlim([-1 1] * 0.15)
            ylim([-1 1] * 0.15)
        end
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
        if setLims
            ylim([-1 1]*0.15)
        end
        subplot(3,2,6)
        plot(1e-3.*timeVec(cropT(1):cropT(2)), yCentresM.*1e6,'.')
        xlabel('Time (s)')
        ylabel('Y centre (\mu m)')
        title('Y centre position time trace')
        if setLims
            ylim([-1 1]*0.15)
        end
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
%% Try high-pass filtering position
fpass = 0.5; % Pass frequency
fs = length(xCentres)./ (1e-3 * (max(timeVec) - min(timeVec))); % Sampling frequency in Hz
cropT = cropTs{fileIdx};
[aVarX, Tau, xCentresHP] = func_bead_hp_allan_var(xCentres, timeVec, fpass, cropT, 1, fName, true);
[aVarY, ~, yCentresHP] = func_bead_hp_allan_var(yCentres, timeVec, fpass, cropT, 1, fName, true);
%% Look at mean-square displacement (for cell-bead expts)
idx = [1 2e6];
msd = msdanalyzer(1, 'um', 'ms');
msd = msd.addAll({[timeVec(idx(1):idx(2))' 1e6.*xCentresM(idx(1):idx(2))']});
tic
msd = msd.computeMSD;
toc
fh = figure;
fh.Name = num2str(diff(idx)+1);
msd.plotMSD
%%
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
    str = 'fitPoly = [';
    str = [str repmat('0 ', 1, length(struct))];
    str = [str '];'];
    disp(str);
    error('CropTs is the wrong size, copy the above lines into the script');
end
end

function changeIm(sld, ims, ax, imCentres, xyCentres) 
fr = round(sld.Value);
cla(ax);
imagesc(ax, ims{1}{fr,1});
plot(ax, imCentres(1, fr), imCentres(2, fr), 'kx');
plot(ax, xyCentres(1, fr), xyCentres(2, fr), 'k.');
axis(ax, 'image');
legend(ax, 'MATLAB calculated centre','Live calculated centre')
diffCentres = imCentres - xyCentres(:, 1:length(imCentres));
title(ax,['Frame ' num2str(fr) ' x_{diff} = ' num2str(diffCentres(1,fr)) ...
    ' y_{diff} = ' num2str(diffCentres(2,fr))]);
end
