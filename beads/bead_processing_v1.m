%% Process multiple sets of bead data sequentially
% Experiment parameters
mPerPx = 0.07e-6;           % Camera pixel size calibration
laserPowers = 0;      % Laser power in % for the datasets used
ignoreDirs = {}; % Directories to ignore

% Processing parameters
cropTs = {[1 6e4], [3e4 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4] };
fitPoly = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]; % Fit a polynomial to remove drift
fitPolyOrder = 1;       % Order of polynomial to be fitted
calcStiff = 0;          % Calculate trap stiffness from position variance
fpass = 0;              % Pass frequency
msdOffset = 1;          % Offset from start when taking data to calculate mean-square displacements
msdDim = 'all';         % Direction to calculate MSD in - 'x', 'y', or 'all'
doFFT = true;

% Data file parameters
forceRun = false;       % Try to take data from file and reuse as much as possible
saveData = true;        % Save data to file

% Plotting parameters
saveFigs = false;
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

for fileIdx = 27:35%length(dirList)
    %% Load and pre-process
    % Either create a new struct or load one named dataFile
    dataFile = [dirList(fileIdx).name '_processed.mat'];
    if forceRun || ~exist(dataFile, 'file')
        data = struct([]);
        data(1).opts.forceRun = forceRun;
        
        % Set names and load data
        data.dirPath = [dirList(fileIdx).folder '/' dirList(fileIdx).name];
        data.fName = dirList(fileIdx).name;
        data = bead_loadData(data);
        
        % Apply calibration and crop time
        data.opts.cropT = cropTs{fileIdx};
        data.opts.pOrder = fitPolyOrder*fitPoly(fileIdx);
        data.mPerPx = mPerPx;
    else
        data = load(dataFile, 'data');
        data.opts.forceRun = forceRun;
    end
    
    % Apply scale and do polyfit
    data = bead_preProcessCentres(data);
    %% Process data
    % Calculate the stiffnesses and put into data
    if calcStiff
        xStiff = calcStiffness(data.pro.xCentresM);
        yStiff = calcStiffness(data.pro.yCentresM);
        data.pro.stiffXY(:, fileIdx) = [xStiff, yStiff];
    end
    
    % Plot the processed data
    if doPlots
        fh = bead_plotProData(data, setLims); %#ok<*UNRCH>
        if saveFigs
            saveas(fh, [data.fName '_pro.png'])
        end
    end
    
    % Calculate frequency spectrum in physical units
    if doFFT
        data = bead_fft_scaled(data, doPlots);
        if saveFigs
            saveas(gcf, [data.fName '_fft.png'])
        end
    end
    
    % Compare MATLAB calculated centres with live (Java) calculated centres
    if compCentres
        bead_plotCompareCentres(data)
    end
    
    % Open the Imstack file using ImageJ
    if showStack
        % System can be called with an & in there to run in background
        command = ['imagej ' data.dirPath '/images_and_metadata/images_and_metadata_MMStack_Default.ome.tif '];
        if isfield(data, 'ImstackFullFoV')
            command = [command data.dirPath '/full_images_and_metadata/full_images_and_metadata_MMStack_Default.ome.tif '];
        end
        system([command '&']);
    end
    
    % High-pass filter and calculate Allan variance if pass frequency is
    % positive
    if fpass > 0
        %% Need to adapt to run on both dims and use subplots/etc
        data = bead_hp_allan_var(data, 'xCentresPx', ...
            fpass, 5e3, doPlots);
%         data = bead_hp_allan_var(data, 'yCentresPx', ...
%             fpass, 1, doPlots);
    end
    
    % Look at mean-square displacement (for cell-bead expts)
    if msdOffset
        data = bead_normMSD_polyfit(data, msdDim, msdOffset, 6e4);
        if saveFigs
            saveas(gcf, [data.fName '_MSD.png'])
        end
    end
    
    % Save if requested
    if saveData
        save(dataFile, 'data')
    end

%%

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
