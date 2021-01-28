%% Process multiple sets of bead data sequentially
% Experiment parameters
mPerPx = 0.07e-6;           % Camera pixel size calibration
laserPowers = 30:5:60;      % Laser power in % for the datasets used
ignoreDirs = {}; % Directories to ignore

% Processing parameters
cropTs = {[1 6e4], [3e4 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4], [1 6e4] };
fitPoly = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ]; % Fit a polynomial to remove drift
fitPolyOrder = 1;       % Order of polynomial to be fitted
calcStiff = 0;          % Calculate trap stiffness from position variance
fpass = 0;              % Pass frequency
msdOffset = 1;          % Offset from start when taking data to calculate mean-square displacements
freshStart = false;

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

% Preallocate 
if freshStart 
    out = {};
end

for fileIdx = 18:26%length(dirList)
    %% Load and pre-process
    % Load all the data and the metadata
    if freshStart || ~isstruct(out{fileIdx})
        data = struct([]);
        data.forceRun = freshStart;
        
        % Set names and load data
        data(1).dirPath = [dirList(fileIdx).folder '/' dirList(fileIdx).name];
        data.fName = dirList(fileIdx).name;
        data = bead_loadData(data);
        
        % Apply calibration and crop time
        data.opts.cropT = cropTs{fileIdx};
        data.opts.pOrder = fitPolyOrder*fitPoly(fileIdx);
        data.mPerPx = mPerPx;
    else
        data = out{fileIdx};
        data.forceRun = false;
    end
        
    data = bead_preProcessCentres(data);
    data = bead_fft_scaled(data, doPlots);
    
    if saveFigs
        saveas(gcf, [data.fName '_fft.png'])
    end
    %% Process data
    % Calculate the stiffnesses and put into data
    if calcStiff
        xStiff = calcStiffness(data.pro.xCentresM);
        yStiff = calcStiffness(data.pro.yCentresM);
        data.pro.stiffXY(:, fileIdx) = [xStiff, yStiff];
    end
    
    % Compare MATLAB calculated centres with live (Java) calculated centres
    if compCentres
        bead_plotCompareCentres(data)
    end
    
    % Plot the processed data
    if doPlots
        fh = bead_plotProData(data, setLims); %#ok<*UNRCH>
        if saveFigs
            saveas(fh, [data.fName '_raw.png'])
        end
    end
    
    % Open the Imstack file using ImageJ (kinda redundant since I have
    % compCentres)
    if showStack
        disp('MATLAB is locked until you close imagej')
        system(['imagej ' dirPath '/images_and_metadata/images_and_metadata_MMStack_Default.ome.tif']);
    end
    
    % High-pass filter and calculate Allan variance if pass frequency is
    % positive
    if fpass > 0
        data = bead_hp_allan_var(data, 'xCentresPx', ...
            fpass, 5e3, doPlots);
%         data = bead_hp_allan_var(data, 'yCentresPx', ...
%             fpass, 1, doPlots);
    end
    
    % Look at mean-square displacement (for cell-bead expts)
    if msdOffset
        data = bead_normMSD_polyfit(data, 'xCentresPx', msdOffset, 6e4);
        if saveFigs
            saveas(gcf, [data.fName '_MSD.png'])
        end
%         data = bead_normMSD(data, 'yCentresM', msdOffset);
    end
    
    out{fileIdx} = data; %#ok<SAGROW>
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
