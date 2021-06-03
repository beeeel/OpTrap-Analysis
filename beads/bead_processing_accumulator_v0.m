%% Process multiple sets of bead data sequentially
close all

% Data storage parameters
masterDir = '~/Data/phd/OpTrap/';           % Where all the data is
dayDirs = {'2021_05_17/' '2021_05_13/' '2021_04_22/'};    % Which day's data to load
setIdxs = {{1:3 5:7 10:12} {3:6 7:9 10:13 16:19 20:22 23:25} {1:5 6:8 }};
% Experiment parameters
mPerPx = 0.07e-6;           % Camera pixel size calibration
laserPowers = 0;      % Laser power in % for the datasets used
ignoreDirs = {'focal_sweep_with_bead'}; % Directories to ignore (ones without data)

% Processing parameters
cropTs = {};
fitPoly = 1; % Fit a polynomial to remove drift.

fitPolyOrder = 18;      % Order of polynomial to be fitted
calcStiff = 1;          % Calculate trap stiffness from position variance
fpass = 1;              % Pass frequency
cropTHPval = 0;         % Frames to crop after HP filter
msdOffset = [1];        % Offset from start when taking data to calculate mean-square displacements

msdDim = 'all';         % Direction to calculate MSD in - 'x', 'y', or 'all'
msdCentresRow = 1;      % Row of centres array to use for MSDs (empty for default)
msdNumT = [];           % Number of time points to use for MSDs (empty for all)
msdUseRaw = false;      % Use raw or processed data for MSDs (empty for default)
msdDoNorm = true;       % Normalize MSDs by position variance (empty for default)
doFFT = false;           % Calculate FFT and maybe plot

% Data file parameters
forceRun = true;       % Try to take data from file and reuse as much as possible
saveData = false;        % Save data to file
dataSuff = '_simple';       % Suffix for filename when saving/loading

% Prepare for trouble
out = struct();
% And make it double!
out(1).stiff = double(nan);

% For each day chosen
for dayIdx = 1:length(dayDirs)
    % Go to the correct directory
    dayDir = [masterDir dayDirs{dayIdx}];
    cd(dayDir);
    
    % Get all the children directories in a struct
    dirList = dir;
    dirList = dirList([dirList.isdir]);
    dirList = dirList(3:end);
    for d = 1:length(ignoreDirs)
        dirList = dirList(~strcmp({dirList.name},ignoreDirs{d}));
    end
    
    for fileIdx = setIdxs{dayIdx}
        
        %% Load and pre-process
        % Either create a new struct or load one named dataFile
        dataFile = [dirList(fileIdx).name '_processed' dataSuff '.mat'];
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
            tmp = load(dataFile, 'data');
            data = tmp.data;
            clear tmp
            data.opts.forceRun = forceRun;
        end
        
        % Apply scale and do polyfit
        data = bead_preProcessCentres(data);
        %% Process data
        % Calculate the stiffnesses and put into data
        if calcStiff
            stiffIdx = 1;
            out(fileIdx).suffix = data.raw.suffixes{stiffIdx};
            
            xStiff = calcStiffness(data.pro.xCentresM);
            yStiff = calcStiffness(data.pro.yCentresM);
            data.pro.stiffXYpro = [xStiff, yStiff];
            out(fileIdx).stiff = data.pro.stiffXYpro(stiffIdx,:);
            
            xStiff = calcStiffness(data.raw.xCentresPx, data.mPerPx);
            yStiff = calcStiffness(data.raw.yCentresPx, data.mPerPx);
            data.pro.stiffXYraw = [xStiff, yStiff];
            out(fileIdx).stiffraw = data.pro.stiffXYraw(stiffIdx,:);
        end
        
        %% High-pass filter and calculate Allan variance if pass frequency is
        % positive
        if fpass > 0
            data.opts.fpass = fpass;
            data.opts.cropTHPval = cropTHPval;
            % Need to adapt to run on both dims and use subplots/etc
            data = bead_hp_allan_var(data, 'xCentresPx', ...
                fpass, cropTHPval, false);
            data = bead_hp_allan_var(data, 'yCentresPx', ...
                fpass, cropTHPval, false);
            
            xStiff = calcStiffness(data.pro.xCentresHP);
            yStiff = calcStiffness(data.pro.yCentresHP);
            data.pro.stiffXYHP = [xStiff, yStiff];
            out(fileIdx).stiffHP = data.pro.stiffXYHP(2,:);
            
        end
        
        % Look at mean-square displacement (for cell-bead expts)
        if ~isempty(msdOffset)
            data = bead_normMSD_polyfit(data, msdDim, msdOffset, msdNumT, false, msdUseRaw, msdCentresRow, msdDoNorm);
        end
        
        % Save if requested
        if saveData
            save(dataFile, 'data')
        end
        
    end
end

