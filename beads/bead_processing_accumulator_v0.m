%% Process multiple sets of bead data sequentially
%
% This file accumulates: Control data

%% Data parameters
close all

% Data storage parameters
masterDir = '~/Data/phd/OpTrap/';           % Where all the data is
saveDir = '2021_07_19';        % Where to save processed data

% Which day's data to load
dayDirs = {'2021_07_19'};

% This cell contains 1 cell per dayDir above. Within that needs to be
% indexes for that day's dirList, to choose the correct datasets.
setIdxs = {
    ... % Sets going up in 5μm steps
        ...{[5 6 8 10 11 13 14]}... 16 17 19 20 22 23 25 26]} ...
    ... % Sets going back down in 10μm steps
        {[7 9 12 15 18 21 24 27]} ...
    };

% Ain't nobody got time (be)for this!
% Later I will record time in a filename when recording data
thymes = { 
    ... % Going up
        ...{[0:5:30]}...70]} ...
    ... % Going down
        {[0:10:70]}...
    }; % Actually for this data these are depths

%% Experiment parameters
mPerPx = 0.065e-6;           % Camera pixel size calibration
ignoreDirs = {'focal_sweep_with_bead'}; % Directories to ignore (ones without data)

%% Processing parameters
cropTs = {[]};
fitPoly = 1; % Fit a polynomial to remove drift.

fitPolyOrder = 1;      % Order of polynomial to be fitted
calcStiff = 1;          % Calculate trap stiffness from position variance
stiffIdx = 1;           % Row of centres to store stiffness from
fpass = 0;              % Pass frequency
cropTHPval = 0;         % Frames to crop after HP filter
msdOffset = [1];        % Offset from start when taking data to calculate mean-square displacements

msdDim = 'all';         % Direction to calculate MSD in - 'x', 'y', or 'all'
msdCentresRow = 1;      % Row of centres array to use for MSDs (empty for default)
msdNumT = [];           % Number of time points to use for MSDs (empty for all)
msdUseRaw = true;      % Use raw or processed data for MSDs (empty for default)
msdDoNorm = true;       % Normalize MSDs by position variance (empty for default)
doFFT = false;           % Calculate FFT and maybe plot

% Data file parameters
forceRun = true;       % Try to take data from file and reuse as much as possible
saveData = false;        % Save data to file
dataSuff = '_simple';       % Suffix for filename when saving/loading

% ARE YOU READY??
accumulated = cell(size(dayDirs));
if ~exist([masterDir saveDir], 'dir')
    mkdir([masterDir saveDir])
end
cd([masterDir saveDir])
fprintf('Working in %s\n',pwd)
fileCount = 1;
%%
% For each day chosen
for dayIdx = 1:length(dayDirs)
    % Get the correct directory
    dayDir = [masterDir dayDirs{dayIdx}];
    
    % Get all the children directories in a struct
    dirList = dir(dayDir);
    dirList = dirList([dirList.isdir]);
    dirList = dirList(3:end);
    for d = 1:length(ignoreDirs)
        dirList = dirList(~strcmp({dirList.name},ignoreDirs{d}));
    end
    
    % Prepare for trouble
    accumulated{dayIdx} = {};
    
    for cellIdx = 1:length(setIdxs{dayIdx})
        % Prepare for trouble
        accumulated{dayIdx}{1,cellIdx} = struct();
        % And make it double
        accumulated{dayIdx}{1,cellIdx}.stiff = double(nan);
        % Get the times array for this set
        accumulated{dayIdx}{2,cellIdx} = thymes{dayIdx}{cellIdx};
        
        for fIdx = 1:length(setIdxs{dayIdx}{cellIdx})
            fileIdx = setIdxs{dayIdx}{cellIdx}(fIdx);
            
            %% Load and pre-process
            % Either create a new struct or load one named dataFile
            dataFile = [dirList(fileIdx).folder '/' dirList(fileIdx).name '_processed' dataSuff '.mat'];
            if forceRun || ~exist(dataFile, 'file')
                data = struct([]);
                data(1).opts.forceRun = forceRun;
                
                % Set names and load data
                data.dirPath = [dirList(fileIdx).folder '/' dirList(fileIdx).name];
                data.fName = dirList(fileIdx).name;
                data = bead_loadData(data, false);
                
                % Apply calibration and crop time
                data.opts.cropT = cropTs{1};
                data.opts.pOrder = fitPolyOrder*fitPoly;
                data.mPerPx = mPerPx;
            else
                tmp = load(dataFile, 'data');
                data = tmp.data;
                clear tmp
                data.opts.forceRun = forceRun;
            end
            
            if dayIdx == 1 && fIdx == 1 
                data.raw.xCentresPx = data.raw.xCentresPx(2,:);
                data.raw.yCentresPx = data.raw.yCentresPx(2,:);
                disp('Replaced centres once')
%             elseif dIdx == 1
%                 data.raw.xCentresPx = data.raw.xCentresPx(3,:);
%                 data.raw.yCentresPx = data.raw.yCentresPx(3,:);
            end
            
            % Store filenames
            accumulated{dayIdx}{1,cellIdx}(fIdx).fName = data.fName;
            
            % Apply scale and do polyfit
            data = bead_preProcessCentres(data);
            %% Process data
            % Calculate the stiffnesses and put into data
            if calcStiff
                accumulated{dayIdx}{1,cellIdx}(fIdx).suffix = data.raw.suffixes{stiffIdx};
                
                xStiff = calcStiffness(data.pro.xCentresM);
                yStiff = calcStiffness(data.pro.yCentresM);
                data.pro.stiffXYpro = [xStiff, yStiff];
                accumulated{dayIdx}{1,cellIdx}(fIdx).stiff = data.pro.stiffXYpro(stiffIdx,:);
                
                xStiff = calcStiffness(data.raw.xCentresPx, data.mPerPx);
                yStiff = calcStiffness(data.raw.yCentresPx, data.mPerPx);
                data.pro.stiffXYraw = [xStiff, yStiff];
                accumulated{dayIdx}{1,cellIdx}(fIdx).stiffraw = data.pro.stiffXYraw(stiffIdx,:);
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
                accumulated{dayIdx}{1,cellIdx}(fIdx).stiffHP = data.pro.stiffXYHP;
                
            end
            
            % Look at mean-square displacement (for cell-bead expts)
            if ~isempty(msdOffset)
                data = bead_normMSD_polyfit(data, msdDim, msdOffset, msdNumT, false, msdUseRaw, msdCentresRow, msdDoNorm);
                accumulated{dayIdx}{1,cellIdx}(fIdx).msd = data.pro.amsdObj;
            end
            
            % Save if requested
            if saveData
                save(dataFile, 'data')
            end
            
        end
    end
end
