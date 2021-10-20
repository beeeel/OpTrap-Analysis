%% Process multiple sets of bead data sequentially
%
% V1: Perform iterative short time processing as per Manlio's suggestion:
%       Perform initial processing as normal, then find time of long time
%       uptick/low frequency intercept. Use this to apply highpass
%       filtering before recalculating MSD and FTs.
%
% This file accumulates: Control data

%% Data parameters
close all

% Data storage parameters
masterDir = '../';             % Where all the data is
saveDir = 'accumulator_2021_07_21-latB/';        % Where to save processed data

% Which day's data to load
dayDirs = {'2021_07_27' '2021_07_26' '2021_07_23'};% '2021_07_13' '2021_07_16' '2021_05_17'};

% This cell contains 1 cell per dayDir above. Within that needs to be
% indexes for that day's dirList, to choose the correct datasets.
setIdxs = {
    ... % Sets from 2021_07_27
        {1:9 [10:13 15:18]} ...
    ... % Sets from 2021_07_26
        {[5:10 12 14]} ...
    ... % Sets from 2021_07_23
        {1:10 11:20 [22:27 29:31]} ...
    ... % Sets from 2021_07_13
        {6:17} ...
    ... % Sets from 2021_07_16
        {8:14} ... {1:2 3:5 6:12} ...
    ... % Sets from 2021_05_17
        {5:9 10:14} };

% Ain't nobody got time (be)for this!
% Later I will record time in a filename when recording data
thymes = { 
    ... % Times from 2021_07_27
        {[-40 -20 -10 7 17 32 42 56 82] [-26 -15 1 12 37 48 61 67]} ...
    ... % Times from 2021_07_26
        {[-14 -7 7 14 27 32 46 61]} ...
    ... % Times from 2021_07_23
        {[-30 -14 2 8 15 22 28 35 43 60] [-18 -12 -6 1 7 16 39 52 69 96] [-18 -6 0 13 22 26 46 54 60]} ...
    ... % Times from 2021_07_13
        {[49 63 66 71 77 83 88 94 101 107 113 120]-60} ...
    ... % Times from 2021_07_16
        {[2 8 11 15 18 31 48]-8} ... {[0 6] [0 4 10] [0 2 8 11 15 18 31]} ...
    ... % Times from 2021_05_17
        {([80 103 142 165 203]-100) ([91 107 156 172 211]-100)}};
%% Experiment parameters
mPerPx = 0.065e-6;           % Camera pixel size calibration
ignoreDirs = {'focal_sweep_with_bead'}; % Directories to ignore (ones without data)

%% Processing parameters
cropTs = {[]};
fitPoly = 1; % Fit a polynomial to remove drift.
fitPolyOrder = 1;      % Order of polynomial to be fitted

angleCorrection = true; % Transform to (r, rθ) co-ordinates
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
saveAccu = true;        % Save data to file
dataSuff = '_simple';       % Suffix for filename when saving/loading
accuFile = 'accumulated_3days';

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
                data = bead_loadData(data, angleCorrection);
                
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

            data.opts.angleCorrection = angleCorrection;

            % Store filenames
            accumulated{dayIdx}{1,cellIdx}(fIdx).fName = data.fName;

            % Apply scale and do polyfit
            data = bead_preProcessCentres(data);
            
            % Record where the data is coming from
            accumulated{dayIdx}{1,cellIdx}(fIdx).dirPath = data.dirPath;
            accumulated{dayIdx}{1,cellIdx}(fIdx).fName = data.fName;

            % Update times array to show time halfway through acquisition
            accumulated{dayIdx}{2,cellIdx}(fIdx) = accumulated{dayIdx}{2,cellIdx}(fIdx) ...
                + data.raw.timeVecMs(end/2)/60e3;
            %% Process data

            % Look at mean-square displacement (for cell-bead expts)
            if ~isempty(msdOffset)
                data = bead_normMSD_polyfit(data, msdDim, msdOffset, msdNumT, false, msdUseRaw, msdCentresRow, msdDoNorm);
                accumulated{dayIdx}{1,cellIdx}(fIdx).msd = data.pro.amsdObj;
            end
            
            % High-pass filter and calculate Allan variance if pass frequency is positive
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
            
            
            
            % Save if requested
            if saveData
                save(dataFile, 'data')
            end
            
        end
    end
end

if saveAccu
    save(accuFile, 'accumulated', 'dayDirs')
end