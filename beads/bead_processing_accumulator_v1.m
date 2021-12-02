%% Process multiple sets of bead data sequentially
%
% V1: Perform iterative short time processing as per Manlio's suggestion:
%       Perform initial processing as normal, then find time of long time
%       uptick/low frequency intercept. Use this to apply highpass
%       filtering before recalculating MSD and FTs.
%
% This file accumulates: Latrunculin B data
%
%% TO DO:
%
% Better warnings when no corners: Filename and which corner
% Quieter warnings when plotting etc - there's always negative Y
% Finding α, D from MSDs
% Store α,D from fitting to tRanges - for high static error, is D constant?
% Storing MSD(1./fpass)
% Better MSD field handling - need to use centresHP when possible and
%  centresM otherwise
% data.opts.nSkip - use this in any function which handles MSDs

%% Data parameters
close all

% Data storage parameters
masterDir = '../';             % Where all the data is
saveDir = 'accumulator_2021_07_21-latB/';        % Where to save processed data

% Which day's data to load
dayDirs = {'2021_07_27' '2021_07_26' '2021_07_23' '2021_07_13' '2021_07_16' '2021_05_17'};

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
    
wRanges = get_wRanges_struct(dayDirs);
tRanges = get_tRanges_struct(dayDirs);
tRangesLow = get_low_tRanges_struct(dayDirs);
%% Experiment parameters
mPerPx = 0.065e-6;           % Camera pixel size calibration
ignoreDirs = {'focal_sweep_with_bead'}; % Directories to ignore (ones without data)

%% Processing parameters
cropTs = {[]};
fitPoly = 1; % Fit a polynomial to remove drift.
fitPolyOrder = 0;      % Order of polynomial to be fitted

angleCorrection = true; % Transform to (r, rθ) co-ordinates
calcStiff = 1;          % Calculate trap stiffness from position variance
stiffIdx = 1;           % Row of centres to store stiffness from
fpass = 0;              % Pass frequency
cropTHPval = 0;         % Frames to crop after HP filter
msdOffset = [1];        % Offset from start when taking data to calculate mean-square displacements
loadPics = false;       % Load images from TIFs. May be overriden by angleCorrection

msdDim = 'all';         % Direction to calculate MSD in - 'x', 'y', or 'all'
centresRow = 1;      % Row of centres array to use for MSDs (empty for default)
msdNumT = [];           % Number of time points to use for MSDs (empty for all)
msdUseRaw = false;      % Use raw or processed data for MSDs (empty for default)
msdDoNorm = false;       % Normalize MSDs by position variance (empty for default)
doFFT = false;           % Calculate FFT and maybe plot

% Data file parameters
forceRun = false;       % Try to take data from file and reuse as much as possible
saveDataPro = false;    % Save processed data to file (probably only makes things slower for now)
saveDataRaw = false;     % Save raw data to file (should speed up loading)
saveAccu = true;        % Save data to file
dataSuff = '_simple';       % Suffix for filename when saving/loading
accuFile = 'accumulated_all_c';

% ARE YOU READY??
accumulated = cell(size(dayDirs));
if ~exist([masterDir saveDir], 'dir')
    mkdir([masterDir saveDir])
end
cd([masterDir saveDir])
fprintf('Working in %s\n',pwd)
fileCount = 1;
%%
fh = figure(69);
fh2 = figure(99);
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
                data = bead_loadData(data, any([loadPics, angleCorrection]));
                fprintf('\n\t\t Loaded:\t\t%s\n\n',data.fName)
                
                % Save short data if requested
                if saveDataRaw
                    save(dataFile, 'data')
                end
            else
                tmp = load(dataFile, 'data');
                data = tmp.data;
                clear tmp
                data.opts.forceRun = forceRun;
            end

            % Record options and calibration
            data.opts.cropT = cropTs{1};
            data.opts.pOrder = fitPolyOrder*fitPoly;
            data.mPerPx = mPerPx;
            data.opts.angleCorrection = angleCorrection;
            data.opts.centresRow = centresRow;

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
                data = bead_normMSD_polyfit(data, msdDim, msdOffset, msdNumT, false, msdUseRaw, centresRow, false);
                accumulated{dayIdx}{1,cellIdx}(fIdx).msd = data.pro.amsdObj;
            end
            
            % Get wRange from struct
            wR = Range_getter(wRanges, dayDirs{dayIdx}, cellIdx, fIdx);
            
            % Do the fourier transform to find intercept frequency
            [FT, oC] = msd_fourier_transformator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), ...
                'wRange',wR, 'fh', fh);
            oC = cat(1,oC{:});
            
            % Get tRange from struct
            tR = Range_getter(tRanges, dayDirs{dayIdx}, cellIdx, fIdx);
            % Also fit to the MSD to find corner time equivalent frequency
            tC = msd_cornerator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), tR);
            
            if any([ isnan( tC ) & isnan( oC' ) , isempty(tC)])
                tr = {'tangential','radial'};
                if sum( isnan( tC ) & isnan( oC' ) ) == 1
                    % Could give name of file in warning
                    warning('No corner in %s direction, highpass only using', tr{ sum( isnan( tC ) & isnan( oC' ) ) })
                else
                    warning('No corner in either direction, skipping highpass.')
                end
                % If you've come here looking for answers, I'm sorry.
                warning('The last warning was not well defined because stuff changes and it''s all a mess')
                
                % Could give better warnings:
%             else
%                 if any( isnan( tC ) )
%                     warning('Got NaN back from msd_cornerator. Highpass only using corner frequency')
%                 end
%                 if any( isnan( oC ) )
%                     warning('Got NaN back from msd_fourier_transformator. Highpass only using corner time')
%                 end
            end
            if isempty(tC)
                tC = nan(1,2);
            end
            
            % Do something with the corners
            accumulated{dayIdx}{1,cellIdx}(fIdx).corners = [1./tC; oC'];
            
            % Use mean of 1/corner time and intercept frequency. There's
            % definitely scope for mistakes to occur if you want to be
            % getting smarter, but that will have to wait. 
            fpass = mean([1./tC; oC'],'omitnan');
            
            % This has assumed you've only given one intercept frequency
            % range and one corner range
            
            % High-pass filter and calculate Allan variance if pass frequency is positive
            if any(fpass > 0)
                data.opts.fpass = fpass;
                data.opts.cropTHPval = cropTHPval;
                % Need to adapt to run on both dims and use subplots/etc
                if fpass(1) > 0
                    data = bead_hp_allan_var(data, 'xCentresM', ...
                        fpass(1), cropTHPval, false);
                    xStiff = calcStiffness(data.pro.xCentresHP);
                else
                    xStiff = calcStiffness(data.pro.xCentresM(centresRow,:));
                end
                
                if fpass(2) > 0
                    data = bead_hp_allan_var(data, 'yCentresM', ...
                        fpass(2), cropTHPval, false);
                    yStiff = calcStiffness(data.pro.yCentresHP);
                else
                    yStiff = calcStiffness(data.pro.yCentresM(centresRow,:));
                end
                
                data.pro.stiffXY = [xStiff, yStiff];
                accumulated{dayIdx}{1,cellIdx}(fIdx).stiff = data.pro.stiffXY;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fpass = fpass;
                accumulated{dayIdx}{1,cellIdx}(fIdx).opts = data.opts;
                
                % Force normMSD to run processing (hacky but it works)
                fr = data.opts.forceRun;
                data.opts.forceRun = true;
                data = bead_normMSD_polyfit(data, msdDim, msdOffset, msdNumT, true, false, centresRow, false, false, 'CentresHP');
                % Restore previous setting
                data.opts.forceRun = fr;
                
                if isfield(accumulated{dayIdx}{1,cellIdx}(fIdx), 'msd')
                    accumulated{dayIdx}{1,cellIdx}(fIdx).msdRaw = accumulated{dayIdx}{1,cellIdx}(fIdx).msd;
                end
                accumulated{dayIdx}{1,cellIdx}(fIdx).msd = data.pro.amsdObj;
            else
                xStiff = calcStiffness(data.pro.xCentresM(1,:));
                yStiff = calcStiffness(data.pro.yCentresM(1,:));
                data.pro.stiffXY = [xStiff, yStiff];
                accumulated{dayIdx}{1,cellIdx}(fIdx).stiff = data.pro.stiffXY;
            end
            
            % Get tRange from struct
            tR = Range_getter(tRangesLow, dayDirs{dayIdx}, cellIdx, fIdx);
            % Measure the gradient (this function is a bit overkill but
            % it's easiest to implement)
            [tC, fps] = msd_cornerator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), tR);
            % Store the results
            accumulated{dayIdx}{1,cellIdx}(fIdx).fits = fps;
            
            % Save if requested
            if saveDataPro
                save(dataFile, 'data')
            end
            
        end
    end
end

if saveAccu
    save(accuFile, 'accumulated', 'dayDirs', '-v7.3')
end

function Range = Range_getter(Ranges, day, cIdx, fIdx)
if ~isstruct(Ranges)
    tmp = whos('wRanges');
    error('wRanges needs to be a struct, instead got: %s\n', tmp.class)
end
day = ['d' day];
if isfield(Ranges, day)
    c = ['c' num2str(cIdx)];
    if isfield(Ranges.(day), c)
        if size(Ranges.(day).(c), 1) >= fIdx
            Range = Ranges.(day).(c)(fIdx,:);
        else
            Range = {};
        end
    else
        Range = {};
    end
else
    Range = {};
end
end

function [wRanges] = get_wRanges_struct(dayDirs)
% This has numbers from LatB data
wRanges = struct(strcat('d',dayDirs{end}), struct('c1', {{}}));

%% wRanges with 1 low-frequency corner
wRanges.d2021_07_27.c2 = {...
    {[2e-2 1]} {[2e-2 1]};
    {[2e-2 1]} {[2e-2 1]};
    {[2e-2 0.8]} {[8 80]};
    {[2e-2 0.8]} {[8 80]};
    {[2 10]} {[2 10]};
    {[2 10]} {[2 10]};
    {[1 7]} {[1 10]};
    {[1e-2 1.5e-1]} {[2 10]};};
wRanges.d2021_07_27.c1 = {...
    {[1e-2 2]} {[1e-2 5]};
    {[1e-2 5]} {[1e-2 6]};
    {[1e-2 2]} {[1e-2 1]};
    {[1e-2 2]} {[1e-1 2]};
    {[1e-2 1]} {[1e-2 1]};
    {[1e-1 3]} {[1e-2 1]};
    {[1e-2 1]} {[1e-2 1]};
    {[1e-2 1]} {[1e-2 1]};
    {[1e-2 1]} {[1e-2 1]}};
wRanges.d2021_07_26.c1 = {...
    {[1e-2 10]} {[1e-2 10]};
    {[1e-2 5]} {[1e-2 6]};
    {[1e-2 10]} {[1e-2 10]};
    {[6e-2 2]} {[6e-2 2]};
    {[1e-2 1]} {[1e-2 1]};
    {[1e-1 3]} {[1e-2 3]};
    {[1e-1 3]} {[1e-2 3]};
    {[1e-1 3]} {[1e-2 3]}};
wRanges.d2021_07_23.c3 = {...
    {[100 1e4]} {[20 8e2]};
    {[100 1e4]} {[20 8e2]};
    {[100 1e4]} {[20 8e2]};
    {[100 1e4]} {[20 8e2]};
    {[100 1e4]} {[20 8e2]};
    {[100 1e4]} {[20 8e2]};
    {[1e-2 0.8]} {[100 1e4]};
    {[100 1e4]} {[20 8e2]};
    {[100 1e4]} {[20 1e3]}};
wRanges.d2021_07_23.c2 = {...
    {[1e-2 1]} {[5 80]};
    {[3 50]} {[2 30]};
    {[3 50]} {[0.5 5]};
    {[3 50]} {[0.5 5]};
    {[0.5 4]} {[0.5 4]};
    {[3e-1 5]} {[3e-1 3]};
    {[1 10]} {[1 10]}; %
    {[2e-1 10]} {[2e-1 8]};
    {[4e-1 3]} {[1e-1 3]};
    {[1 5]} {[1 5]};};
wRanges.d2021_07_23.c1 = {...
    {[1e-1 5]} {[1e-1 5]};
    {[5e-1 500]} {[1e-1 5]};
    {[1e-1 5]} {[1e-1 5]};
    {[1e-1 3]} {[1e-1 3]};
    {[1e-1 3]} {[1e-1 3]};
    {[1e-1 5]} {[1e-1 5]};
    {[1e-2 2]} {[1e-2 2]}; %
    {[1e-2 30]} {[1e-2 8]};
    {[1e-2 1]} {[1e-2 2]};
    {[1e-2 1]} {[1e-2 2]}};
wRanges.d2021_05_17.c2 = {
    {[0.1 1]} {[0.1 1]}
    {[0.1 1]} {[0.1 1]}
    {[0.1 1]} {[0.1 1]}
    {[0.1 1]} {[0.1 1]}
    {[0.1 1]} {[0.1 1]}};

wRanges.d2021_05_17.c1 = {
    {} {}
    {} {}
    {[0.1 1]} {[0.1 1]}
    {[0.1 1]} {[0.1 1]}
    {[0.1 1]} {[0.1 1]}};

wRanges.d2021_07_16.c1 = {
    {} {}
    {[0.1 1]} {[0.1 1]}
    {[0.1 1]} {[0.1 1]}
    {[0.1 1]} {[0.1 1]}
    {[0.1 1]} {[0.1 1]}
    {[0.1 1]} {[0.1 1]}
    {[0.1 1]} {[0.1 1]}};

wRanges.d2021_07_13.c1 = {...
    {[0.1 3]} {[0.1 2]}
    {[0.1 3]} {[0.1 2]}
    {[0.1 3]} {[0.1 2]}
    {[0.1 3]} {[0.1 2]}
    {[0.1 3]} {[0.1 2]}
    {[0.1 3]} {[0.1 2]}
    {[0.1 1]} {[0.1 2]}
    {[0.1 3]} {[0.1 2]}
    {[0.1 3]} {[0.1 2]}
    {[0.1 3]} {[0.1 2]}
    {[0.1 3]} {[0.1 2]}
    {[0.1 3]} {[0.1 2]}};
%% wRanges with 2 corners
% wRanges.d2021_07_27.c2 = {...
%     {[2e-2 1] [5 1e2]} {[2e-2 1] [5 1e2]};
%     {[2e-2 1] [5 1e2]} {[2e-2 1] [5 1e2]};
%     {[2e-2 0.8] [8e2 2e3]} {[8 80] []};
%     {[2e-2 0.8] [8e2 2e3]} {[8 80] []};
%     {[2 10]} {[2 10]};
%     {[2 10]} {[2 10]};
%     {[1 7]} {[1 10]};
%     {[1e-2 1.5e-1] [30 1e4]} {[2 10]};};
% wRanges.d2021_07_27.c1 = {...
%     {[1e-2 2] [200 2e3]} {[1e-2 5] [60 8e2]};
%     {[1e-2 5] [500 5e4]} {[1e-2 6] [200 1e4]};
%     {[1e-2 2] [8e2 2e3]} {[1e-2 1] [200 1e4]};
%     {[1e-2 2] [8e2 2e3]} {[1e-1 2] [200 1e4]};
%     {[1e-2 1] [8e2 2e3]} {[1e-2 1] [200 1e4]};
%     {[1e-1 3] [20 1e2]} {[1e-2 1] [200 1e4]};
%     {[1e-2 1] [6e2 2e3]} {[1e-2 1] [200 1e4]};
%     {[1e-2 1] [6e2 2e3]} {[1e-2 1] [200 1e4]};
%     {[1e-2 1] [6e2 2e3]} {[1e-2 1] [200 1e4]}};
% wRanges.d2021_07_26.c1 = {...
%     {[1e-2 10] [50 3e2]} {[1e-2 10] [20 8e2]};
%     {[1e-2 5] [50 5e2]} {[1e-2 6] [20 8e2]};
%     {[1e-2 10] [50 3e2]} {[1e-2 10] [20 8e2]};
%     {[6e-2 2] [50 3e2]} {[6e-2 2] [20 8e2]};
%     {[1e-2 1] [50 3e2]} {[1e-2 1] [20 8e2]};
%     {[1e-1 3] [20 1e2]} {[1e-2 3] [20 2e2]};
%     {[1e-1 3] [50 3e2]} {[1e-2 3] [50 3e2]};
%     {[1e-1 3] [10 2e2]} {[1e-2 3] [10 2e2]}};
% wRanges.d2021_07_23.c3 = {...
%     {[100 1e4]} {[20 8e2]};
%     {[100 1e4]} {[20 8e2]};
%     {[100 1e4]} {[20 8e2]};
%     {[100 1e4]} {[20 8e2]};
%     {[100 1e4]} {[20 8e2]};
%     {[100 1e4]} {[20 8e2]};
%     {[1e-2 0.8] [100 1e4]} {[100 1e4]};
%     {[100 1e4]} {[20 8e2]};
%     {[100 1e4]} {[20 1e3]}};
% wRanges.d2021_07_23.c2 = {...
%     {[1e-2 1] [70 1e4]} {[5 80]};
%     {[3 50]} {[2 30]};
%     {[3 50]} {[0.5 5]};
%     {[3 50]} {[0.5 5]};
%     {[0.5 4]} {[0.5 4]};
%     {[3e-1 5]} {[3e-1 3]};
%     {[1 10]} {[1 10]}; %
%     {[2e-1 10]} {[2e-1 8]};
%     {[4e-1 3]} {[1e-1 3]};
%     {[1 5]} {[1 5]};};
% wRanges.d2021_07_23.c1 = {...
%     {[1e-1 5]} {[1e-1 5]};
%     {[5e-1 500]} {[1e-1 5]};
%     {[1e-1 5]} {[1e-1 5]};
%     {[1e-1 3]} {[1e-1 3]};
%     {[1e-1 3]} {[1e-1 3]};
%     {[1e-1 5]} {[1e-1 5]};
%     {[1e-2 2]} {[1e-2 2]}; %
%     {[1e-2 30]} {[1e-2 8]};
%     {[1e-2 1] [20 200] [300 1e3]} {[1e-2 2]};
%     {[1e-2 1] [20 200] [300 1e3]} {[1e-2 2]}};

end

function [tRanges] = get_tRanges_struct(dayDirs)
% This has numbers from LatB data
tRanges = struct(strcat('d',dayDirs{end}), struct('c1', {{}}));

%% Long time corners

tRanges.d2021_05_17.c2 = {
    {[1e-2 0.5] [4 20]} {[1e-2 1] [4 20]}
    {[1e-2 0.5] [4 20]} {[1e-2 1] [4 20]}
    {[1e-2 0.5] [4 20]} {[1e-2 1] [4 20]}
    {[1e-2 0.5] [4 20]} {[1e-2 1] [4 20]}
    {[1e-2 0.5] [3 20]} {[1e-2 1] [4 20]}
    {[1e-2 0.5] [3 20]} {[1e-2 1] [4 20]}};

tRanges.d2021_05_17.c1 = {
    {} {}
    {} {}
    {[1e-1 1] [6 20]} {}
    {[1e-1 1] [6 20]} {[1.4e-2 1] [5 20]}
    {} {}};

tRanges.d2021_07_16.c1 = {
    {} {}
    {[3e-1 3] [10 20]} {[3e-1 2] [10 20]}
    {[2e-1 3/5] [10 20]} {[3e-1 2] [10 20]}
    {[2e-2 1] [5 20]} {[2e-2 1] [10 20]}
    {[2e-2 0.25] [5 20]} {[2e-2 1] [10 20]}
    {[2e-2 0.25] [5 20]} {[2e-2 1] [10 20]}
    {[2e-2 0.25] [5 20]} {}};

tRanges.d2021_07_13.c1 = {...
    {[1e-2 2e-1] [5 20]} {[1e-2 2e-1] [5 20]}
    {[1e-2 2e-1] [5 20]} {[1e-2 2e-1] [5 20]}
    {[1e-2 2e-1] [5 20]} {[1e-2 2e-1] [5 20]}
    {[1e-2 2e-1] [5 20]} {[1e-2 2e-1] [5 20]}
    {[1e-2 2e-1] [5 20]} {[1e-2 2e-1] [5 20]}
    {[1e-2 2e-1] [5 20]} {[1e-2 2e-1] [5 20]}
    {[1e-2 2e-1] [5 20]} {}
    {[1e-2 2e-1] [5 20]} {[1e-2 2e-1] [5 20]}
    {[1e-2 2e-1] [5 9]} {[1e-2 2e-1] [5 20]}
    {[1e-2 2e-1] [5 9]} {[1e-2 2e-1] [5 20]}
    {[1e-2 2e-1] [5 20]} {}
    {[1e-2 2e-1] [5 20]} {}};

tRanges.d2021_07_23.c3 = {...
    {[1e-2 1] [10 30]} {}
    {[1e-2 1] [10 30]} {[1e-1 2] [10 30]}
    {[1e-2 1] [10 30]} {[1e-1 2] [10 30]}
    {[1e-2 1] [10 30]} {[1e-1 2] [10 30]}
    {[1e-2 1] [10 30]} {[1e-1 2] [10 30]}
    {[1e-2 1] [10 30]} {}
    {[1e-2 1] [10 30]} {}
    {[1e-2 1] [5 30]} {[1e-1 2] [5 30]}
    {[1e-2 1] [10 30]} {}};
tRanges.d2021_07_23.c2 = {...
    { } { }
    { } { }
    { } { }
    { } { }
    { } { }
    { } { }
    { } { }
    { } { }
    { } { }
    { } { }};

tRanges.d2021_07_23.c1 = {...
    { } { }
    { } { }
    { } { }
    { } { }
    { } { }
    { } { }
    { } { }
    { } { }
    { } { }
    { } { }};

tRanges.d2021_07_26.c1 = {...
    {[0.08 0.3] [3 20]} {[0.05 0.3] [3 20]}
    {[0.08 0.3] [3 10]} {[0.05 0.3] [3 10]}
    {[0.08 0.3] [3 20]} {[0.05 0.3] [3 20]}
    {[0.1 0.8] [5 30]} {[0.1 0.8] [5 30]}
    {[0.1 0.8] [5 30]} {[0.1 0.8] [5 30]}
    {[0.1 0.8] [5 15]} {[0.1 0.8] [5 30]}
    {[0.1 0.8] [5 15]} {[0.1 0.8] [5 30]}
    {[0.2 1] [3 15]} {[0.2 0.8] [7 30]}};

tRanges.d2021_07_27.c2 = {...
    {[0.2 3] [20 50]} {[0.1 1] [20 50]}
    {[0.4 1.4] [8 20]} {[0.4 1.4] [8 20]}
    {[0.4 2] [8 40]} { }
    {[0.4 4] [15 40]} { }
    { } { }
    { } { }
    { } { }
    {[5 20] [70 200]} { }};

tRanges.d2021_07_27.c1 = {...
    {[3e-2 1/2] [2 7]} {[3e-2 1/2] [2 7]}
    {[1e-2 0.1] [4 20]} {[1e-2 0.4] [5 20]}
    {[1e-2 0.1] [4 20]} {[1e-2 0.4] [5 20]}
    {[1e-2 0.1] [8 30]} {[1e-2 0.4] [8 30]}
    {[1e-2 1] [7 30]} {[1e-2 1] [7 30]}
    {[1e-2 1] [7 30]} {[1e-2 1] [7 30]}
    {[1e-2 1] [7 30]} {[1e-2 1] [7 30]}
    {[1e-2 1] [7 30]} {[1e-2 1] [7 30]}
    {[1e-2 1] [7 30]} {[1e-2 1] [7 30]}};

end
function [tRanges] = get_low_tRanges_struct(dayDirs)
% This has numbers from LatB data
tRanges = struct(strcat('d',dayDirs{end}), struct('c1', {{}}));

%% Time ranges for elastic plateau onset

tRanges.d2021_05_17.c2 = {
    {[1e-2 1e-1] [1 100]} {[1e-2 1e-1] [1 100]}
    {[1e-2 1e-1] [1 100]} {[1e-2 1e-1] [1 100]}
    {[1e-2 1e-1] [1 100]} {[1e-2 1e-1] [1 100]}
    {[1e-2 1e-1] [1 100]} {[1e-2 1e-1] [1 100]}
    {[1e-2 1e-1] [1 100]} {[1e-2 1e-1] [1 100]}};

tRanges.d2021_05_17.c1 = {
    {[1e-2 1] [10 100]} {[1e-2 0.5] [10 100]}
    {[3e-3 0.08] [1 10]} {[3e-3 0.08] [1 100]}
    {[1e-2 1e-1] [1 100]} {[1e-2 1e-1] [1 100]}
    {[1e-2 1e-1] [1 100]} {[1e-2 1e-1] [1 100]}
    {[1e-2 1e-1] [1 100]} {[1e-2 1e-1] [1 100]}
    {[1e-2 1e-1] [1 100]} {[1e-2 1e-1] [1 100]}};

tRanges.d2021_07_16.c1 = {
    {[1e-4 1e-2] [2e-1 2]} {[1e-4 1e-2] [10e-2 100]}
    {[1e-4 1e-2] [2e-1 2]} {[1e-4 1e-2] [10e-2 100]}
    {[1e-4 1e-2] [2e-1 2]} {[1e-4 1e-2] [10e-2 100]}
    {[1e-4 1e-2] [2e-1 2]} {[1e-4 1e-2] [10e-2 100]}%
    {[1e-4 1e-2] [2e-1 2]} {[1e-4 5e-3] [10e-2 100]}
    {[1e-4 1e-2] [2e-1 2]} {[1e-4 5e-3] [10e-2 100]}
    {[1e-4 1e-2] [2e-1 2]} {[1e-4 5e-3] [10e-2 100]}};

tRanges.d2021_07_13.c1 = {...
    {[1e-2 4e-1] [1 100]} {[1e-2 4e-1] [1 100]}
    {[1e-2 2e-1] [1 100]} {[1e-2 2e-1] [1 100]}
    {[1e-2 2e-1] [1 100]} {[1e-2 2e-1] [1 100]}
    {[1e-2 2e-1] [1 100]} {[1e-2 2e-1] [1 100]}
    {[1e-2 2e-1] [1 100]} {[1e-2 2e-1] [1 100]}
    {[1e-2 2e-1] [1 100]} {[1e-2 2e-1] [1 100]}
    {[1e-2 2e-1] [1 100]} {[1e-2 2e-1] [1 100]}
    {[1e-2 2e-1] [1 100]} {[1e-2 2e-1] [1 100]}
    {[1e-2 2e-1] [1 100]} {[1e-2 2e-1] [1 100]}
    {[1e-2 2e-1] [1 100]} {[1e-2 2e-1] [1 100]}
    {[1e-2 2e-1] [1 100]} {[1e-2 2e-1] [1 100]}
    {[1e-2 2e-1] [1 100]} {[1e-2 2e-1] [1 100]}};

tRanges.d2021_07_23.c3 = {...
    {} {}
    {} {}
    {} {}
    {} {}
    {} {}
    {} {}
    {[1e-2 1e-1] [5 50]} {}
    {[1e-2 1e-1] [5 50]} {}
    {[1e-2 1e-1] [5 50]} {}};

tRanges.d2021_07_23.c2 = {...
    {[1e-4 1e-2] [0.5 20]} {[1e-4 1e-2] [0.1 20]}
    {[1e-3 1e-2] [0.1 2]} {[1e-4 1e-2] [1 20]}
    {[1e-4 4e-3] [0.1 100]} {[1e-3 5e-2] [2 20]}
    {[1e-4 4e-3] [0.1 100]} {[1e-3 5e-2] [2 20]}
    {[1e-3 1e-2] [1 20]} {[1e-3 5e-2] [2 20]}
    {[1e-3 5e-2] [1 40]} {[1e-3 5e-2] [5 40]}
    {[1e-4 1e-2] [1 40]} {[1e-4 1e-2] [1 40]}
    {[1e-3 1e-1] [3 40]} {[1e-3 5e-2] [5 40]}
    {[1e-3 1e-1] [3 40]} {[1e-3 5e-2] [3 40]}
    {[1e-3 1e-1] [3 40]} {[1e-3 5e-2] [3 40]}};

tRanges.d2021_07_23.c1 = {...
    {[1e-4 5e-2] [1 50]} {[1e-4 5e-2] [1 50]}
    {[1e-4 5e-2] [1 50]} {[1e-4 5e-2] [1 50]}
    {[1e-4 5e-2] [1 50]} {[1e-4 5e-2] [1 50]}
    {[1e-4 5e-2] [1 50]} {[1e-4 5e-2] [1 50]}
    {[1e-4 5e-2] [1 50]} {[1e-4 5e-2] [1 50]}
    {[1e-4 5e-2] [1 50]} {[1e-4 5e-2] [1 50]}
    {[1e-4 5e-2] [1 50]} {[1e-4 5e-2] [1 50]}
    {[1e-4 1e-2] [1 50]} {[1e-4 5e-2] [1 50]}
    {[1e-4 5e-2] [0.1 0.8] [1 50]} {[1e-4 5e-2] [1 50]}
    {[1e-4 1e-2] [0.04 0.4] [1 50]} {[1e-4 5e-2] [1 50]}};

tRanges.d2021_07_26.c1 = {...
    {[1e-4 1e-2] [1 100]} {[1e-4 0.4e-2] [1 100]}
    {[1e-4 1e-2] [1 100]} {[1e-4 0.4e-2] [1 100]}
    {[1e-4 1e-2] [1 100]} {[1e-4 0.4e-2] [1 100]}
    {[1e-4 1e-2] [1 100]} {[1e-4 0.4e-2] [1 100]}
    {[1e-4 1e-2] [1 100]} {[1e-4 0.4e-2] [1 100]}
    {[1e-4 1e-2] [1 100]} {[1e-4 0.4e-2] [1 100]}
    {[1e-4 1e-2] [1 100]} {[1e-4 0.4e-2] [1 100]}
    {[1e-4 1e-2] [1 100]} {[1e-4 0.4e-2] [1 100]}};

tRanges.d2021_07_27.c2 = {...
    {[1e-4 1e-2] [5 50]} {[3e-4 1e-2] [5 50]}
    {[1e-4 1e-2] [5 50]} {[3e-4 1e-2] [5 50]}
    {[1e-4 1e-2] [5 50]} {[3e-4 1e-2] [5 50]}
    {[1e-4 1e-2] [5 50]} {[3e-4 1e-2] [5 50]}
    {[1e-4 1e-2] [5 50]} {[3e-4 1e-2] [5 50]}
    {[1e-4 1e-2] [5 50]} {[3e-4 1e-2] [5 50]}
    {[1e-4 1e-2] [5 50]} {[3e-4 1e-2] [5 50]}
    {[1e-4 1e-2] [5 50]} {[3e-4 1e-2] [5 50]}};

tRanges.d2021_07_27.c1 = {...
    {[1e-2 1e-1] [1 50]} {[1e-2 1e-1] [1 50]}
    {[1e-2 1e-1] [1 50]} {[1e-2 1e-1] [1 50]}
    {[1e-2 1e-1] [1 50]} {[1e-2 1e-1] [1 50]}
    {[1e-2 1e-1] [1 50]} {[1e-2 1e-1] [1 50]}
    {[1e-2 1e-1] [1 50]} {[1e-2 1e-1] [1 50]}
    {[1e-2 1e-1] [1 50]} {[1e-2 1e-1] [1 50]}
    {[1e-2 1e-1] [1 50]} {[1e-2 1e-1] [1 50]}
    {[1e-2 1e-1] [1 50]} {[1e-2 1e-1] [1 50]}
    {[1e-2 1e-1] [1 50]} {[1e-2 1e-1] [1 50]}};

end