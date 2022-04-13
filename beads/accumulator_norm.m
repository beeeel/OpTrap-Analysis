%% Process multiple sets of bead data sequentially
%
% V1: Perform iterative short time processing as per Manlio's suggestion:
%       Perform initial processing as normal, then find time of long time
%       uptick/low frequency intercept. Use this to apply highpass
%       filtering before recalculating MSD and FTs.
%
% This file accumulates: Good Latrunculin B data
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
% Copy data.nPoints and exposure time
%  data.metadata.FrameKey_0_0_0.UserData.Camera_1_Exposure.scalar into
%  accumulated

%% Data parameters
close all

% Data storage parameters
masterDir = '../';             % Where all the data is
saveDir = 'accumulator_2022_01_14_latB_good/';        % Where to save processed data

% Which day's data to load
dayDirs = {'2022_02_01','2022_01_18', '2022_01_13','2022_01_11', '2021_07_27' '2021_07_26' '2021_07_23'};

% This cell contains 1 cell per dayDir above. Within that needs to be
% indexes for that day's dirList, to choose the correct datasets.
setIdxs = {
    ... % Sets from 2022_02_01
        {[11 14]} ...
    ... % Sets from 2022_01_18
        {[8 11]} ...
    ... % Sets from 2022_01_13
        {[3 5] [11 13]} ...
    ... % Sets from 2022_01_11
        {[3 6] [10 12]} ...
    ... % Sets from 2021_07_27
        {[11 15]} ...
    ... % Sets from 2021_07_26
        {[6 10]} ...
    ... % Sets from 2021_07_23
        {[2 7] [12 17]} };

% Ain't nobody got time (be)for this!
% Later I will record time in a filename when recording data
thymes = { 
    ... % Times from 2022_02_01
        {[-14 28]} ...
    ... % Times from 2022_01_18
        {[-19 31]} ...
    ... % Times from 2022_01_13
        {[-11 33] [-27 23]} ...
    ... % Times from 2022_01_11
        {[-23 42] [-13 32]} ...
    ... % Times from 2021_07_27
        {[-15 37]} ...
    ... % Times from 2021_07_26
        {[-7 32]} ...
    ... % Times from 2021_07_23
        {[-14 28] [-12 39] }};
    
wRanges = get_wRanges_struct(dayDirs);
tRanges = get_tRanges_struct(dayDirs);
tRangesLow = get_low_tRanges_struct(dayDirs);
wRangesLow = get_low_wRanges_struct(dayDirs);
%% Experiment parameters
mPerPx = 0.065e-6;           % Camera pixel size calibration
ignoreDirs = {'focal_sweep_with_bead'}; % Directories to ignore (ones without data)

%% Processing parameters
doPro = false;      % Actually do the processing (useful if you've not found corners yet)

cropThack = [14e4 1e4]; % I should create another variable to store which sets cropThack will be applied to
cropTs = {[]};
fitPoly = 1; % Fit a polynomial to remove drift.
fitPolyOrder = 0;      % Order of polynomial to be fitted

timeReg = true;         % Fix non-uniform time vector (caused by new version of fast_acq)
angleCorrection = true; % Transform to (r, rθ) co-ordinates
calcStiff = 1;          % Calculate trap stiffness from position variance
stiffIdx = 1;           % Row of centres to store stiffness from
cropTHPval = 0;         % Frames to crop after HP filter
msdOffset = [1];        % Offset from start when taking data to calculate mean-square displacements
loadPics = false;       % Load images from TIFs. May be overriden by angleCorrection
doPlots = false;

msdDim = 'all';         % Direction to calculate MSD in - 'x', 'y', or 'all'
centresRow = 1;         % Row of centres array to use for MSDs (empty for default)
msdNumT = [];           % Number of time points to use for MSDs (empty for all)
msdUseRaw = false;      % Use raw or processed data for MSDs (empty for default)
msdDoNorm = false;      % Normalize MSDs by position variance (empty for default)
doFFT = false;          % Calculate FFT and maybe plot
msdFTsmoothF = [];      % Length of smoothing kernel used before first FT

% Data file parameters
forceRun = true;        % Don't try to take data from file and reuse as much as possible
saveDataPro = false;    % Save processed data to file (probably only makes things slower for now)
saveDataRaw = true;     % Save raw data to file (should speed up loading)
saveAccu = true;        % Save data to file
dataSuff = '_simple';       % Suffix for filename when saving/loading
accuFile = 'accumulated_norm';

% ARE YOU READY??
% accumulated = cell(size(dayDirs));
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

            % Record options and calibration - needs to use a variable to store which sets to apply to
            if false %dayIdx == 3 && cellIdx == 2 && fIdx == 3
                data.opts.cropT = [cropThack(1) data.nPoints];
            elseif false %dayIdx == 6 && cellIdx == 1 && fIdx == 3
                data.opts.cropT = [cropThack(2) data.nPoints];
            else
                data.opts.cropT = cropTs{1};
            end
            data.opts.pOrder = fitPolyOrder*fitPoly;
            data.mPerPx = mPerPx;
            data.opts.angleCorrection = angleCorrection;
            data.opts.centresRow = centresRow;
            data.opts.timeRegularisation = timeReg;

            % Apply scale and do polyfit
            data = bead_preProcessCentres(data);
            
            % Record where the data is coming from
            accumulated{dayIdx}{1,cellIdx}(fIdx).dirPath = data.dirPath;
            accumulated{dayIdx}{1,cellIdx}(fIdx).fName = data.fName;
            accumulated{dayIdx}{1,cellIdx}(fIdx).nP = data.nPoints;
            
            if isfield(data,'metadata')
                if isfield(data.metadata.FrameKey_0_0_0.UserData, 'Camera_1_Exposure')
                    accumulated{dayIdx}{1,cellIdx}(fIdx).exp = str2double(data.metadata.FrameKey_0_0_0.UserData.Camera_1_Exposure.scalar);
                elseif isfield(data.metadata.FrameKey_0_0_0.UserData, 'optiMOSSCMOS_Exposure')
                    accumulated{dayIdx}{1,cellIdx}(fIdx).exp = str2double(data.metadata.FrameKey_0_0_0.UserData.optiMOSSCMOS_Exposure.scalar);
                end
            end
            
            % Record whether Imstack was loaded (for angular correction)
            accumulated{dayIdx}{1,cellIdx}(fIdx).angCorr = isfield(data,'ImstackFullFoV');
            
            % Update times array to show time halfway through acquisition
            accumulated{dayIdx}{2,cellIdx}(fIdx) = accumulated{dayIdx}{2,cellIdx}(fIdx) ...
                + data.raw.timeVecMs(end/2)/60e3;
            %% Process data
            % Look at mean-square displacement (for cell-bead expts)
            if ~isempty(msdOffset)
                data = bead_normMSD(data, 'direction', msdDim, 'offset', msdOffset, ...
                    'numT', msdNumT, 'doPlots', false, 'useRaw', msdUseRaw, ...
                    'doNorm', msdDoNorm);
                accumulated{dayIdx}{1,cellIdx}(fIdx).msd = data.pro.amsdObj;
            end
            if doPro
                % Store max time
                accumulated{dayIdx}{1,cellIdx}(fIdx).tMax = msd(end,1);
                
                % G'_0 proportional to 1/MSD at gradient minimum.
                p = nan(2,1,2);
                for dim = 1:2
                    % Calculate derivative
                    tau = msd(:,1,dim);
                    m = msd(:,2,dim);
                    [dydx, tout] = msd_gradientor(tau, m, 'lsq', 15);
                    % Take minima
                    [~, idx] = min(dydx);
                    [~, ind] = min(abs(tau - tout(idx)));
                    p(dim,1,:) = [2*kBT(293)./m(ind) dydx(idx)];
                end
                % Store G'_0 and gradient minimum 
                accumulated{dayIdx}{1,cellIdx}(fIdx).stiff2 = p;
                
                % Get wRange from struct
                wR = Range_getter(wRangesLow, dayDirs{dayIdx}, cellIdx, fIdx);

                % Do the fourier transform to find intercept frequency
                [FT, oC] = msd_fourier_transformator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), ...
                    'wRange',wR, 'fh', fh, 'trunc', 'FF','show_int',true);
                
                % Extract and store the corners
                oC = cat(1,oC{:});
                accumulated{dayIdx}{1,cellIdx}(fIdx).oC = oC;
                
                % Get tRange from struct
                tR = Range_getter(tRangesLow, dayDirs{dayIdx}, cellIdx, fIdx);
                % Also fit to the MSD to find corner time equivalent frequency
                [tC1, fps1, fitErr1] = msd_cornerator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), tR);
                
                accumulated{dayIdx}{1,cellIdx}(fIdx).tR = tR;
                accumulated{dayIdx}{1,cellIdx}(fIdx).tC = tC1;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fps = fps1;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fiterr = fitErr1;

                if isempty(tC1)
                    tC1 = nan(1,2);
                end
                
                % Use mean of corner time and 1/intercept frequency.
                % There's definitely scope for mistakes to occur if you
                % want to be getting smarter, but that will have to wait.
                tnorm = mean([tC1; 1./oC'],'omitnan');
                accumulated{dayIdx}{1,cellIdx}(fIdx).tnorm = tnorm;
                
                % This has assumed you've only given one intercept
                % frequency range and one corner range (two time ranges)
                
                % Get the time-normalized subdiffusion parameters
                [tC2, fps2, fitErr2] = msd_cornerator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), tR, ...
                    'normT', tnorm);
                accumulated{dayIdx}{1,cellIdx}(fIdx).tCtau = tC2;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fpstau = fps2;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrtau = fitErr2;
                
                % Get time-and-space-normalized subdiffusion parameters
                [tC3, fps3, fitErr3] = msd_cornerator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), tR, ...
                    'normT', tnorm, 'normR', p(1:2));
                accumulated{dayIdx}{1,cellIdx}(fIdx).tCtauG = tC3;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fpstauG = fps3;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrtauG = fitErr3;
                
                % Long-time corner, time domain first
                tR = Range_getter(tRanges, dayDirs{dayIdx}, cellIdx, fIdx);
                % Measure the gradient (this function is a bit overkill but
                % it's easiest to implement)
                [tC4, fps4, fitErr4] = msd_cornerator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), tR);
                % Store the results
                accumulated{dayIdx}{1,cellIdx}(fIdx).tRH = tR;
                accumulated{dayIdx}{1,cellIdx}(fIdx).tCH = tC4;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fpsH = fps4;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrH = fitErr4;
                
                % Frequency domain
                wR = Range_getter(wRanges, dayDirs{dayIdx}, cellIdx, fIdx);
                % Do the fourier transform to find intercept frequency
                [FT, oC] = msd_fourier_transformator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), ...
                    'wRange',wR, 'fh', fh, 'trunc', 'FF','show_int',true);
                
                % Extract and store the corners
                oC = cat(1,oC{:});
                % Average the frequency and time domain measurements
                tnorm = mean([tC4; 1./oC'],'omitnan');
                accumulated{dayIdx}{1,cellIdx}(fIdx).tnormH = tnorm;
                
                % Use the corners to get long-time normalized fits
                [tC5, fps5, fitErr5] = msd_cornerator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), tR, ...
                    'normT', tnorm);
                accumulated{dayIdx}{1,cellIdx}(fIdx).tCtauH = tC5;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fpstauH = fps5;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrtauH = fitErr5;
                
                % Get time-and-space-normalized superdiffusion parameters
                [tC6, fps6, fitErr6] = msd_cornerator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), tR, ...
                    'normT', tnorm, 'normR', p(1:2));
                accumulated{dayIdx}{1,cellIdx}(fIdx).tCtauGH = tC6;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fpstauGH = fps6;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrtauGH = fitErr6;
                
                % Save if requested
                if saveDataPro
                    save(dataFile, 'data')
                end
            else
                warning('Did not process data because doPro is false')
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
function [wRanges] = get_low_wRanges_struct(dayDirs)
% This has numbers from LatB data
wRanges = struct(strcat('d',dayDirs{end}), struct('c1', {{}}));

%% wRanges with 1 high-frequency corner

wRanges.d2022_02_01.c1 = {
    {[0.1 10]} {[0.1 10]}
    {[0.1 3]} {[0.1 3]}};

wRanges.d2022_01_18.c1 = {
    {[50 300]} {[50 500]}
    {[10 1000]} {[10 1000]}};

wRanges.d2022_01_13.c1 = {
    {[20 150]} {[10 300]}
    {[7 30]} {[10 200]}};

wRanges.d2022_01_13.c2 = {
    {[30 200]} {[20 200]}
    {} {[4 80]}};

wRanges.d2022_01_11.c1 = {
    {} {[5 50]}
    {} {}};

wRanges.d2022_01_11.c2 = {
    {[10 100]} {[1 20]}
    {[10 100]} {[1 20]}};

wRanges.d2021_07_23.c1 = {...
    {[5e-1 500]} {[1e-1 5]};
    {[1 40]} {[1e-2 2]}};

wRanges.d2021_07_23.c2 = {...
    {[3 100]} {[2 30]};
    {[1 20]} {[1 20]}};

wRanges.d2021_07_26.c1 = {
    {[50 5e2]} {[20 8e2]};
    {[20 1e2]} {[20 2e2]}};
% {[1e-2 5] [50 5e2]} {[1e-2 6] [20 8e2]};
%     {[1e-1 3] [20 1e2]} {[1e-2 3] [20 2e2]}};

wRanges.d2021_07_27.c1 = {...
    {[10 100]} {[10 100]}
    {[2 15]} {[1 10]}};
%     {[1e-2 5] [500 5e4]} {[1e-2 6] [200 1e4]};
%     {[1e-2 1] [8e2 2e3]} {[1e-2 1] [200 1e4]}};

end

function [tRanges] = get_low_tRanges_struct(dayDirs)
% This has numbers from LatB data
tRanges = struct(strcat('d',dayDirs{end}), struct('c1', {{}}));

%% Time ranges for elastic plateau onset

tRanges.d2022_02_01.c1 = {
    {[0.1 1] [4 20]} {[0.1 1] [4 20]}
    {[0.1 1] [20 100]} {[0.1 1] [7 40]}};

tRanges.d2022_01_18.c1 = {
    {[2e-3 2e-2] [0.1 1]} {[2e-3 2e-2] [0.1 1]}
    {[2e-3 1e-2] [0.1 1]} {[2e-3 1e-2] [0.1 1]}};

tRanges.d2022_01_13.c1 = {
    {[1e-3 1e-2] [0.1 1]} {[1e-4 1e-2] [0.1 1]}
    {[1e-3 1e-2] [0.2 2]} {[1e-4 1e-2] [0.2 2]}};

tRanges.d2022_01_13.c2 = {
    {[1e-4 3e-3] [6e-2 1]} {[1e-4 3e-3] [6e-2 1]}
    {[2e-3 4e-2] [4e-1 10]} {[2e-3 4e-2] [4e-1 10]}};

tRanges.d2022_01_11.c1 = {
    {[6e-3 1e-1] [4e-1 2]} {[6e-3 1e-1] [3e-1 2]}
    {[6e-3 1e-1] [1 10]} {[5e-3 5e-2] [1 10]}};

tRanges.d2022_01_11.c2 = {
    {[2e-3 2e-2] [1 10]} {[6e-3 6e-2] [1 10]}
    {[3e-3 3e-2] [0.3 3]} {[3e-3 3e-2] [0.3 3]}};
    
tRanges.d2021_07_27.c1 = {...
    {[3e-3 3e-2] [0.1 1]} {[3e-3 3e-2] [0.3 0.9]}
    {[1e-2 1e-1] [1 50]} {[1e-2 1e-1] [1 50]}};

tRanges.d2021_07_26.c1 = {...
    {[1e-4 1e-2] [0.06 0.3]} {[1e-4 0.4e-2] [0.04 0.2]}
    {[1e-4 1e-2] [0.06 0.8]} {[1e-4 0.4e-2] [0.1 1]}};

tRanges.d2021_07_23.c1 = {...
    {[1e-4 5e-2] [8 100]} {[1e-4 5e-2] [4 100]}
    {[1e-4 5e-2] [0.2 50]} {[1e-2 5e-1] [8 80]}};

tRanges.d2021_07_23.c2 = {...
    {[1e-3 1e-2] [0.1 2]} {[1e-4 1e-2] [1 20]}
    {[1e-4 1e-2] [0.7 10]} {[1e-4 1e-2] [0.3 13]}};

end

function [wRanges] = get_wRanges_struct(dayDirs)
% This has numbers from LatB data
wRanges = struct(strcat('d',dayDirs{end}), struct('c1', {{}}));

%% wRanges with 1 low-frequency corner

wRanges.d2022_02_01.c1 = {
    {[0.1 10]} {[0.1 10]}
    {[2e-2 2]} {[2e-2 2]}};

wRanges.d2022_01_18.c1 = {
    {[1e-1 6]} {[1e-1 4]}
    {[1e-2 1]} {[1e-2 1]}};

wRanges.d2022_01_13.c1 = {
    {[6e-1 8]} {[2e-1 2]}
    {[2e-2 1]} {[2e-2 1]}};

wRanges.d2022_01_13.c2 = {
    {[2e-1 3]} {[7e-2 1]}
    {} {[3e-2 7e-1]}};

wRanges.d2022_01_11.c1 = {
    {[3e-2 0.5]} {}
    {} {}};

wRanges.d2022_01_11.c2 = {
    {} {}
    {[3e-2 1]} {[3e-2 1]}};

wRanges.d2021_07_27.c1 = {...
    {[1e-2 5]} {[1e-2 6]};
    {} {}};

wRanges.d2021_07_26.c1 = {...
    {[1e-2 5]} {[1e-2 6]};
    {[1e-1 3]} {[1e-2 3]}};

wRanges.d2021_07_23.c2 = {...
    {[0.01 1]} {[2 30]};
    {} {}};

wRanges.d2021_07_23.c1 = {...
    {} {};
    {} {[1e-2 2]}};

end

function [tRanges] = get_tRanges_struct(dayDirs)
% This has numbers from LatB data
tRanges = struct(strcat('d',dayDirs{end}), struct('c1', {{}}));

%% tRanges used for determining highpass frequency
tRanges.d2022_02_01.c1 = {
    {[0.1 1] [4 20]} {[0.1 1] [4 20]}
    {[0.1 1] [20 100]} {[0.1 1] [7 40]}};

tRanges.d2022_01_18.c1 = {
    {[0.1 1] [10 100]} {[0.1 1] [10 100]}
    {[0.1 1] [60 200]} {[0.1 1] [60 200]}};

tRanges.d2022_01_13.c1 = {
    {[0.1 1] [25 150]} {[0.1 1] [25 150]}
    {[0.2 2] [40 100]} {[0.2 2] [40 100]}};

tRanges.d2022_01_13.c2 = {
    {[6e-2 1] [40 200]} {[6e-2 1] [10 100]}
    {[0.5 10] [120 400]} {[4e-1 10] [100 400]}};

tRanges.d2022_01_11.c1 = {
    {[4e-1 2] [120 400]} {[0.5 10] [60 200]}
    {[1 10] [100 400]}   {[1 10] [100 400]}};

tRanges.d2022_01_11.c2 = {
    {[1 10] [20 50]} {[1 10] [30 100]}
    {[3e-1 2] [40 100]} {[4e-1 6] [90 400]}};
    
tRanges.d2021_07_26.c1 = {...
    {[0.06 0.3] [5 20]} {[0.05 0.3] [3 100]}
    {[0.08 0.8] [5 100]} {[0.09 0.7] [20 100]}};

tRanges.d2021_07_27.c1 = {...
    {[0.2 1] [8 20]} {[0.4 1] [8 20]}
    {[1 50] [50 100]} {}};

tRanges.d2021_07_23.c1 = {...
    {[1e-4 5e-2] [8 100]} {[1e-4 5e-2] [4 100]}
    {[0.2 30] [30 100]} {[0.1 5] [8 80]}};

tRanges.d2021_07_23.c2 = {...
    {[0.1 2] [60 100]} {[1 20] [60 100]}
    {[0.7 10] [60 200]} {[0.3 13] [60 200]}};
end