%% Process multiple sets of bead data sequentially
%
% bead_processing_accumulator_v2 - do normalisations and shiz
%
% This file accumulates: trapped bead data
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
masterDir = '../';           % Where all the data is
saveDir = 'MAMRheo_reprocess/';        % Where to save processed data

% Which day's data to load
dayDirs = {'2022_06_21'};    

% This cell contains 1 cell per dayDir above. Within that needs to be
% indexes for that day's dirList, to choose the correct datasets.
setIdxs = {
    ... % Sets from 2022_06_21
        {[4 5 6 9 12]} ...
        }; 

% Ain't nobody got time (be)for this!
% Later I will record time in a filename when recording data
thymes = { 
    ... % Times from 2022_06_21
        {[0 0 0 0 0]} ...
    };
    
wRanges = get_wRanges_struct(dayDirs);
tRanges = get_tRanges_struct(dayDirs);
tRangesLow = get_low_tRanges_struct(dayDirs);
wRangesLow = get_low_wRanges_struct(dayDirs);
%% Experiment parameters
mPerPx = 0.065e-6;           % Camera pixel size calibration
ignoreDirs = {'focal_sweep_with_bead'}; % Directories to ignore (ones without data)

%% Processing parameters
doPro = true;      % Actually do the processing (useful if you've not found corners yet)

cropThack = []; % I should create another variable to store which sets cropThack will be applied to
cropTs = {[]};
fitPoly = 1; % Fit a polynomial to remove drift.
fitPolyOrder = 0;      % Order of polynomial to be fitted

timeReg = true;         % Fix non-uniform time vector (caused by new version of fast_acq)
angleCorrection = false; % Transform to (r, rθ) co-ordinates
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
forceRun = false;        % Don't try to take data from file and reuse as much as possible
saveDataPro = false;    % Save processed data to file (probably only makes things slower for now)
saveDataRaw = true;     % Save raw data to file (should speed up loading)
saveAccu = true;        % Save data to file
dataSuff = '_simple';       % Suffix for filename when saving/loading
accuFile = 'accumulated_trap';

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
                msd = cat(3,accumulated{dayIdx}{1,cellIdx}(fIdx).msd.msd{:});
                % Store max time
                accumulated{dayIdx}{1,cellIdx}(fIdx).tMax = msd(end,1);
                
                % G'_0 proportional to 1/MSD at gradient minimum.
                p = nan(2,1,3);
                for dim = 1:2
                    % Calculate derivative (discard end of MSD first)
                    tau = msd(:,1,dim);
                    tau = tau(tau < max(tau)/20);
                    m = msd(1:length(tau),2,dim);
                    [dydx, tout] = msd_gradientor(tau, m, 'lsq', 15);
                    % Take minima
                    [~, idx] = min(dydx);
                    [~, ind] = min(abs(tau - tout(idx)));
                    p(dim,1,:) = [2*kBT(293)./m(ind) dydx(idx) tout(idx)];
                end
                % Store G'_0 and gradient minimum 
                accumulated{dayIdx}{1,cellIdx}(fIdx).stiff2 = p;
                
                % Get wRange from struct
                wR = Range_getter(wRangesLow, dayDirs{dayIdx}, cellIdx, fIdx);

                if isempty(wR) || isempty([wR{:}])
                    warning('Empty wR on %s', data.dirPath)
                end
                
                % Do the fourier transform to find intercept frequency
                [FT, oC] = msd_fourier_transformator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), ...
                    'wRange',wR, 'trunc', 'FF','show_int',true, 'doPlot', false);
                
                % Extract and store the corners
                oC = cat(1,oC{:});
                accumulated{dayIdx}{1,cellIdx}(fIdx).oC = oC;
                
                % Get tRange from struct
                tR = Range_getter(tRangesLow, dayDirs{dayIdx}, cellIdx, fIdx);
                
                if isempty(tR) || isempty([tR{:}])
                    warning('Empty tR on %s', data.dirPath)
                end
                
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
                    'normT', tnorm, 'normR', 2*kBT(293)./p(1:2));
                accumulated{dayIdx}{1,cellIdx}(fIdx).tCtauG = tC3;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fpstauG = fps3;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrtauG = fitErr3;
                
                
                % Get space-normalized subdiffusion parameters
                [tC31, fps31, fitErr31] = msd_cornerator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), tR, ...
                    'normR', 2*kBT(293)./p(1:2));
                accumulated{dayIdx}{1,cellIdx}(fIdx).tCG = tC31;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fpsG = fps31;
                accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrG = fitErr31;
                
                %% Long-time corner, time domain first
                tR = Range_getter(tRanges, dayDirs{dayIdx}, cellIdx, fIdx);
                
                if isempty(tR) || isempty([tR{:}])
                    warning('Empty tR on %s', data.dirPath)
                end
                
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
                
                if isempty(wR) || isempty([wR{:}])
                    warning('Empty wR on %s', data.dirPath)
                end
                
                % Do the fourier transform to find intercept frequency
                [FT, oC] = msd_fourier_transformator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), ...
                    'wRange',wR, 'trunc', 'FF','show_int',true, 'doPlot',false);
                
                % Extract and store the corners
                oC = cat(1,oC{:});
                % Average the frequency and time domain measurements
                tnorm = mean([tC4; 1./oC'],'omitnan');
                accumulated{dayIdx}{1,cellIdx}(fIdx).tnormH = tnorm;
                
                if ~all(isnan(tnorm))
                    % Use the corners to get long-time normalized fits
                    [tC5, fps5, fitErr5] = msd_cornerator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), tR, ...
                        'normT', tnorm);
                    accumulated{dayIdx}{1,cellIdx}(fIdx).tCtauH = tC5;
                    accumulated{dayIdx}{1,cellIdx}(fIdx).fpstauH = fps5;
                    accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrtauH = fitErr5;
                    
                    % Get time-and-space-normalized superdiffusion parameters
                    [tC6, fps6, fitErr6] = msd_cornerator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), tR, ...
                        'normT', tnorm, 'normR', 2*kBT(293)./p(1:2));
                    accumulated{dayIdx}{1,cellIdx}(fIdx).tCtauGH = tC6;
                    accumulated{dayIdx}{1,cellIdx}(fIdx).fpstauGH = fps6;
                    accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrtauGH = fitErr6;
                    
                    % Get space-normalized superdiffusion parameters
                    [tC61, fps61, fitErr61] = msd_cornerator(data.pro.amsdObj, accumulated{dayIdx}{2,cellIdx}(fIdx), tR, ...
                        'normR', 2*kBT(293)./p(1:2));
                    accumulated{dayIdx}{1,cellIdx}(fIdx).tCGH = tC61;
                    accumulated{dayIdx}{1,cellIdx}(fIdx).fpsGH = fps61;
                    accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrGH = fitErr61;
                else
                    tcNaN = nan(1,2);
                    fNaN = nan(2,2,2);
                    accumulated{dayIdx}{1,cellIdx}(fIdx).tCtauH = tcNaN;
                    accumulated{dayIdx}{1,cellIdx}(fIdx).fpstauH = fNaN;
                    accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrtauH = fNaN;
                    
                    accumulated{dayIdx}{1,cellIdx}(fIdx).tCtauGH = tcNaN;
                    accumulated{dayIdx}{1,cellIdx}(fIdx).fpstauGH = fNaN;
                    accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrtauGH = fNaN;
                    
                    accumulated{dayIdx}{1,cellIdx}(fIdx).tCGH = tcNaN;
                    accumulated{dayIdx}{1,cellIdx}(fIdx).fpsGH = fNaN;
                    accumulated{dayIdx}{1,cellIdx}(fIdx).fiterrGH = fNaN;
                end
                
                % Save if requested
                if saveDataPro
                    save(dataFile, 'data')
                end
            else
%                 warning('Did not process data because doPro is false')
            end
        end
    end
end

if saveAccu
    accumulated = [accumulated; dayDirs];
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
% This hasn't numbers 
wRanges = struct(strcat('d',dayDirs{end}), struct('c1', {{}}));
%% wRanges with 1 high-frequency corner - Free diffusion end/elastic plateau onset

wRanges.d2022_06_21.c1 = {
    {[1 10]} {[1 100]}
    {[1 10]} {[1 100]}
    {[1 10]} {[1 100]}
    {[1 10]} {[1 100]}
    {[1 10]} {[1 100]}};


end


function [tRanges] = get_low_tRanges_struct(dayDirs)
% This hasn't numbers 
tRanges = struct(strcat('d',dayDirs{end}), struct('c1', {{}}));

%% tRanges with 1 low time corner - Free diffusion end/elastic plateau onset
tRanges.d2022_06_21.c1 = {
    {[1e-3 1e-1] [2 20]} {[1e-3 1e-1] [2 20]} ; 
    {[1e-3 5e-2] [2 20]} {[1e-3 5e-2] [2 20]} ; 
    {[1e-3 2e-2] [2 20]} {[1e-3 2e-2] [2 20]} ; 
    {[1e-3 1e-2] [2 20]} {[1e-3 1e-2] [2 20]} ;
    {[1e-3 1e-2] [0.1 5]} {[1e-3 1e-2] [0.1 5]}
    };
    
end

function [wRanges] = get_wRanges_struct(dayDirs)
% This hasn't numbers 
wRanges = struct(strcat('d',dayDirs{end}), struct('c1', {{}}));

%% wRanges with 1 low-frequency corner - Superdiffuse onset
wRanges.d2022_06_21.c1 = {
    {} {}
    {} {}
    {} {}
    {} {}
    {} {}
    };
end

function [tRanges] = get_tRanges_struct(dayDirs)
% This hasn't numbers 
tRanges = struct(strcat('d',dayDirs{end}), struct('c1', {{}}));

%% tRanges with 1 long-time corner - Superdiffuse onset
tRanges.d2022_06_21.c1 = {
    {} {}
    {} {}
    {} {}
    {} {}
    {} {}
    };


end
