%% Process multiple sets of bead data sequentially
close all
% Experiment parameters
mPerPx = 0.065e-6;           % Camera pixel size calibration
laserPowers = 0;      % Laser power in % for the datasets used
ignoreDirs = {}; % Directories to ignore (ones without data)

% Processing parameters
angleCorrection = false; % Convert cartesian co-ordinates into polar (r, r.*theta)
timeRegularise = true;  % Regularise time vector (for data acquired with fast_acq_v9 onwards)
cropTs = {[]};
fitPoly = [1]; % Fit a polynomial to remove drift
fitPolyOrder = 1;       % Order of polynomial to be fitted
calcStiff = 1;          % Calculate trap stiffness from position variance
fpass = 0;              % High-pass frequency
fstop = [124.5 125.5];  % Bandstop frequencies - one range per row.
cropTHPval = 1;       % Frames to crop after HP filter
msdOffset = 1;+[0 3e5];          % Offset from start when taking data to calculate mean-square displacements
msdDim = 'all';         % Direction to calculate MSD in - 'x', 'y', or 'all'
msdCentresRow = [];      % Row of centres array to use for MSDs (empty for default)
msdNumT = [];7e5;           % Number of time points to use for MSDs (empty for all)
msdUseRaw = false;      % Use raw or processed data for MSDs (empty for default)
msdDoNorm = false;       % Normalize MSDs by position variance (empty for default)
msdErrorBars = false;    % Plot errorbars on MSD (empty for default)
doFFT = false;           % Calculate FFT and maybe plot

% Data file parameters
forceRun = false;       % Try to take data from file and reuse as much as possible
saveData = false;        % Save data to file
dataSuff = '_120k_28min';       % Suffix for filename when saving/loading

% Plotting parameters
loadPics = false;    % Load the TIF images
saveFigs = false;
showStack = false;   % Open the image data in ImageJ
doPlots = true;      % Plot the centres data
compCentres = false; % Show the Imstack with live calculated and offline calculated centres
setLims = [];[-1 1] * 0.2;     % Set axis limits in um on position plots, empty for auto limit
fftYlim = [0 1];    % Set limits in units amplitude for FFT plots, empty for auto limit
plotDC = false;
plotRaw = false;
setDCLims = [0 0.1];  % Set limits in units pixel brightness for DC plots, empty for auto limit

%while true

% Get all the children directories in a struct
dirList = dir;
dirList = dirList([dirList.isdir]);
dirList = dirList(3:end);
for d = 1:length(ignoreDirs)
    dirList = dirList(~strcmp({dirList.name},ignoreDirs{d}));
end
% Check cropTs has been made with the right number of elements. Throws an
% error and gives an empty list if not.
% checkCropTs(cropTs, dirList);

% Output names
tmp = {dirList.name}'; for idx = 1:length(tmp); tmp{idx} = sprintf('%i: %s', idx, tmp{idx}); end
disp(tmp)

out = struct();
out(1).stiff = nan;
newIdx = find(max([dirList.datenum]) == [dirList.datenum], 1, 'first');
%%
for    fileIdx = [2 5 7 8]%length(dirList)-2
    %% Load and pre-process
    % Either create a new struct or load one named dataFile
    dataFile = [dirList(fileIdx).name '_processed' dataSuff '.mat'];
    if forceRun || ~exist(dataFile, 'file')
        data = struct([]);
        data(1).opts.forceRun = forceRun;
        
        % Set names and load data
        data.dirPath = [dirList(fileIdx).folder '/' dirList(fileIdx).name];
        data.fName = dirList(fileIdx).name;
        data = bead_loadData(data, any([loadPics angleCorrection]));
        
        % Apply calibration and crop time
        data.opts.cropT = cropTs{1};%fileIdx};
        data.opts.pOrder = fitPolyOrder*fitPoly(1);%fileIdx);
        data.mPerPx = mPerPx;
    else
        tmp = load(dataFile, 'data');
        data = tmp.data;
        clear tmp
        data.opts.forceRun = forceRun;
    end
    
    data.opts.timeRegularisation = timeRegularise;
    % Apply scale and do polyfit
    data = bead_preProcessCentres(data);
    %% Process data
    % Calculate the stiffnesses and put into data
    if calcStiff
        stiffIdx = 1;
        out(fileIdx).suffix = data.raw.suffixes{stiffIdx};
        out(fileIdx).name = data.fName;
        
        xStiff = calcStiffness(data.pro.xCentresM);
        yStiff = calcStiffness(data.pro.yCentresM);
        data.pro.stiffXYpro = [xStiff, yStiff];
        out(fileIdx).stiff = data.pro.stiffXYpro(stiffIdx,:);
        
        xStiff = calcStiffness(data.raw.xCentresPx, data.mPerPx);
        yStiff = calcStiffness(data.raw.yCentresPx, data.mPerPx);
        data.pro.stiffXYraw = [xStiff, yStiff];
        out(fileIdx).stiffraw = data.pro.stiffXYraw(stiffIdx,:);
    end
    
    % Plot the processed data
    if doPlots
        if plotRaw
            fh = bead_plotRawData(data, setLims);
        end

        fh = bead_plotProData(data, setLims); %#ok<*UNRCH>
        if saveFigs
            saveas(fh, [data.fName '_pro.png'])
        end
        
        % If there's DC data, plot it
        if isfield(data.raw,'dcAvg') && plotDC
            fh = bead_plotDCData(data, setDCLims);
            if saveFigs
                saveas(fh, [data.fName '_DC.png'])
            end
        end
    end
    %%
    % Calculate frequency spectrum in physical units
    if doFFT
        data = bead_fft_scaled(data, doPlots, fftYlim);
        if saveFigs
            saveas(gcf, [data.fName '_fft.png'])
        end
    end
    
    if ~isempty(fstop)
        data.opts.bandstop = fstop;
        data = bead_filter_bandstop(data);
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
        data.opts.fpass = fpass;
        data.opts.cropTHPval = cropTHPval;
        %% Need to adapt to run on both dims and use subplots/etc
        data = bead_hp_allan_var(data, 'xCentresPx', ...
            fpass, cropTHPval, false);
        data = bead_hp_allan_var(data, 'yCentresPx', ...
            fpass, cropTHPval, false);
        fh = bead_plotHPData(data, setLims);
        
        if saveFigs
            saveas(fh, [data.fName '_HP.png'])
        end
        
        xStiff = calcStiffness(data.pro.xCentresHP);
        yStiff = calcStiffness(data.pro.yCentresHP);
        data.pro.stiffXYHP = [xStiff, yStiff];
        out(fileIdx).stiffHP = data.pro.stiffXYHP(1,:);

    end
    
    % Look at mean-square displacement (for cell-bead expts)
    if ~isempty(msdOffset)
%         data = bead_normMSD_polyfit(data, msdDim, msdOffset, msdNumT, doPlots, msdUseRaw, msdCentresRow, msdDoNorm);
        data = bead_normMSD(data, 'direction', msdDim, 'offset', msdOffset, 'numT', msdNumT, 'doPlots', doPlots, 'useRaw', msdUseRaw, 'centresRow', msdCentresRow, 'doNorm', msdDoNorm, 'errorBars', msdErrorBars);
        out(fileIdx).msdObj = data.pro.amsdObj;
        if saveFigs
            saveas(gcf, [data.fName '_MSD.png'])
        end
%         fh = gcf;
%         fh.Position = [404          42        1296         954];
%         for idx = 1:length(fh.Children)
%             ax = fh.Children(idx);
%             if isa(ax, 'matlab.graphics.axis.Axes') && isscalar(ax.Children)
%                 h = ax.Children;
%                 x = ax.Children.XData';
%                 y = ax.Children.YData';
%                 [dydx, tout] = msd_gradientor(x, y, 'lsq', 5);
%                 yyaxis(ax, 'right')
%                 semilogx(ax, tout, dydx, 'k--','LineWidth',2)
%                 hold(ax, 'on')
%                 plot(ax, xlim(ax), [1 1], ':','Color',0.75*[1 1 1], 'LineWidth', 3)
%                 ylim(ax, [0 2])
%                 
%             end
%         end
    end
    
    % Save if requested
    if saveData
        save(dataFile, 'data')
    end
    %pause(60)
    %close all
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
