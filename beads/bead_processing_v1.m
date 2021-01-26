%% Process multiple sets of bead data sequentially
% Experiment parameters
mPerPx = 0.07e-6;           % Camera pixel size calibration
laserPowers = 30:5:60;      % Laser power in % for the datasets used
ignoreDirs = {'hela_s3_w_bead_1','hela_s3_w_bead_2',...
    'hela_s3_w_bead_3','hela_s3_w_bead_z_3'}; % Directories to ignore

% Processing parameters
cropTs = {[1 10e5], [1 20e5], [1 20e5], [1 20e5], [1 10e5] };
fitPoly = [0 0 0 0 0 ]; % Fit a polynomial to remove drift. Only do this for calibration sets!
fitPolyOrder = 1;       % Order of polynomial to be fitted
calcStiff = 0;          % Calculate trap stiffness from position variance
fpass = 0;              % Pass frequency


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
data = struct([]);
out = {};

for fileIdx = 2%1:8%length(dirList)
    % Load all the data and the metadata
    data(1).dirPath = [dirList(fileIdx).folder '/' dirList(fileIdx).name];
    data.fName = dirList(fileIdx).name;
    data = bead_loadData(data);
    
    % Apply calibration and crop time
    data.opts.cropT = cropTs{fileIdx};
    data.opts.pOrder = fitPolyOrder*fitPoly(fileIdx);
    data.mPerPx = mPerPx;
    data = bead_preProcessCentres(data, mPerPx);
    %%
    
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
        fh = bead_plotRawData(data, setLims, fileIdx); %#ok<*UNRCH>
        fh.Name = data.fName;
    end
    
    % Open the Imstack file using ImageJ (kinda redundant since I have
    % compCentres)
    if showStack
        disp('MATLAB is locked until you close imagej')
        system(['imagej ' dirPath '/images_and_metadata/images_and_metadata_MMStack_Default.ome.tif']);
    end
    
    if fpass > 0
        data = bead_hp_allan_var(data, 'xCentresPx', ...
            fpass, doPlots, true);
        data = bead_hp_allan_var(data, 'yCentresPx', ...
            fpass, doPlots, true);
    end
    
    out{fileIdx} = data; %#ok<SAGROW>
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