function data = bead_loadData(data, varargin)
%% data = loadBeadData(data, [loadImages])
% Load centres, times, images and metadata for a given dataset

loadImages = true;
if nargin > 1
    loadImages = varargin{1};
end

% New method: Find anything with a suffix and then record it and what
% suffix
if ~isfield(data, 'raw') 
    data.raw.suffixes = {};
elseif ~isfield(data.raw, 'suffixes')
    data.raw.suffixes = {};
end
idx = 1;
errMsg = '';

if ~exist(data.dirPath,'dir')
    error('Folder %s does not exist!',data.dirPath)
end

for suff = {'Centres', '', 'simple', 'quart', 'gauss', 'l', 'r', 'b', 's', 'v', 'lq', 'rq'}
    fNameX = [data.dirPath '/X' suff{:} '.dat'];
    fNameY = [data.dirPath '/Y' suff{:} '.dat'];
    if exist(fNameX, 'file') && exist(fNameY, 'file')
        % Screen for NaNs
        tmp = byteStreamToDouble(fNameX);
        warnNaN(tmp, fNameX);
        tmp(isnan(tmp)) = mean(tmp(~isnan(tmp)));
        data.raw.xCentresPx(idx,:) = tmp;
        
        tmp = byteStreamToDouble(fNameY);
        warnNaN(tmp, fNameY);
        tmp(isnan(tmp)) = mean(tmp(~isnan(tmp)));
        data.raw.yCentresPx(idx,:) = tmp;
        
        data.raw.suffixes = [data.raw.suffixes, suff];
        idx = idx + 1;
    else
        errMsg = [errMsg 'File ' fNameX ' does not exist\n'];
    end
end

if idx == 1
    fprintf(errMsg);
    error('Could not find any centres data');
end

if exist([data.dirPath '/I.dat'], 'file')
    data.raw.dcAvg = byteStreamToDouble([data.dirPath '/I.dat']);
end

if exist([data.dirPath '/Times.dat'], 'file')
    data.raw.timeVecMs  = byteStreamToDouble([data.dirPath '/Times.dat']);
elseif exist([data.dirPath '/T.dat'], 'file')
    data.raw.timeVecMs  = byteStreamToDouble([data.dirPath '/T.dat']);
else
    disp(['File ' data.dirPath '/Times.dat does not exist'])
    disp(['File ' data.dirPath '/T.dat does not exist'])
    error('Could not find any times data');
end

if loadImages
    % More sensible method: If the folder exists, find all .tif files
    % within it. If there's only one, load it. Likewise for metadata.txt.
    if exist([data.dirPath '/images_and_metadata/'], 'dir')
        dirList = dir([data.dirPath '/images_and_metadata/']);
        imList = dirList(endsWith({dirList.name}, '.tif'));
        if isscalar(imList)
            fPath = strjoin({imList.folder imList.name}, '/');
            data.Imstack = bfopen(fPath);
        else
            error('Found multiple (or 0) TIFs in ROI folder')
        end
        txtList = dirList(endsWith({dirList.name}, 'metadata.txt'));
        
        if isscalar(txtList)
            fPath = strjoin({txtList.folder txtList.name}, '/');
            metadata = fileread(fPath);
            data.metadata = jsondecode(metadata);
            
            % Get ROI position
            roi = str2double( strsplit( data.metadata.FrameKey_0_0_0.ROI, '-' ) );
            data.opts.roi = roi;
        else
            error('Found multiple (or 0) files ending with metadata.txt in ROI folder')
        end
    else
        warning('Could not find ROI images or metadata')
    end
    
    if exist([data.dirPath '/full_images_and_metadata/'], 'dir')
        dirList = dir([data.dirPath '/full_images_and_metadata/']);
        dirList = dirList(endsWith({dirList.name}, 'tif'));
        if isscalar(dirList)
            fPath = strjoin({dirList.folder dirList.name}, '/');
            data.ImstackFullFoV  = bfopen(fPath);
        else
            error('Found multiple (or 0) TIFs in full FoV folder')
        end
    else
        warning('Could not find directory for full FoV images')
    end
else
    warning('Instructed to not load images')
end

data.nPoints = length(data.raw.xCentresPx);

function warnNaN(arr, name)
nNans = sum(isnan(arr));
if nNans
    warning(['Loaded and replaced ' num2str(nNans) ' from ' name])
end