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
    if exist([data.dirPath '/images_and_metadata/images_and_metadata_MMStack_Default.ome.tif'], 'file')
        data.Imstack  = bfopen([data.dirPath '/images_and_metadata/images_and_metadata_MMStack_Default.ome.tif']);
        metadata = fileread([data.dirPath '/images_and_metadata/images_and_metadata_MMStack_Default_metadata.txt']);
        data.metadata = jsondecode(metadata);
    elseif exist([data.dirPath '/ROI/ROI_MMStack_Default.ome.tif'], 'file')
        data.Imstack  = bfopen([data.dirPath '/ROI/ROI_MMStack_Default.ome.tif']);
    else
        warning('Could not find ROI images or metadata')
    end
    
    if exist([data.dirPath '/full_images_and_metadata/full_images_and_metadata_MMStack_Default.ome.tif'], 'file')
        data.ImstackFullFoV  = bfopen([data.dirPath '/full_images_and_metadata/full_images_and_metadata_MMStack_Default.ome.tif']);
    elseif exist([data.dirPath '/full/full_MMStack_Default.ome.tif'], 'file')
        data.ImstackFullFoV  = bfopen([data.dirPath '/full/full_MMStack_Default.ome.tif']);
    else
        warning('Could not find full FoV images')
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