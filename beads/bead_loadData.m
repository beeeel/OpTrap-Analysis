function data = bead_loadData(data, varargin)
%% data = loadBeadData(data/dirPath, [loadImages])
% Load centres, times, images and metadata for a given dataset

loadImages = true;
if nargin > 1
    loadImages = varargin{1};
end

if ~isstruct(data)
    data = struct('dirPath',data);
end

data.opts = bead_loadOpts(data);

if ~isfield(data.opts, 'skipSuffixes')
    data.opts.skipSuffixes = {};
else
    data.opts.skipSuffixes = strsplit(data.opts.skipSuffixes,',');
end

% New method: Find anything with a suffix and then record it and what
% suffix
if ~isfield(data, 'raw') 
    data.raw.suffixes = {};
elseif ~isfield(data.raw, 'suffixes')
    data.raw.suffixes = {};
end
errMsg = '';

if ~exist(data.dirPath,'dir')
    error('Folder %s does not exist!',data.dirPath)
end

dl = dir(data.dirPath);
dl = dl(endsWith({dl.name}, '.dat') & startsWith({dl.name}, {'X', 'Y','T','Times'}));
nx = sum(startsWith({dl.name}, 'X'));
ny = sum(startsWith({dl.name}, 'Y'));
nt = sum(startsWith({dl.name}, {'T','Times'}));

checkSizes;

tl = dl(startsWith({dl.name}, {'T','Times'}));
[~, tidx] = min([tl.bytes]);
data.raw.timeVecMs  = byteStreamToDouble([data.dirPath '/' tl(tidx).name]);
nP = length(data.raw.timeVecMs);

nNans = zeros(nx,1);
xdx = 1;
ydx = 1;
for idx = 1:length(dl)
    fName = sprintf('%s/%s', data.dirPath, dl(idx).name);
    suff = strsplit(dl(idx).name, {'X', 'Y','.dat'});
    suff = suff{end-1};
    if any(strcmp(data.opts.skipSuffixes, suff))
        warning('Skipping file %s because skipSuffixes',dl(idx).name)
    else
        % Screen for NaNs
        tmp = byteStreamToDouble(fName);
        if length(tmp) ~= nP
            warning('Skipping file %s with %i points (Times has %i points)', fName, length(tmp), nP)
        else
            nNans(xdx) = warnNaN(tmp, fName);
            tmp(isnan(tmp)) = mean(tmp(~isnan(tmp)));
            if startsWith(dl(idx).name, 'X')
                data.raw.xCentresPx(xdx,:) = tmp;
                if length(data.raw.suffixes)<xdx || strcmp(data.raw.suffixes{xdx}, suff)
                    data.raw.suffixes{xdx} = suff;
                else
                    error('Suffix mismatch like you said wouldn''t happen')
                end
                xdx = xdx + 1;
            elseif startsWith(dl(idx).name, 'Y')
                data.raw.yCentresPx(ydx,:) = tmp;
                if length(data.raw.suffixes)<ydx || strcmp(data.raw.suffixes{ydx}, suff)
                    data.raw.suffixes{ydx} = suff;
                else
                    error('Suffix mismatch like you said wouldn''t happen')
                end
                ydx = ydx + 1;
            end
        end
    end
end
data.raw.nNans = nNans;

if exist([data.dirPath '/subWidth.dat'],'file')
    data.raw.subWidth = byteStreamToDouble([data.dirPath '/subWidth.dat']);
end

if exist([data.dirPath '/I.dat'], 'file')
    data.raw.dcAvg = byteStreamToDouble([data.dirPath '/I.dat']);
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

function nNans = warnNaN(arr, name)
nNans = sum(isnan(arr));
if nNans
    warning(['Loaded and replaced ' num2str(nNans) ' from ' name])
end

end

function checkSizes

if nx ~= ny
    error('Found mismatched number of X (%i) and Y (%i) data', nx, ny)
elseif nt ~= 1 && mod(nx, nt) ~= 0
    error('Found mismatched number of T (%i) and XY (%i) data', nt, nx)
elseif all([nx ny nt] == 0)
    fprintf(errMsg);
    error('Could not find any centres data');
end

end

end
